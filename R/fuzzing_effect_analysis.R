#!/usr/bin/env Rscript
# R/fuzzing_effect_analysis.R
# Visualize and quantify how FIA coordinate fuzzing affects spatial predictions
# Author: Soren Donisvitch
# Date: December 2024

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sf)
  library(terra)
  library(xgboost)
  library(patchwork)
  library(viridis)
  library(fs)
})

# =============================================================================
# MAIN FUNCTION
# =============================================================================

analyze_fuzzing_effects <- function(
    nefin_path = "data/processed/nefin_processed.csv",
    nefin_climate_path = "data/processed/climate_at_plots/nefin_climate.csv",
    nefin_ndvi_path = "data/processed/ndvi_at_plots/nefin_ndvi.csv",
    ndvi_raster_path = "data/processed/ndvi/modis/MODIS_NDVI_5yr_blocked_2020_2024.tif",
    tmean_raster_path = "data/processed/prism/prism_tmean_ne_2020_2024.tif",
    ppt_raster_path = "data/processed/prism/prism_ppt_ne_2020_2024.tif",
    states_path = "data/boundaries/states_5070.geojson",
    model_dir = "runs/spatial_model_comparison",
    output_dir = "runs/fuzzing_effect_analysis",
    n_fuzz_replicates = 100,
    fuzz_radius_m = 1609  # 1 mile in meters (FIA fuzzing radius)
) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════╗\n")
  cat("║  Fuzzing Effect Analysis                                              ║\n")
  cat("║  How does coordinate uncertainty affect spatial predictions?          ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  fs::dir_create(output_dir, recurse = TRUE)
  
  # =========================================================================
  # 1. LOAD DATA AND MODELS
  # =========================================================================
  
  cat("═══════════════════════════════════════════════════════════════════════\n")
  cat("STEP 1: Loading Data\n")
  cat("═══════════════════════════════════════════════════════════════════════\n\n")
  
  # Load NEFIN data with covariates
  nefin <- read_csv(nefin_path, show_col_types = FALSE)
  cat("  NEFIN plots:", nrow(nefin), "\n")
  
  # Load covariates
  if (fs::file_exists(nefin_climate_path)) {
    nefin_climate <- read_csv(nefin_climate_path, show_col_types = FALSE) %>%
      distinct(CN, .keep_all = TRUE)
    nefin <- nefin %>% 
      left_join(nefin_climate %>% select(CN, tmean, ppt), by = "CN")
  }
  
  if (fs::file_exists(nefin_ndvi_path)) {
    nefin_ndvi <- read_csv(nefin_ndvi_path, show_col_types = FALSE) %>%
      distinct(CN, .keep_all = TRUE)
    nefin <- nefin %>% 
      left_join(nefin_ndvi %>% select(CN, ndvi_modis), by = "CN")
  }
  
  # Filter to complete cases
  nefin_complete <- nefin %>%
    filter(!is.na(aglb_Mg_per_ha), !is.na(ndvi_modis), !is.na(tmean), !is.na(ppt))
  
  cat("  NEFIN with complete covariates:", nrow(nefin_complete), "\n")
  
  # Find coordinate columns
  lon_col <- names(nefin)[grepl("lon", names(nefin), ignore.case = TRUE)][1]
  lat_col <- names(nefin)[grepl("lat", names(nefin), ignore.case = TRUE)][1]
  
  # Load rasters for covariate extraction
  cat("\n  Loading rasters...\n")
  
  r_ndvi <- if (fs::file_exists(ndvi_raster_path)) {
    terra::rast(ndvi_raster_path)
  } else NULL
  
  r_tmean <- if (fs::file_exists(tmean_raster_path)) {
    terra::rast(tmean_raster_path)
  } else NULL
  
  r_ppt <- if (fs::file_exists(ppt_raster_path)) {
    terra::rast(ppt_raster_path)
  } else NULL
  
  if (is.null(r_ndvi) || is.null(r_tmean) || is.null(r_ppt)) {
    cat("  ⚠ Missing raster files - some analyses will be skipped\n")
  } else {
    cat("  ✓ All rasters loaded\n")
  }
  
  # Load state boundaries
  states <- if (fs::file_exists(states_path)) {
    st_read(states_path, quiet = TRUE)
  } else NULL
  
  # =========================================================================
  # 2. MONTE CARLO FUZZING SENSITIVITY
  # =========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════\n")
  cat("STEP 2: Monte Carlo Fuzzing Sensitivity Analysis\n")
  cat("═══════════════════════════════════════════════════════════════════════\n\n")
  
  cat("  Simulating FIA-style fuzzing (±", fuzz_radius_m, "m) for", n_fuzz_replicates, "replicates\n\n")
  
  # Sample subset of NEFIN plots for fuzzing analysis
  set.seed(42)
  n_sample <- min(500, nrow(nefin_complete))
  sample_idx <- sample(1:nrow(nefin_complete), n_sample)
  nefin_sample <- nefin_complete[sample_idx, ]
  
  cat("  Using", n_sample, "plots for fuzzing sensitivity\n")
  
  # Convert to sf
  nefin_sf <- st_as_sf(nefin_sample, 
                       coords = c(lon_col, lat_col), 
                       crs = 4326) %>%
    st_transform(5070)  # Albers Equal Area for distance-based fuzzing
  
  # Storage for fuzzed extractions
  fuzz_results <- list()
  
  if (!is.null(r_ndvi) && !is.null(r_tmean) && !is.null(r_ppt)) {
    
    cat("  Extracting covariates at fuzzed locations...\n")
    pb <- txtProgressBar(min = 0, max = n_fuzz_replicates, style = 3)
    
    for (rep in 1:n_fuzz_replicates) {
      # Generate random fuzz offsets
      angles <- runif(n_sample, 0, 2 * pi)
      distances <- runif(n_sample, 0, fuzz_radius_m)
      
      dx <- distances * cos(angles)
      dy <- distances * sin(angles)
      
      # Apply fuzzing
      coords_orig <- st_coordinates(nefin_sf)
      coords_fuzz <- cbind(coords_orig[, 1] + dx, coords_orig[, 2] + dy)
      
      nefin_fuzz <- st_as_sf(data.frame(id = 1:n_sample, x = coords_fuzz[, 1], y = coords_fuzz[, 2]),
                             coords = c("x", "y"), crs = 5070)
      
      # Transform to raster CRS and extract
      nefin_fuzz_wgs <- st_transform(nefin_fuzz, crs(r_ndvi))
      
      ndvi_vals <- terra::extract(r_ndvi, terra::vect(nefin_fuzz_wgs))[[2]]
      
      nefin_fuzz_prism <- st_transform(nefin_fuzz, crs(r_tmean))
      tmean_vals <- terra::extract(r_tmean, terra::vect(nefin_fuzz_prism))[[2]]
      ppt_vals <- terra::extract(r_ppt, terra::vect(nefin_fuzz_prism))[[2]]
      
      fuzz_results[[rep]] <- data.frame(
        plot_idx = 1:n_sample,
        replicate = rep,
        ndvi_fuzz = ndvi_vals,
        tmean_fuzz = tmean_vals,
        ppt_fuzz = ppt_vals,
        fuzz_distance = distances
      )
      
      setTxtProgressBar(pb, rep)
    }
    close(pb)
    
    # Combine results
    fuzz_df <- bind_rows(fuzz_results)
    
    # Add true values
    fuzz_df <- fuzz_df %>%
      left_join(
        data.frame(
          plot_idx = 1:n_sample,
          ndvi_true = nefin_sample$ndvi_modis,
          tmean_true = nefin_sample$tmean,
          ppt_true = nefin_sample$ppt,
          biomass_true = nefin_sample$aglb_Mg_per_ha
        ),
        by = "plot_idx"
      )
    
    # Calculate covariate errors due to fuzzing
    fuzz_df <- fuzz_df %>%
      mutate(
        ndvi_error = ndvi_fuzz - ndvi_true,
        tmean_error = tmean_fuzz - tmean_true,
        ppt_error = ppt_fuzz - ppt_true
      )
    
    # =========================================================================
    # 3. QUANTIFY COVARIATE UNCERTAINTY FROM FUZZING
    # =========================================================================
    
    cat("\n═══════════════════════════════════════════════════════════════════════\n")
    cat("STEP 3: Quantifying Covariate Uncertainty\n")
    cat("═══════════════════════════════════════════════════════════════════════\n\n")
    
    # Summary statistics per plot
    plot_uncertainty <- fuzz_df %>%
      group_by(plot_idx) %>%
      summarise(
        ndvi_true = first(ndvi_true),
        ndvi_fuzz_mean = mean(ndvi_fuzz, na.rm = TRUE),
        ndvi_fuzz_sd = sd(ndvi_fuzz, na.rm = TRUE),
        ndvi_fuzz_range = max(ndvi_fuzz, na.rm = TRUE) - min(ndvi_fuzz, na.rm = TRUE),
        tmean_true = first(tmean_true),
        tmean_fuzz_sd = sd(tmean_fuzz, na.rm = TRUE),
        ppt_true = first(ppt_true),
        ppt_fuzz_sd = sd(ppt_fuzz, na.rm = TRUE),
        biomass_true = first(biomass_true),
        .groups = "drop"
      )
    
    # Overall summary
    cat("  Covariate uncertainty due to ±1.6km fuzzing:\n\n")
    
    cat("  NDVI:\n")
    cat("    Mean absolute error:", round(mean(abs(fuzz_df$ndvi_error), na.rm = TRUE), 4), "\n")
    cat("    SD of fuzzed values:", round(mean(plot_uncertainty$ndvi_fuzz_sd, na.rm = TRUE), 4), "\n")
    cat("    Mean range across replicates:", round(mean(plot_uncertainty$ndvi_fuzz_range, na.rm = TRUE), 4), "\n")
    
    cat("\n  Temperature (°C):\n")
    cat("    Mean SD across replicates:", round(mean(plot_uncertainty$tmean_fuzz_sd, na.rm = TRUE), 2), "°C\n")
    
    cat("\n  Precipitation (mm):\n")
    cat("    Mean SD across replicates:", round(mean(plot_uncertainty$ppt_fuzz_sd, na.rm = TRUE), 1), "mm\n")
    
    # Save uncertainty summary
    uncertainty_summary <- data.frame(
      variable = c("NDVI", "Tmean", "PPT"),
      mean_abs_error = c(
        mean(abs(fuzz_df$ndvi_error), na.rm = TRUE),
        mean(abs(fuzz_df$tmean_error), na.rm = TRUE),
        mean(abs(fuzz_df$ppt_error), na.rm = TRUE)
      ),
      mean_sd = c(
        mean(plot_uncertainty$ndvi_fuzz_sd, na.rm = TRUE),
        mean(plot_uncertainty$tmean_fuzz_sd, na.rm = TRUE),
        mean(plot_uncertainty$ppt_fuzz_sd, na.rm = TRUE)
      ),
      pct_of_mean = c(
        100 * mean(plot_uncertainty$ndvi_fuzz_sd, na.rm = TRUE) / mean(plot_uncertainty$ndvi_true, na.rm = TRUE),
        100 * mean(plot_uncertainty$tmean_fuzz_sd, na.rm = TRUE) / mean(abs(plot_uncertainty$tmean_true), na.rm = TRUE),
        100 * mean(plot_uncertainty$ppt_fuzz_sd, na.rm = TRUE) / mean(plot_uncertainty$ppt_true, na.rm = TRUE)
      )
    )
    
    cat("\n  Uncertainty as % of mean:\n")
    print(uncertainty_summary, row.names = FALSE)
    
    write_csv(uncertainty_summary, fs::path(output_dir, "covariate_uncertainty_summary.csv"))
    
    # =========================================================================
    # 4. PREDICTION UNCERTAINTY FROM FUZZING
    # =========================================================================
    
    cat("\n═══════════════════════════════════════════════════════════════════════\n")
    cat("STEP 4: Prediction Uncertainty from Fuzzing\n")
    cat("═══════════════════════════════════════════════════════════════════════\n\n")
    
    # Train a simple model on NEFIN (true coords)
    nefin_train_data <- nefin_complete %>%
      filter(!is.na(ndvi_modis), !is.na(tmean), !is.na(ppt))
    
    X_train <- as.matrix(nefin_train_data %>% select(ndvi_modis, tmean, ppt))
    y_train <- nefin_train_data$aglb_Mg_per_ha
    
    dtrain <- xgb.DMatrix(data = X_train, label = y_train)
    
    cat("  Training XGBoost model on NEFIN true coordinates...\n")
    
    xgb_model <- xgb.train(
      data = dtrain,
      nrounds = 50,
      params = list(
        objective = "reg:squarederror",
        eta = 0.1,
        max_depth = 4
      ),
      verbose = 0
    )
    
    # Predict at true locations
    X_sample <- as.matrix(nefin_sample %>% select(ndvi_modis, tmean, ppt))
    pred_true <- predict(xgb_model, xgb.DMatrix(X_sample))
    
    # Predict at all fuzzed locations
    fuzz_df$pred_fuzz <- NA
    
    for (rep in 1:n_fuzz_replicates) {
      rep_data <- fuzz_df %>% filter(replicate == rep)
      X_fuzz <- as.matrix(rep_data %>% select(ndvi_fuzz, tmean_fuzz, ppt_fuzz))
      
      # Handle NAs
      valid_idx <- complete.cases(X_fuzz)
      if (sum(valid_idx) > 0) {
        preds <- rep(NA, nrow(X_fuzz))
        preds[valid_idx] <- predict(xgb_model, xgb.DMatrix(X_fuzz[valid_idx, , drop = FALSE]))
        fuzz_df$pred_fuzz[fuzz_df$replicate == rep] <- preds
      }
    }
    
    # Summarize prediction uncertainty per plot
    pred_uncertainty <- fuzz_df %>%
      group_by(plot_idx) %>%
      summarise(
        biomass_true = first(biomass_true),
        pred_true = pred_true[plot_idx[1]],
        pred_fuzz_mean = mean(pred_fuzz, na.rm = TRUE),
        pred_fuzz_sd = sd(pred_fuzz, na.rm = TRUE),
        pred_fuzz_range = max(pred_fuzz, na.rm = TRUE) - min(pred_fuzz, na.rm = TRUE),
        pred_fuzz_min = min(pred_fuzz, na.rm = TRUE),
        pred_fuzz_max = max(pred_fuzz, na.rm = TRUE),
        .groups = "drop"
      )
    
    cat("\n  Prediction uncertainty due to fuzzing:\n")
    cat("    Mean prediction SD:", round(mean(pred_uncertainty$pred_fuzz_sd, na.rm = TRUE), 2), "Mg/ha\n")
    cat("    Mean prediction range:", round(mean(pred_uncertainty$pred_fuzz_range, na.rm = TRUE), 2), "Mg/ha\n")
    cat("    Max prediction range:", round(max(pred_uncertainty$pred_fuzz_range, na.rm = TRUE), 2), "Mg/ha\n")
    cat("    Mean biomass:", round(mean(pred_uncertainty$biomass_true, na.rm = TRUE), 2), "Mg/ha\n")
    cat("    Prediction SD as % of mean:", 
        round(100 * mean(pred_uncertainty$pred_fuzz_sd, na.rm = TRUE) / mean(pred_uncertainty$biomass_true, na.rm = TRUE), 1), "%\n")
    
    write_csv(pred_uncertainty, fs::path(output_dir, "prediction_uncertainty_by_plot.csv"))
    
    # =========================================================================
    # 5. CREATE FIGURES
    # =========================================================================
    
    cat("\n═══════════════════════════════════════════════════════════════════════\n")
    cat("STEP 5: Creating Figures\n")
    cat("═══════════════════════════════════════════════════════════════════════\n\n")
    
    # --- Figure 1: Covariate error distributions ---
    fig1a <- ggplot(fuzz_df, aes(x = ndvi_error)) +
      geom_histogram(bins = 50, fill = "#3498db", alpha = 0.7) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      labs(title = "NDVI Error", x = "NDVI (fuzzed - true)", y = "Count") +
      theme_minimal()
    
    fig1b <- ggplot(fuzz_df, aes(x = tmean_error)) +
      geom_histogram(bins = 50, fill = "#27ae60", alpha = 0.7) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      labs(title = "Temperature Error", x = "Tmean (fuzzed - true) °C", y = "Count") +
      theme_minimal()
    
    fig1c <- ggplot(fuzz_df, aes(x = ppt_error)) +
      geom_histogram(bins = 50, fill = "#9b59b6", alpha = 0.7) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      labs(title = "Precipitation Error", x = "PPT (fuzzed - true) mm", y = "Count") +
      theme_minimal()
    
    fig1 <- fig1a + fig1b + fig1c +
      plot_annotation(
        title = "Covariate Extraction Errors Due to ±1.6km Coordinate Fuzzing",
        subtitle = paste0("Based on ", n_fuzz_replicates, " Monte Carlo replicates of ", n_sample, " NEFIN plots"),
        theme = theme(plot.title = element_text(face = "bold", size = 14))
      )
    
    ggsave(fs::path(output_dir, "fig1_covariate_errors.png"), fig1, 
           width = 14, height = 5, dpi = 300)
    cat("  ✓ Saved fig1_covariate_errors.png\n")
    
    # --- Figure 2: Prediction uncertainty distribution ---
    fig2a <- ggplot(pred_uncertainty, aes(x = pred_fuzz_sd)) +
      geom_histogram(bins = 30, fill = "#e74c3c", alpha = 0.7) +
      geom_vline(xintercept = mean(pred_uncertainty$pred_fuzz_sd, na.rm = TRUE), 
                 linetype = "dashed", color = "black", linewidth = 1) +
      labs(
        title = "A) Distribution of Prediction SD",
        x = "Prediction SD (Mg/ha)", y = "Number of Plots"
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))
    
    fig2b <- ggplot(pred_uncertainty, aes(x = pred_fuzz_range)) +
      geom_histogram(bins = 30, fill = "#e67e22", alpha = 0.7) +
      labs(
        title = "B) Distribution of Prediction Range",
        subtitle = "Max - Min across 100 fuzz replicates",
        x = "Prediction Range (Mg/ha)", y = "Number of Plots"
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))
    
    fig2 <- fig2a + fig2b +
      plot_annotation(
        title = "Biomass Prediction Uncertainty Due to Coordinate Fuzzing",
        theme = theme(plot.title = element_text(face = "bold", size = 14))
      )
    
    ggsave(fs::path(output_dir, "fig2_prediction_uncertainty.png"), fig2, 
           width = 12, height = 5, dpi = 300)
    cat("  ✓ Saved fig2_prediction_uncertainty.png\n")
    
    # --- Figure 3: Prediction envelope for sample plots ---
    # Show a few example plots with prediction ranges
    example_plots <- pred_uncertainty %>%
      arrange(desc(pred_fuzz_range)) %>%
      head(20) %>%
      mutate(plot_rank = row_number())
    
    fig3 <- ggplot(example_plots, aes(x = reorder(factor(plot_idx), -pred_fuzz_range))) +
      geom_errorbar(aes(ymin = pred_fuzz_min, ymax = pred_fuzz_max), width = 0.3, color = "#3498db") +
      geom_point(aes(y = pred_true), color = "#27ae60", size = 3) +
      geom_point(aes(y = biomass_true), color = "#e74c3c", size = 3, shape = 17) +
      labs(
        title = "Prediction Range from Coordinate Fuzzing (Top 20 Most Variable Plots)",
        subtitle = "Blue bars = range across 100 fuzz replicates | Green = true coord prediction | Red = observed biomass",
        x = "Plot", y = "Biomass (Mg/ha)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        axis.text.x = element_blank()
      )
    
    ggsave(fs::path(output_dir, "fig3_prediction_envelopes.png"), fig3, 
           width = 12, height = 6, dpi = 300)
    cat("  ✓ Saved fig3_prediction_envelopes.png\n")
    
    # --- Figure 4: NDVI uncertainty vs prediction uncertainty ---
    plot_uncertainty <- plot_uncertainty %>%
      left_join(pred_uncertainty %>% select(plot_idx, pred_fuzz_sd, pred_fuzz_range), by = "plot_idx")
    
    fig4 <- ggplot(plot_uncertainty, aes(x = ndvi_fuzz_sd, y = pred_fuzz_sd)) +
      geom_point(alpha = 0.5, color = "#3498db") +
      geom_smooth(method = "lm", color = "#e74c3c") +
      labs(
        title = "Does NDVI Uncertainty Drive Prediction Uncertainty?",
        subtitle = paste0("Correlation: r = ", 
                          round(cor(plot_uncertainty$ndvi_fuzz_sd, plot_uncertainty$pred_fuzz_sd, use = "complete.obs"), 3)),
        x = "NDVI SD (across fuzz replicates)",
        y = "Prediction SD (Mg/ha)"
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(fs::path(output_dir, "fig4_ndvi_vs_pred_uncertainty.png"), fig4, 
           width = 8, height = 6, dpi = 300)
    cat("  ✓ Saved fig4_ndvi_vs_pred_uncertainty.png\n")
    
    # --- Figure 5: Relationship attenuation ---
    # Compare NDVI-biomass relationship at true vs fuzzed coords
    
    # True relationship
    true_fit <- lm(biomass_true ~ ndvi_true, data = plot_uncertainty)
    true_slope <- coef(true_fit)[2]
    true_r2 <- summary(true_fit)$r.squared
    
    # Fuzzed relationship (using mean of fuzzed NDVI)
    fuzz_fit <- lm(biomass_true ~ ndvi_fuzz_mean, data = plot_uncertainty)
    fuzz_slope <- coef(fuzz_fit)[2]
    fuzz_r2 <- summary(fuzz_fit)$r.squared
    
    attenuation_pct <- 100 * (true_slope - fuzz_slope) / true_slope
    
    cat("\n  Relationship Attenuation:\n")
    cat("    True NDVI slope:", round(true_slope, 2), "(R² =", round(true_r2, 4), ")\n")
    cat("    Fuzzed NDVI slope:", round(fuzz_slope, 2), "(R² =", round(fuzz_r2, 4), ")\n")
    cat("    Attenuation:", round(attenuation_pct, 1), "%\n")
    
    scatter_data <- plot_uncertainty %>%
      select(plot_idx, biomass_true, ndvi_true, ndvi_fuzz_mean) %>%
      pivot_longer(cols = c(ndvi_true, ndvi_fuzz_mean), 
                   names_to = "type", values_to = "ndvi") %>%
      mutate(type = ifelse(type == "ndvi_true", 
                           paste0("True Coords (slope=", round(true_slope, 1), ")"),
                           paste0("Fuzzed Mean (slope=", round(fuzz_slope, 1), ")")))
    
    fig5 <- ggplot(scatter_data, aes(x = ndvi, y = biomass_true, color = type)) +
      geom_point(alpha = 0.4) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
      scale_color_manual(values = c("#27ae60", "#e74c3c")) +
      labs(
        title = "NDVI-Biomass Relationship: True vs Fuzzed Coordinates",
        subtitle = paste0("Fuzzing attenuates the slope by ", round(attenuation_pct, 1), "%"),
        x = "NDVI", y = "Biomass (Mg/ha)",
        color = "Coordinate Type"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "bottom"
      )
    
    ggsave(fs::path(output_dir, "fig5_relationship_attenuation.png"), fig5, 
           width = 10, height = 7, dpi = 300)
    cat("  ✓ Saved fig5_relationship_attenuation.png\n")
    
    # Save attenuation summary
    attenuation_df <- data.frame(
      metric = c("True slope", "Fuzzed slope", "Attenuation %", "True R²", "Fuzzed R²"),
      value = c(true_slope, fuzz_slope, attenuation_pct, true_r2, fuzz_r2)
    )
    write_csv(attenuation_df, fs::path(output_dir, "relationship_attenuation.csv"))
    
  } else {
    cat("  ⚠ Skipping Monte Carlo analysis - missing raster data\n")
  }
  
  # =========================================================================
  # 6. SUMMARY
  # =========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════╗\n")
  cat("║  ANALYSIS COMPLETE                                                    ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  cat("Output directory:", output_dir, "\n\n")
  
  cat("Key findings:\n")
  if (exists("uncertainty_summary")) {
    cat("  • NDVI uncertainty from fuzzing:", round(uncertainty_summary$mean_sd[1], 4), 
        "(", round(uncertainty_summary$pct_of_mean[1], 1), "% of mean)\n")
  }
  if (exists("pred_uncertainty")) {
    cat("  • Mean prediction SD from fuzzing:", round(mean(pred_uncertainty$pred_fuzz_sd, na.rm = TRUE), 2), "Mg/ha\n")
    cat("  • Max prediction range:", round(max(pred_uncertainty$pred_fuzz_range, na.rm = TRUE), 2), "Mg/ha\n")
  }
  if (exists("attenuation_pct")) {
    cat("  • NDVI-biomass slope attenuation:", round(attenuation_pct, 1), "%\n")
  }
  
  cat("\nFigures:\n")
  cat("  • fig1_covariate_errors.png - Distribution of covariate extraction errors\n")
  cat("  • fig2_prediction_uncertainty.png - Prediction uncertainty from fuzzing\n")
  cat("  • fig3_prediction_envelopes.png - Example plots with prediction ranges\n")
  cat("  • fig4_ndvi_vs_pred_uncertainty.png - NDVI uncertainty drives prediction uncertainty\n")
  cat("  • fig5_relationship_attenuation.png - Slope attenuation from fuzzing\n")
  
  invisible(list(
    uncertainty_summary = if (exists("uncertainty_summary")) uncertainty_summary else NULL,
    pred_uncertainty = if (exists("pred_uncertainty")) pred_uncertainty else NULL,
    attenuation_pct = if (exists("attenuation_pct")) attenuation_pct else NULL
  ))
}

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  analyze_fuzzing_effects()
}
