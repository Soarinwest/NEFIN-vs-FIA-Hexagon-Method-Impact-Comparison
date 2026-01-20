#!/usr/bin/env Rscript
# R/landscape_heterogeneity_analysis.R
# Does coordinate fuzzing hurt more in heterogeneous landscapes?
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
# MAIN ANALYSIS
# =============================================================================

analyze_heterogeneity_effect <- function(
    nefin_path = "data/processed/nefin_processed.csv",
    ndvi_path = "data/processed/ndvi/modis/MODIS_NDVI_5yr_blocked_2020_2024.tif",
    tmean_path = "data/processed/prism/prism_tmean_ne_2020_2024.tif",
    ppt_path = "data/processed/prism/prism_ppt_ne_2020_2024.tif",
    states_path = "data/boundaries/states_5070.geojson",
    output_dir = "runs/heterogeneity_analysis",
    fuzz_radius_m = 1609,
    n_replicates = 50
) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  Landscape Heterogeneity Analysis                                         ║\n")
  cat("║  Does fuzzing hurt more in variable landscapes?                           ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  fs::dir_create(output_dir, recurse = TRUE)
  
  # ===========================================================================
  # 1. LOAD DATA
  # ===========================================================================
  
  cat("═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 1: Loading Data\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Load NEFIN
  nefin <- read_csv(nefin_path, show_col_types = FALSE)
  lon_col <- names(nefin)[grepl("lon", names(nefin), ignore.case = TRUE)][1]
  lat_col <- names(nefin)[grepl("lat", names(nefin), ignore.case = TRUE)][1]
  
  cat("  NEFIN plots:", nrow(nefin), "\n")
  
  # Load NDVI raster
  r_ndvi <- rast(ndvi_path)
  r_tmean <- rast(tmean_path)
  r_ppt <- rast(ppt_path)
  
  cat("  NDVI resolution:", res(r_ndvi)[1], "m\n")
  
  # ===========================================================================
  # 2. COMPUTE LOCAL HETEROGENEITY
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 2: Computing Local NDVI Heterogeneity\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Convert plots to sf
  nefin_sf <- st_as_sf(nefin, coords = c(lon_col, lat_col), crs = 4326) %>%
    st_transform(crs(r_ndvi))
  
  # For each plot, compute NDVI statistics in 1.6km radius (fuzzing zone)
  cat("  Computing NDVI statistics in 1.6km radius around each plot...\n")
  
  # Create circular buffer template
  buffer_radius <- fuzz_radius_m
  
  # Sample for speed
  set.seed(42)
  n_sample <- min(1000, nrow(nefin))
  sample_idx <- sample(1:nrow(nefin_sf), n_sample)
  nefin_sample <- nefin_sf[sample_idx, ]
  
  cat("  Using", n_sample, "plots\n")
  
  # Compute heterogeneity metrics
  pb <- txtProgressBar(min = 0, max = n_sample, style = 3)
  
  heterogeneity <- list()
  
  for (i in 1:n_sample) {
    # Create buffer
    pt <- nefin_sample[i, ]
    buffer <- st_buffer(pt, buffer_radius)
    
    # Extract NDVI within buffer
    ndvi_vals <- terra::extract(r_ndvi, vect(buffer), fun = NULL)[[2]]
    ndvi_vals <- ndvi_vals[!is.na(ndvi_vals)]
    
    if (length(ndvi_vals) > 5) {
      heterogeneity[[i]] <- data.frame(
        idx = i,
        CN = nefin_sample$CN[i],
        ndvi_mean = mean(ndvi_vals),
        ndvi_sd = sd(ndvi_vals),
        ndvi_cv = sd(ndvi_vals) / mean(ndvi_vals),  # Coefficient of variation
        ndvi_range = max(ndvi_vals) - min(ndvi_vals),
        ndvi_iqr = IQR(ndvi_vals),
        n_pixels = length(ndvi_vals)
      )
    }
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  het_df <- bind_rows(heterogeneity)
  cat("\n  Computed heterogeneity for", nrow(het_df), "plots\n")
  
  # Add biomass
  het_df <- het_df %>%
    left_join(nefin %>% select(CN, aglb_Mg_per_ha), by = "CN")
  
  # Summary
  cat("\n  NDVI Heterogeneity Summary (1.6km radius):\n")
  cat("    Mean SD:", round(mean(het_df$ndvi_sd, na.rm = TRUE), 4), "\n")
  cat("    Mean CV:", round(mean(het_df$ndvi_cv, na.rm = TRUE), 3), "\n")
  cat("    Mean Range:", round(mean(het_df$ndvi_range, na.rm = TRUE), 3), "\n")
  
  # ===========================================================================
  # 3. MONTE CARLO FUZZING BY HETEROGENEITY CLASS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 3: Fuzzing Analysis by Heterogeneity Class\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Classify plots by heterogeneity
  het_df <- het_df %>%
    mutate(
      het_class = cut(ndvi_sd, 
                      breaks = quantile(ndvi_sd, probs = c(0, 0.33, 0.67, 1), na.rm = TRUE),
                      labels = c("Low", "Medium", "High"),
                      include.lowest = TRUE)
    )
  
  cat("  Heterogeneity classes:\n")
  print(table(het_df$het_class))
  
  # Get coordinates for fuzzing
  het_sf <- nefin_sample[het_df$idx, ] %>%
    st_transform(5070)
  
  coords <- st_coordinates(het_sf)
  
  # Extract true NDVI and climate at original locations
  het_sf_ndvi <- st_transform(het_sf, crs(r_ndvi))
  het_df$ndvi_true <- terra::extract(r_ndvi, vect(het_sf_ndvi))[[2]]
  
  het_sf_prism <- st_transform(het_sf, crs(r_tmean))
  het_df$tmean <- terra::extract(r_tmean, vect(het_sf_prism))[[2]]
  het_df$ppt <- terra::extract(r_ppt, vect(het_sf_prism))[[2]]
  
  # Monte Carlo fuzzing
  cat("\n  Running", n_replicates, "fuzzing replicates...\n")
  
  fuzz_results <- list()
  pb <- txtProgressBar(min = 0, max = n_replicates, style = 3)
  
  for (rep in 1:n_replicates) {
    # Random offsets
    angles <- runif(nrow(het_df), 0, 2 * pi)
    distances <- runif(nrow(het_df), 0, fuzz_radius_m)
    
    dx <- distances * cos(angles)
    dy <- distances * sin(angles)
    
    coords_fuzz <- cbind(coords[, 1] + dx, coords[, 2] + dy)
    
    fuzz_sf <- st_as_sf(
      data.frame(x = coords_fuzz[, 1], y = coords_fuzz[, 2]),
      coords = c("x", "y"), crs = 5070
    ) %>% st_transform(crs(r_ndvi))
    
    ndvi_fuzz <- terra::extract(r_ndvi, vect(fuzz_sf))[[2]]
    
    fuzz_results[[rep]] <- data.frame(
      idx = 1:nrow(het_df),
      replicate = rep,
      ndvi_fuzz = ndvi_fuzz,
      ndvi_error = ndvi_fuzz - het_df$ndvi_true
    )
    
    setTxtProgressBar(pb, rep)
  }
  close(pb)
  
  fuzz_df <- bind_rows(fuzz_results)
  
  # Add heterogeneity class
  fuzz_df <- fuzz_df %>%
    left_join(het_df %>% select(idx, het_class, ndvi_sd, aglb_Mg_per_ha, tmean, ppt, ndvi_true),
              by = "idx")
  
  # ===========================================================================
  # 4. ANALYZE ERROR BY HETEROGENEITY
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 4: Analyzing Error by Heterogeneity Class\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Error statistics by class
  error_by_class <- fuzz_df %>%
    group_by(het_class) %>%
    summarise(
      n_plots = n_distinct(idx),
      mean_ndvi_sd = mean(ndvi_sd, na.rm = TRUE),
      mean_abs_error = mean(abs(ndvi_error), na.rm = TRUE),
      rmse_error = sqrt(mean(ndvi_error^2, na.rm = TRUE)),
      .groups = "drop"
    )
  
  cat("  NDVI Extraction Error by Heterogeneity Class:\n\n")
  print(error_by_class, row.names = FALSE)
  
  # Error per plot
  error_per_plot <- fuzz_df %>%
    group_by(idx, het_class, ndvi_sd, aglb_Mg_per_ha, tmean, ppt, ndvi_true) %>%
    summarise(
      mae = mean(abs(ndvi_error), na.rm = TRUE),
      ndvi_fuzz_sd = sd(ndvi_fuzz, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Correlation between heterogeneity and error
  cor_het_error <- cor(error_per_plot$ndvi_sd, error_per_plot$mae, use = "complete.obs")
  cat("\n  Correlation (NDVI heterogeneity vs MAE): r =", round(cor_het_error, 3), "\n")
  
  # ===========================================================================
  # 5. MODEL PERFORMANCE BY HETEROGENEITY
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 5: Model Performance by Heterogeneity Class\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Train model on true NDVI
  train_data <- error_per_plot %>%
    filter(!is.na(aglb_Mg_per_ha), !is.na(ndvi_true), !is.na(tmean), !is.na(ppt))
  
  X_train <- as.matrix(train_data %>% select(ndvi_true, tmean, ppt))
  y_train <- train_data$aglb_Mg_per_ha
  
  model <- xgb.train(
    data = xgb.DMatrix(X_train, label = y_train),
    nrounds = 50,
    params = list(objective = "reg:squarederror", eta = 0.1, max_depth = 4),
    verbose = 0
  )
  
  # Predict with true and fuzzed NDVI
  train_data$pred_true <- predict(model, xgb.DMatrix(X_train))
  
  # For fuzzed predictions, use mean across replicates
  fuzz_means <- fuzz_df %>%
    group_by(idx) %>%
    summarise(ndvi_fuzz_mean = mean(ndvi_fuzz, na.rm = TRUE), .groups = "drop")
  
  train_data <- train_data %>%
    left_join(fuzz_means, by = "idx")
  
  X_fuzz <- as.matrix(train_data %>% select(ndvi_fuzz_mean, tmean, ppt))
  train_data$pred_fuzz <- predict(model, xgb.DMatrix(X_fuzz))
  
  # Compute errors
  train_data <- train_data %>%
    mutate(
      error_true = aglb_Mg_per_ha - pred_true,
      error_fuzz = aglb_Mg_per_ha - pred_fuzz,
      delta_error = abs(error_fuzz) - abs(error_true)  # Positive = fuzzing hurts
    )
  
  # RMSE by heterogeneity class
  rmse_by_class <- train_data %>%
    group_by(het_class) %>%
    summarise(
      n = n(),
      rmse_true = sqrt(mean(error_true^2, na.rm = TRUE)),
      rmse_fuzz = sqrt(mean(error_fuzz^2, na.rm = TRUE)),
      delta_rmse = rmse_fuzz - rmse_true,
      mean_delta_error = mean(delta_error, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("  Model RMSE by Heterogeneity Class:\n\n")
  print(rmse_by_class, row.names = FALSE)
  
  # ===========================================================================
  # 6. CREATE FIGURES
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 6: Creating Figures\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # --- Figure 1: NDVI Error vs Heterogeneity Scatter ---
  fig1 <- ggplot(error_per_plot, aes(x = ndvi_sd, y = mae)) +
    geom_point(aes(color = het_class), alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    scale_color_viridis_d(option = "D") +
    labs(
      title = "NDVI Extraction Error vs Local Heterogeneity",
      subtitle = paste0("r = ", round(cor_het_error, 3), 
                       " | Higher heterogeneity → larger fuzzing error"),
      x = "Local NDVI SD (1.6km radius)",
      y = "Mean Absolute NDVI Error",
      color = "Heterogeneity\nClass"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(fs::path(output_dir, "fig1_error_vs_heterogeneity.png"), fig1,
         width = 10, height = 7, dpi = 300)
  cat("  ✓ fig1_error_vs_heterogeneity.png\n")
  
  # --- Figure 2: NDVI Error Distributions by Class ---
  fig2 <- ggplot(fuzz_df, aes(x = ndvi_error, fill = het_class)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~het_class, ncol = 1) +
    scale_fill_viridis_d(option = "D") +
    labs(
      title = "NDVI Error Distribution by Heterogeneity Class",
      subtitle = "High heterogeneity shows wider error distribution",
      x = "NDVI Error (fuzzed - true)", y = "Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  ggsave(fs::path(output_dir, "fig2_error_distributions_by_class.png"), fig2,
         width = 10, height = 8, dpi = 300)
  cat("  ✓ fig2_error_distributions_by_class.png\n")
  
  # --- Figure 3: Model ΔRMSE by Heterogeneity Class ---
  fig3 <- ggplot(rmse_by_class, aes(x = het_class, y = delta_rmse, fill = het_class)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = sprintf("+%.1f", delta_rmse)), vjust = -0.5, size = 5) +
    scale_fill_viridis_d(option = "D") +
    labs(
      title = "Model ΔRMSE by Landscape Heterogeneity",
      subtitle = "ΔRMSE = RMSE(fuzzed) - RMSE(true) | Higher = more degradation",
      x = "Heterogeneity Class", y = "ΔRMSE (Mg/ha)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  ggsave(fs::path(output_dir, "fig3_delta_rmse_by_heterogeneity.png"), fig3,
         width = 8, height = 6, dpi = 300)
  cat("  ✓ fig3_delta_rmse_by_heterogeneity.png\n")
  
  # --- Figure 4: Map of Heterogeneity ---
  het_map_sf <- nefin_sample[het_df$idx, ]
  het_map_sf$ndvi_sd <- het_df$ndvi_sd
  het_map_sf$het_class <- het_df$het_class
  
  states <- if (fs::file_exists(states_path)) {
    st_read(states_path, quiet = TRUE)
  } else NULL
  
  fig4 <- ggplot() +
    geom_sf(data = het_map_sf, aes(color = ndvi_sd), size = 2, alpha = 0.7) +
    scale_color_viridis_c(option = "C") +
    labs(
      title = "Spatial Distribution of NDVI Heterogeneity",
      subtitle = "Higher values = more variable landscape in 1.6km radius",
      color = "NDVI SD"
    ) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  if (!is.null(states)) {
    fig4 <- fig4 + geom_sf(data = states, fill = NA, color = "black", linewidth = 0.4)
  }
  
  ggsave(fs::path(output_dir, "fig4_heterogeneity_map.png"), fig4,
         width = 10, height = 8, dpi = 300)
  cat("  ✓ fig4_heterogeneity_map.png\n")
  
  # --- Figure 5: Summary Dashboard ---
  fig5a <- fig1 + theme(legend.position = "none")
  fig5b <- fig3
  
  fig5 <- fig5a + fig5b +
    plot_annotation(
      title = "Fuzzing Impact Increases with Landscape Heterogeneity",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
  
  ggsave(fs::path(output_dir, "fig5_heterogeneity_summary.png"), fig5,
         width = 14, height = 6, dpi = 300)
  cat("  ✓ fig5_heterogeneity_summary.png\n")
  
  # ===========================================================================
  # 7. SAVE RESULTS
  # ===========================================================================
  
  write_csv(het_df, fs::path(output_dir, "plot_heterogeneity.csv"))
  write_csv(error_by_class, fs::path(output_dir, "error_by_heterogeneity_class.csv"))
  write_csv(rmse_by_class, fs::path(output_dir, "rmse_by_heterogeneity_class.csv"))
  
  cat("\n  ✓ Saved CSV results\n")
  
  # ===========================================================================
  # SUMMARY
  # ===========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  ANALYSIS COMPLETE                                                        ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  cat("KEY FINDINGS:\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  cat("  Correlation (heterogeneity vs NDVI error): r =", round(cor_het_error, 3), "\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  cat("  Model ΔRMSE by heterogeneity class:\n")
  for (i in 1:nrow(rmse_by_class)) {
    cat("    ", as.character(rmse_by_class$het_class[i]), ": +", 
        round(rmse_by_class$delta_rmse[i], 1), " Mg/ha\n")
  }
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  
  cat("\nOutput directory:", output_dir, "\n")
  
  invisible(list(
    het_df = het_df,
    error_by_class = error_by_class,
    rmse_by_class = rmse_by_class,
    correlation = cor_het_error
  ))
}

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  analyze_heterogeneity_effect()
}
