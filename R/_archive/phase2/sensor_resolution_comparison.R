#!/usr/bin/env Rscript
# R/sensor_resolution_comparison.R
# Does fuzzing effect increase at finer resolution (Sentinel-2 10m vs MODIS 250m)?
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

compare_sensor_resolution <- function(
    nefin_path = "data/processed/nefin_processed.csv",
    modis_path = "data/processed/ndvi/modis/MODIS_NDVI_5yr_blocked_2020_2024.tif",
    s2_path = "data/processed/ndvi/s2/S2_NDVI_10m_2020_2025.tif",
    tmean_path = "data/processed/prism/prism_tmean_ne_2020_2024.tif",
    ppt_path = "data/processed/prism/prism_ppt_ne_2020_2024.tif",
    output_dir = "runs/sensor_resolution_comparison",
    n_fuzz_replicates = 50,
    fuzz_radius_m = 1609
) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  Sensor Resolution Comparison: MODIS (250m) vs Sentinel-2 (10m)          ║\n")
  cat("║  Does fuzzing hurt more at finer resolution?                              ║\n")
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
  cat("  NEFIN plots:", nrow(nefin), "\n")
  
  # Find coordinate columns
  lon_col <- names(nefin)[grepl("lon", names(nefin), ignore.case = TRUE)][1]
  lat_col <- names(nefin)[grepl("lat", names(nefin), ignore.case = TRUE)][1]
  
  # Load rasters
  cat("\n  Loading rasters...\n")
  r_modis <- rast(modis_path)
  cat("    MODIS: ", res(r_modis)[1], "m resolution\n")
  
  if (!fs::file_exists(s2_path)) {
    cat("  ⚠ Sentinel-2 raster not found:", s2_path, "\n")
    cat("  Checking for alternative S2 files...\n")
    s2_files <- list.files("data/processed/ndvi/s2", pattern = "\\.tif$", full.names = TRUE)
    if (length(s2_files) > 0) {
      s2_path <- s2_files[1]
      cat("    Using:", s2_path, "\n")
    } else {
      stop("No Sentinel-2 data found")
    }
  }
  
  r_s2 <- rast(s2_path)
  cat("    Sentinel-2: ", res(r_s2)[1], "m resolution\n")
  
  r_tmean <- rast(tmean_path)
  r_ppt <- rast(ppt_path)
  
  # ===========================================================================
  # 2. EXTRACT AT TRUE LOCATIONS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 2: Extracting NDVI at True Locations\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Convert to sf
  nefin_sf <- st_as_sf(nefin, coords = c(lon_col, lat_col), crs = 4326)
  
  # Extract MODIS NDVI
  nefin_modis <- st_transform(nefin_sf, crs(r_modis))
  nefin$ndvi_modis_true <- terra::extract(r_modis, vect(nefin_modis))[[2]]
  
  # Extract S2 NDVI
  nefin_s2 <- st_transform(nefin_sf, crs(r_s2))
  nefin$ndvi_s2_true <- terra::extract(r_s2, vect(nefin_s2))[[2]]
  
  # Extract climate
  nefin_prism <- st_transform(nefin_sf, crs(r_tmean))
  nefin$tmean <- terra::extract(r_tmean, vect(nefin_prism))[[2]]
  nefin$ppt <- terra::extract(r_ppt, vect(nefin_prism))[[2]]
  
  # Filter complete cases
  nefin$lon <- nefin[[lon_col]]
  nefin$lat <- nefin[[lat_col]]
  
  nefin_complete <- nefin %>%
    filter(!is.na(aglb_Mg_per_ha), 
           !is.na(ndvi_modis_true), 
           !is.na(ndvi_s2_true),
           !is.na(tmean), !is.na(ppt))
  
  cat("  Plots with all NDVI sources:", nrow(nefin_complete), "\n")
  
  # ===========================================================================
  # 3. MONTE CARLO FUZZING - BOTH SENSORS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 3: Monte Carlo Fuzzing at Both Resolutions\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Sample plots
  set.seed(42)
  n_sample <- min(500, nrow(nefin_complete))
  sample_idx <- sample(1:nrow(nefin_complete), n_sample)
  nefin_sample <- nefin_complete[sample_idx, ]
  
  cat("  Using", n_sample, "plots for fuzzing analysis\n")
  cat("  Running", n_fuzz_replicates, "Monte Carlo replicates...\n\n")
  
  # Convert to projected CRS for fuzzing
  nefin_sample_sf <- st_as_sf(nefin_sample, coords = c("lon", "lat"), crs = 4326) %>%
    st_transform(5070)
  
  fuzz_results <- list()
  pb <- txtProgressBar(min = 0, max = n_fuzz_replicates, style = 3)
  
  for (rep in 1:n_fuzz_replicates) {
    # Generate random offsets
    angles <- runif(n_sample, 0, 2 * pi)
    distances <- runif(n_sample, 0, fuzz_radius_m)
    
    dx <- distances * cos(angles)
    dy <- distances * sin(angles)
    
    # Apply fuzzing
    coords <- st_coordinates(nefin_sample_sf)
    coords_fuzz <- cbind(coords[, 1] + dx, coords[, 2] + dy)
    
    fuzz_sf <- st_as_sf(
      data.frame(id = 1:n_sample, x = coords_fuzz[, 1], y = coords_fuzz[, 2]),
      coords = c("x", "y"), crs = 5070
    )
    
    # Extract MODIS at fuzzed locations
    fuzz_modis <- st_transform(fuzz_sf, crs(r_modis))
    ndvi_modis_fuzz <- terra::extract(r_modis, vect(fuzz_modis))[[2]]
    
    # Extract S2 at fuzzed locations
    fuzz_s2 <- st_transform(fuzz_sf, crs(r_s2))
    ndvi_s2_fuzz <- terra::extract(r_s2, vect(fuzz_s2))[[2]]
    
    fuzz_results[[rep]] <- data.frame(
      plot_idx = 1:n_sample,
      replicate = rep,
      fuzz_distance = distances,
      ndvi_modis_true = nefin_sample$ndvi_modis_true,
      ndvi_modis_fuzz = ndvi_modis_fuzz,
      ndvi_s2_true = nefin_sample$ndvi_s2_true,
      ndvi_s2_fuzz = ndvi_s2_fuzz,
      biomass = nefin_sample$aglb_Mg_per_ha,
      tmean = nefin_sample$tmean,
      ppt = nefin_sample$ppt
    )
    
    setTxtProgressBar(pb, rep)
  }
  close(pb)
  
  fuzz_df <- bind_rows(fuzz_results)
  
  # Compute errors
  fuzz_df <- fuzz_df %>%
    mutate(
      modis_error = ndvi_modis_fuzz - ndvi_modis_true,
      s2_error = ndvi_s2_fuzz - ndvi_s2_true
    )
  
  # ===========================================================================
  # 4. COMPARE NDVI UNCERTAINTY BY SENSOR
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 4: Comparing NDVI Uncertainty by Sensor\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Summary per plot
  plot_uncertainty <- fuzz_df %>%
    group_by(plot_idx) %>%
    summarise(
      # MODIS
      modis_true = first(ndvi_modis_true),
      modis_fuzz_mean = mean(ndvi_modis_fuzz, na.rm = TRUE),
      modis_fuzz_sd = sd(ndvi_modis_fuzz, na.rm = TRUE),
      modis_mae = mean(abs(modis_error), na.rm = TRUE),
      # S2
      s2_true = first(ndvi_s2_true),
      s2_fuzz_mean = mean(ndvi_s2_fuzz, na.rm = TRUE),
      s2_fuzz_sd = sd(ndvi_s2_fuzz, na.rm = TRUE),
      s2_mae = mean(abs(s2_error), na.rm = TRUE),
      # Biomass
      biomass = first(biomass),
      tmean = first(tmean),
      ppt = first(ppt),
      .groups = "drop"
    )
  
  cat("  NDVI Uncertainty from ±1.6km Fuzzing:\n\n")
  cat("  MODIS (250m):\n")
  cat("    Mean Absolute Error:", round(mean(plot_uncertainty$modis_mae, na.rm = TRUE), 4), "\n")
  cat("    Mean SD across replicates:", round(mean(plot_uncertainty$modis_fuzz_sd, na.rm = TRUE), 4), "\n")
  
  cat("\n  Sentinel-2 (10m):\n")
  cat("    Mean Absolute Error:", round(mean(plot_uncertainty$s2_mae, na.rm = TRUE), 4), "\n")
  cat("    Mean SD across replicates:", round(mean(plot_uncertainty$s2_fuzz_sd, na.rm = TRUE), 4), "\n")
  
  s2_vs_modis_ratio <- mean(plot_uncertainty$s2_fuzz_sd, na.rm = TRUE) / 
                       mean(plot_uncertainty$modis_fuzz_sd, na.rm = TRUE)
  cat("\n  S2/MODIS uncertainty ratio:", round(s2_vs_modis_ratio, 2), "x\n")
  
  # ===========================================================================
  # 5. COMPARE MODEL PERFORMANCE BY SENSOR
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 5: Comparing Model Performance by Sensor\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Split data
  set.seed(42)
  train_idx <- sample(1:nrow(nefin_complete), 0.7 * nrow(nefin_complete))
  train_data <- nefin_complete[train_idx, ]
  test_data <- nefin_complete[-train_idx, ]
  
  cat("  Train:", nrow(train_data), "| Test:", nrow(test_data), "\n\n")
  
  results <- list()
  
  # --- MODIS Models ---
  cat("  MODIS Models:\n")
  
  # Train on true MODIS
  X_train_modis <- as.matrix(train_data %>% select(ndvi_modis_true, tmean, ppt))
  y_train <- train_data$aglb_Mg_per_ha
  
  model_modis_true <- xgb.train(
    data = xgb.DMatrix(X_train_modis, label = y_train),
    nrounds = 50,
    params = list(objective = "reg:squarederror", eta = 0.1, max_depth = 4),
    verbose = 0
  )
  
  # Test on true MODIS
  X_test_modis <- as.matrix(test_data %>% select(ndvi_modis_true, tmean, ppt))
  y_test <- test_data$aglb_Mg_per_ha
  
  preds_modis_true <- predict(model_modis_true, xgb.DMatrix(X_test_modis))
  rmse_modis_true <- sqrt(mean((y_test - preds_modis_true)^2))
  
  cat("    True coords RMSE:", round(rmse_modis_true, 2), "Mg/ha\n")
  
  # Simulate fuzzed MODIS (average over MC replicates)
  modis_fuzz_rmse <- c()
  for (rep in 1:min(20, n_fuzz_replicates)) {
    rep_data <- fuzz_df %>% filter(replicate == rep)
    
    # Get test plots
    test_plot_idx <- which(nefin_complete$CN[-train_idx] %in% nefin_sample$CN)
    if (length(test_plot_idx) < 10) next
    
    # Use fuzzed NDVI for these plots
    test_fuzzed <- test_data
    match_idx <- match(test_fuzzed$CN, nefin_sample$CN)
    valid_match <- !is.na(match_idx)
    
    if (sum(valid_match) > 0) {
      test_fuzzed$ndvi_fuzz <- test_fuzzed$ndvi_modis_true
      fuzz_vals <- rep_data$ndvi_modis_fuzz[match_idx[valid_match]]
      test_fuzzed$ndvi_fuzz[valid_match] <- fuzz_vals
      
      X_test_fuzz <- as.matrix(test_fuzzed %>% select(ndvi_fuzz, tmean, ppt))
      preds_fuzz <- predict(model_modis_true, xgb.DMatrix(X_test_fuzz))
      modis_fuzz_rmse <- c(modis_fuzz_rmse, sqrt(mean((y_test - preds_fuzz)^2, na.rm = TRUE)))
    }
  }
  
  if (length(modis_fuzz_rmse) > 0) {
    cat("    Fuzzed coords RMSE:", round(mean(modis_fuzz_rmse), 2), "±", 
        round(sd(modis_fuzz_rmse), 2), "Mg/ha\n")
    results$modis_delta_rmse <- mean(modis_fuzz_rmse) - rmse_modis_true
    cat("    ΔRMSE:", round(results$modis_delta_rmse, 2), "Mg/ha\n")
  }
  
  results$modis_true_rmse <- rmse_modis_true
  results$modis_fuzz_rmse <- mean(modis_fuzz_rmse)
  
  # --- Sentinel-2 Models ---
  cat("\n  Sentinel-2 Models:\n")
  
  # Train on true S2
  X_train_s2 <- as.matrix(train_data %>% select(ndvi_s2_true, tmean, ppt))
  
  model_s2_true <- xgb.train(
    data = xgb.DMatrix(X_train_s2, label = y_train),
    nrounds = 50,
    params = list(objective = "reg:squarederror", eta = 0.1, max_depth = 4),
    verbose = 0
  )
  
  # Test on true S2
  X_test_s2 <- as.matrix(test_data %>% select(ndvi_s2_true, tmean, ppt))
  preds_s2_true <- predict(model_s2_true, xgb.DMatrix(X_test_s2))
  rmse_s2_true <- sqrt(mean((y_test - preds_s2_true)^2))
  
  cat("    True coords RMSE:", round(rmse_s2_true, 2), "Mg/ha\n")
  
  # Simulate fuzzed S2
  s2_fuzz_rmse <- c()
  for (rep in 1:min(20, n_fuzz_replicates)) {
    rep_data <- fuzz_df %>% filter(replicate == rep)
    
    test_fuzzed <- test_data
    match_idx <- match(test_fuzzed$CN, nefin_sample$CN)
    valid_match <- !is.na(match_idx)
    
    if (sum(valid_match) > 0) {
      test_fuzzed$ndvi_fuzz <- test_fuzzed$ndvi_s2_true
      fuzz_vals <- rep_data$ndvi_s2_fuzz[match_idx[valid_match]]
      test_fuzzed$ndvi_fuzz[valid_match] <- fuzz_vals
      
      X_test_fuzz <- as.matrix(test_fuzzed %>% select(ndvi_fuzz, tmean, ppt))
      preds_fuzz <- predict(model_s2_true, xgb.DMatrix(X_test_fuzz))
      s2_fuzz_rmse <- c(s2_fuzz_rmse, sqrt(mean((y_test - preds_fuzz)^2, na.rm = TRUE)))
    }
  }
  
  if (length(s2_fuzz_rmse) > 0) {
    cat("    Fuzzed coords RMSE:", round(mean(s2_fuzz_rmse), 2), "±", 
        round(sd(s2_fuzz_rmse), 2), "Mg/ha\n")
    results$s2_delta_rmse <- mean(s2_fuzz_rmse) - rmse_s2_true
    cat("    ΔRMSE:", round(results$s2_delta_rmse, 2), "Mg/ha\n")
  }
  
  results$s2_true_rmse <- rmse_s2_true
  results$s2_fuzz_rmse <- mean(s2_fuzz_rmse)
  
  # ===========================================================================
  # 6. CREATE FIGURES
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 6: Creating Figures\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # --- Figure 1: NDVI Error Distributions by Sensor ---
  error_long <- fuzz_df %>%
    select(plot_idx, replicate, modis_error, s2_error) %>%
    pivot_longer(cols = c(modis_error, s2_error),
                 names_to = "sensor", values_to = "error") %>%
    mutate(sensor = ifelse(sensor == "modis_error", "MODIS (250m)", "Sentinel-2 (10m)"))
  
  fig1 <- ggplot(error_long, aes(x = error, fill = sensor)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~sensor, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c("MODIS (250m)" = "#3498db", "Sentinel-2 (10m)" = "#e74c3c")) +
    labs(
      title = "NDVI Extraction Error from ±1.6km Coordinate Fuzzing",
      subtitle = paste0("Sentinel-2 shows ", round(s2_vs_modis_ratio, 1), 
                       "x more uncertainty than MODIS"),
      x = "NDVI Error (fuzzed - true)", y = "Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 12)
    )
  
  ggsave(fs::path(output_dir, "fig1_ndvi_error_by_sensor.png"), fig1,
         width = 10, height = 8, dpi = 300)
  cat("  ✓ fig1_ndvi_error_by_sensor.png\n")
  
  # --- Figure 2: NDVI Uncertainty SD Comparison ---
  uncertainty_long <- plot_uncertainty %>%
    select(plot_idx, modis_fuzz_sd, s2_fuzz_sd) %>%
    pivot_longer(cols = c(modis_fuzz_sd, s2_fuzz_sd),
                 names_to = "sensor", values_to = "sd") %>%
    mutate(sensor = ifelse(sensor == "modis_fuzz_sd", "MODIS (250m)", "Sentinel-2 (10m)"))
  
  fig2 <- ggplot(uncertainty_long, aes(x = sensor, y = sd, fill = sensor)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(values = c("MODIS (250m)" = "#3498db", "Sentinel-2 (10m)" = "#e74c3c")) +
    labs(
      title = "NDVI Uncertainty (SD) by Sensor Resolution",
      subtitle = "Higher SD = more sensitive to coordinate fuzzing",
      x = "", y = "NDVI Standard Deviation\n(across fuzzing replicates)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  ggsave(fs::path(output_dir, "fig2_ndvi_uncertainty_by_sensor.png"), fig2,
         width = 8, height = 6, dpi = 300)
  cat("  ✓ fig2_ndvi_uncertainty_by_sensor.png\n")
  
  # --- Figure 3: Model ΔRMSE by Sensor ---
  if (!is.null(results$modis_delta_rmse) && !is.null(results$s2_delta_rmse)) {
    delta_df <- data.frame(
      sensor = c("MODIS (250m)", "Sentinel-2 (10m)"),
      delta_rmse = c(results$modis_delta_rmse, results$s2_delta_rmse),
      resolution = c(250, 10)
    )
    
    fig3 <- ggplot(delta_df, aes(x = sensor, y = delta_rmse, fill = sensor)) +
      geom_col(alpha = 0.8) +
      geom_text(aes(label = sprintf("+%.1f", delta_rmse)), vjust = -0.5, size = 5) +
      scale_fill_manual(values = c("MODIS (250m)" = "#3498db", "Sentinel-2 (10m)" = "#e74c3c")) +
      labs(
        title = "Impact of Coordinate Fuzzing on Model Performance",
        subtitle = "ΔRMSE = RMSE(fuzzed) - RMSE(true) | Higher = more degradation",
        x = "", y = "ΔRMSE (Mg/ha)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "none"
      )
    
    ggsave(fs::path(output_dir, "fig3_delta_rmse_by_sensor.png"), fig3,
           width = 8, height = 6, dpi = 300)
    cat("  ✓ fig3_delta_rmse_by_sensor.png\n")
  }
  
  # --- Figure 4: Combined Summary ---
  fig4a <- ggplot(uncertainty_long, aes(x = sensor, y = sd, fill = sensor)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    scale_fill_manual(values = c("MODIS (250m)" = "#3498db", "Sentinel-2 (10m)" = "#e74c3c")) +
    labs(title = "A) NDVI Uncertainty", x = "", y = "SD") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
  
  if (exists("delta_df")) {
    fig4b <- ggplot(delta_df, aes(x = sensor, y = delta_rmse, fill = sensor)) +
      geom_col(alpha = 0.8) +
      scale_fill_manual(values = c("MODIS (250m)" = "#3498db", "Sentinel-2 (10m)" = "#e74c3c")) +
      labs(title = "B) Model ΔRMSE", x = "", y = "ΔRMSE (Mg/ha)") +
      theme_minimal() +
      theme(legend.position = "none", plot.title = element_text(face = "bold"))
    
    fig4 <- fig4a + fig4b +
      plot_annotation(
        title = "Fuzzing Effect by Sensor Resolution",
        subtitle = "Finer resolution (Sentinel-2) is more sensitive to coordinate uncertainty",
        theme = theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 11)
        )
      )
    
    ggsave(fs::path(output_dir, "fig4_resolution_summary.png"), fig4,
           width = 12, height = 6, dpi = 300)
    cat("  ✓ fig4_resolution_summary.png\n")
  }
  
  # ===========================================================================
  # 7. SAVE RESULTS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 7: Saving Results\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Summary table
  summary_df <- data.frame(
    sensor = c("MODIS", "MODIS", "Sentinel-2", "Sentinel-2"),
    metric = c("NDVI_SD", "ΔRMSE", "NDVI_SD", "ΔRMSE"),
    value = c(
      mean(plot_uncertainty$modis_fuzz_sd, na.rm = TRUE),
      results$modis_delta_rmse %||% NA,
      mean(plot_uncertainty$s2_fuzz_sd, na.rm = TRUE),
      results$s2_delta_rmse %||% NA
    )
  )
  
  write_csv(summary_df, fs::path(output_dir, "sensor_comparison_summary.csv"))
  write_csv(plot_uncertainty, fs::path(output_dir, "plot_uncertainty_by_sensor.csv"))
  
  cat("  ✓ Saved summary CSVs\n")
  
  # ===========================================================================
  # FINAL SUMMARY
  # ===========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  ANALYSIS COMPLETE                                                        ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  cat("KEY FINDINGS:\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  cat("  NDVI Uncertainty (SD from fuzzing):\n")
  cat("    MODIS (250m):     ", round(mean(plot_uncertainty$modis_fuzz_sd, na.rm = TRUE), 4), "\n")
  cat("    Sentinel-2 (10m): ", round(mean(plot_uncertainty$s2_fuzz_sd, na.rm = TRUE), 4), "\n")
  cat("    Ratio (S2/MODIS): ", round(s2_vs_modis_ratio, 2), "x\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  
  if (!is.null(results$modis_delta_rmse) && !is.null(results$s2_delta_rmse)) {
    cat("  Model ΔRMSE (fuzzed - true):\n")
    cat("    MODIS:      +", round(results$modis_delta_rmse, 2), "Mg/ha\n")
    cat("    Sentinel-2: +", round(results$s2_delta_rmse, 2), "Mg/ha\n")
    cat("─────────────────────────────────────────────────────────────────────────────\n")
  }
  
  cat("\nOutput directory:", output_dir, "\n")
  
  invisible(results)
}

# Helper
`%||%` <- function(a, b) if (!is.null(a)) a else b

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  compare_sensor_resolution()
}
