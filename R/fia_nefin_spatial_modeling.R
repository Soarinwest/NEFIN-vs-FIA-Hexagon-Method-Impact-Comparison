#!/usr/bin/env Rscript
# R/fia_nefin_spatial_modeling.R
# Comprehensive FIA vs NEFIN spatial modeling comparison
# Does using true coordinates improve biomass prediction accuracy?
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
  library(scales)
  library(fs)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

config <- list(
  # Data paths
  nefin_path = "data/processed/nefin_processed.csv",
  fia_covariates_path = "data/processed/fia_covariates.csv",
  fia_plots_path = "data/processed/plot_hex_assignments.csv",
  jitter_dir = "data/processed/mc_jitter_library/replicates",
  
 # Raster paths (by window)
  rasters = list(
    "2015_2019" = list(
      ndvi_modis = "data/processed/ndvi/modis/MODIS_NDVI_5yr_blocked_2015_2019.tif",
      ndvi_s2 = "data/processed/ndvi/s2/S2_NDVI_10m_2016_2019.tif",
      tmean = "data/processed/prism/prism_tmean_ne_2015_2019.tif",
      ppt = "data/processed/prism/prism_ppt_ne_2015_2019.tif"
    ),
    "2020_2024" = list(
      ndvi_modis = "data/processed/ndvi/modis/MODIS_NDVI_5yr_blocked_2020_2024.tif",
      ndvi_s2 = "data/processed/ndvi/s2/S2_NDVI_10m_2020_2025.tif",
      tmean = "data/processed/prism/prism_tmean_ne_2020_2024.tif",
      ppt = "data/processed/prism/prism_ppt_ne_2020_2024.tif"
    )
  ),
  
  # Hex grids
  hex_path = "data/hex/hex_grid_50kha.geojson",
  states_path = "data/boundaries/states_5070.geojson",
  
  # Analysis settings
  n_mc_replicates = 20,  # Use subset of MC replicates for speed
  n_cv_folds = 5,
  holdout_fraction = 0.3,
  time_window = "2020_2024",
  
  # Output
  output_dir = "runs/fia_nefin_spatial_modeling"
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

extract_covariates_at_points <- function(pts_sf, rasters, window) {
  # Extract NDVI, tmean, ppt at point locations
  
  r_list <- rasters[[window]]
  
  result <- pts_sf %>% st_drop_geometry()
  
  # MODIS NDVI
  if (!is.null(r_list$ndvi_modis) && fs::file_exists(r_list$ndvi_modis)) {
    r <- terra::rast(r_list$ndvi_modis)
    pts_reproj <- st_transform(pts_sf, crs(r))
    vals <- terra::extract(r, terra::vect(pts_reproj))
    result$ndvi_modis <- vals[[2]]
  }
  
  # Sentinel-2 NDVI
  if (!is.null(r_list$ndvi_s2) && fs::file_exists(r_list$ndvi_s2)) {
    r <- terra::rast(r_list$ndvi_s2)
    pts_reproj <- st_transform(pts_sf, crs(r))
    vals <- terra::extract(r, terra::vect(pts_reproj))
    result$ndvi_s2 <- vals[[2]]
  }
  
  # Temperature
  if (!is.null(r_list$tmean) && fs::file_exists(r_list$tmean)) {
    r <- terra::rast(r_list$tmean)
    pts_reproj <- st_transform(pts_sf, crs(r))
    vals <- terra::extract(r, terra::vect(pts_reproj))
    result$tmean <- vals[[2]]
  }
  
  # Precipitation
  if (!is.null(r_list$ppt) && fs::file_exists(r_list$ppt)) {
    r <- terra::rast(r_list$ppt)
    pts_reproj <- st_transform(pts_sf, crs(r))
    vals <- terra::extract(r, terra::vect(pts_reproj))
    result$ppt <- vals[[2]]
  }
  
  return(result)
}

create_spatial_folds <- function(pts_sf, n_folds = 5, seed = 42) {
  # Create spatially-blocked CV folds using k-means on coordinates
  set.seed(seed)
  
  coords <- st_coordinates(pts_sf)
  km <- kmeans(coords, centers = n_folds, nstart = 10)
  
  return(km$cluster)
}

fit_xgboost_cv <- function(X, y, folds, params = NULL) {
  # Fit XGBoost with spatial CV, return fold-wise metrics
  
  if (is.null(params)) {
    params <- list(
      objective = "reg:squarederror",
      eta = 0.1,
      max_depth = 4,
      subsample = 0.8,
      colsample_bytree = 0.8,
      min_child_weight = 5
    )
  }
  
  unique_folds <- sort(unique(folds))
  fold_results <- list()
  all_preds <- rep(NA, length(y))
  
  for (fold in unique_folds) {
    train_idx <- which(folds != fold)
    test_idx <- which(folds == fold)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_test <- X[test_idx, , drop = FALSE]
    y_test <- y[test_idx]
    
    dtrain <- xgb.DMatrix(data = X_train, label = y_train)
    dtest <- xgb.DMatrix(data = X_test, label = y_test)
    
    model <- xgb.train(
      data = dtrain,
      nrounds = 50,
      params = params,
      verbose = 0
    )
    
    preds <- predict(model, dtest)
    all_preds[test_idx] <- preds
    
    fold_results[[fold]] <- list(
      fold = fold,
      n_test = length(y_test),
      rmse = sqrt(mean((y_test - preds)^2)),
      mae = mean(abs(y_test - preds)),
      r2 = 1 - sum((y_test - preds)^2) / sum((y_test - mean(y_test))^2),
      bias = mean(preds - y_test)
    )
  }
  
  fold_df <- bind_rows(fold_results)
  
  return(list(
    fold_metrics = fold_df,
    predictions = all_preds,
    overall_rmse = sqrt(mean((y - all_preds)^2, na.rm = TRUE)),
    overall_mae = mean(abs(y - all_preds), na.rm = TRUE),
    overall_r2 = 1 - sum((y - all_preds)^2, na.rm = TRUE) / sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  ))
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

run_fia_nefin_modeling <- function(cfg = config) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  FIA vs NEFIN Spatial Modeling Comparison                                 ║\n")
  cat("║  Does using true coordinates improve biomass prediction accuracy?         ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  fs::dir_create(cfg$output_dir, recurse = TRUE)
  fs::dir_create(fs::path(cfg$output_dir, "figures"), recurse = TRUE)
  fs::dir_create(fs::path(cfg$output_dir, "maps"), recurse = TRUE)
  
  window <- cfg$time_window
  cat("Time window:", window, "\n\n")
  
  # ===========================================================================
  # 1. LOAD AND PREPARE DATA
  # ===========================================================================
  
  cat("═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 1: Loading and Preparing Data\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # --- Load NEFIN (true coordinates) ---
  nefin <- read_csv(cfg$nefin_path, show_col_types = FALSE)
  cat("  NEFIN plots:", nrow(nefin), "\n")
  
  # Find coordinate columns
  lon_col <- names(nefin)[grepl("lon", names(nefin), ignore.case = TRUE)][1]
  lat_col <- names(nefin)[grepl("lat", names(nefin), ignore.case = TRUE)][1]
  
  # Convert to sf
  nefin_sf <- st_as_sf(nefin, coords = c(lon_col, lat_col), crs = 4326)
  
  # Extract covariates at true locations
  cat("  Extracting covariates at NEFIN true locations...\n")
  nefin_data <- extract_covariates_at_points(nefin_sf, cfg$rasters, window)
  nefin_data$lon <- nefin[[lon_col]]
  nefin_data$lat <- nefin[[lat_col]]
  
  # Filter complete cases
  nefin_complete <- nefin_data %>%
    filter(!is.na(aglb_Mg_per_ha), !is.na(ndvi_modis), !is.na(tmean), !is.na(ppt))
  
  cat("  NEFIN complete cases:", nrow(nefin_complete), "\n")
  
  # --- Load FIA plots and compute biomass ---
  cat("\n  Loading FIA data...\n")
  
  fia_plots <- read_csv(cfg$fia_plots_path, show_col_types = FALSE)
  
  # Load tree data to compute biomass
  tree_files <- list.files("data/interim/fia/states", pattern = "tree.csv",
                           recursive = TRUE, full.names = TRUE)
  
  if (length(tree_files) > 0) {
    trees <- lapply(tree_files, function(f) {
      read_csv(f, show_col_types = FALSE) %>%
        select(any_of(c("PLT_CN", "DRYBIO_AG", "TPA_UNADJ")))
    }) %>% bind_rows()
    
    plot_biomass <- trees %>%
      filter(!is.na(DRYBIO_AG), !is.na(TPA_UNADJ)) %>%
      group_by(PLT_CN) %>%
      summarise(
        aglb_Mg_per_ha = sum(DRYBIO_AG * TPA_UNADJ * 0.0004536 * 2.471, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      rename(CN = PLT_CN)
    
    fia_plots <- fia_plots %>%
      left_join(plot_biomass, by = "CN")
  }
  
  # Filter to time window
  year_start <- as.integer(substr(window, 1, 4))
  year_end <- as.integer(substr(window, 6, 9))
  
  fia_plots <- fia_plots %>%
    filter(MEASYEAR >= year_start, MEASYEAR <= year_end,
           !is.na(aglb_Mg_per_ha))
  
  cat("  FIA plots in window:", nrow(fia_plots), "\n")
  
  # --- Load MC jitter replicates ---
  jitter_files <- list.files(cfg$jitter_dir, pattern = "rep_.*\\.csv$", full.names = TRUE)
  n_reps <- min(cfg$n_mc_replicates, length(jitter_files))
  jitter_files <- jitter_files[1:n_reps]
  
  cat("  Using", n_reps, "MC jitter replicates\n")
  
  # ===========================================================================
  # 2. TRAIN MODELS: NEFIN (TRUE COORDS)
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 2: Training Models on NEFIN (True Coordinates)\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Split NEFIN into train/holdout
  set.seed(42)
  n_nefin <- nrow(nefin_complete)
  holdout_idx <- sample(1:n_nefin, size = round(cfg$holdout_fraction * n_nefin))
  
  nefin_train <- nefin_complete[-holdout_idx, ]
  nefin_holdout <- nefin_complete[holdout_idx, ]
  
  cat("  NEFIN train:", nrow(nefin_train), "\n")
  cat("  NEFIN holdout:", nrow(nefin_holdout), "\n")
  
  # Convert to sf for spatial folds
  nefin_train_sf <- st_as_sf(nefin_train, coords = c("lon", "lat"), crs = 4326)
  
  # Create spatial folds
  folds_nefin <- create_spatial_folds(nefin_train_sf, n_folds = cfg$n_cv_folds)
  
  # Prepare matrices
  X_nefin <- as.matrix(nefin_train %>% select(ndvi_modis, tmean, ppt))
  y_nefin <- nefin_train$aglb_Mg_per_ha
  
  # Fit with spatial CV
  cat("\n  Fitting XGBoost with spatial CV...\n")
  nefin_results <- fit_xgboost_cv(X_nefin, y_nefin, folds_nefin)
  
  cat("  NEFIN XGBoost Results:\n")
  cat("    Overall RMSE:", round(nefin_results$overall_rmse, 2), "Mg/ha\n")
  cat("    Overall R²:", round(nefin_results$overall_r2, 4), "\n")
  cat("    Fold RMSE range:", round(min(nefin_results$fold_metrics$rmse), 2), "-", 
      round(max(nefin_results$fold_metrics$rmse), 2), "\n")
  
  # Train final model on all NEFIN train data
  dtrain_nefin <- xgb.DMatrix(data = X_nefin, label = y_nefin)
  model_nefin <- xgb.train(
    data = dtrain_nefin,
    nrounds = 50,
    params = list(objective = "reg:squarederror", eta = 0.1, max_depth = 4),
    verbose = 0
  )
  
  # ===========================================================================
  # 3. TRAIN MODELS: FIA (FUZZED COORDS) - MONTE CARLO
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 3: Training Models on FIA (Fuzzed Coordinates) - Monte Carlo\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  mc_results <- list()
  
  pb <- txtProgressBar(min = 0, max = n_reps, style = 3)
  
  for (rep in 1:n_reps) {
    # Load jittered coordinates
    jitter_data <- read_csv(jitter_files[rep], show_col_types = FALSE)
    
    # Filter to window and plots with biomass
    jitter_data <- jitter_data %>%
      filter(MEASYEAR >= year_start, MEASYEAR <= year_end) %>%
      inner_join(plot_biomass, by = "CN")
    
    if (nrow(jitter_data) < 100) {
      next
    }
    
    # Convert to sf and extract covariates
    jitter_sf <- st_as_sf(jitter_data, 
                          coords = c("lon_jittered", "lat_jittered"), 
                          crs = 4326)
    
    jitter_extracted <- extract_covariates_at_points(jitter_sf, cfg$rasters, window)
    jitter_extracted$aglb_Mg_per_ha <- jitter_data$aglb_Mg_per_ha
    jitter_extracted$lon <- jitter_data$lon_jittered
    jitter_extracted$lat <- jitter_data$lat_jittered
    
    # Filter complete cases
    jitter_complete <- jitter_extracted %>%
      filter(!is.na(aglb_Mg_per_ha), !is.na(ndvi_modis), !is.na(tmean), !is.na(ppt))
    
    if (nrow(jitter_complete) < 100) {
      next
    }
    
    # Create spatial folds
    jitter_sf_complete <- st_as_sf(jitter_complete, coords = c("lon", "lat"), crs = 4326)
    folds_fia <- create_spatial_folds(jitter_sf_complete, n_folds = cfg$n_cv_folds, seed = 42 + rep)
    
    # Fit model
    X_fia <- as.matrix(jitter_complete %>% select(ndvi_modis, tmean, ppt))
    y_fia <- jitter_complete$aglb_Mg_per_ha
    
    fia_cv <- fit_xgboost_cv(X_fia, y_fia, folds_fia)
    
    mc_results[[rep]] <- list(
      replicate = rep,
      n_plots = nrow(jitter_complete),
      overall_rmse = fia_cv$overall_rmse,
      overall_mae = fia_cv$overall_mae,
      overall_r2 = fia_cv$overall_r2,
      fold_metrics = fia_cv$fold_metrics
    )
    
    setTxtProgressBar(pb, rep)
  }
  close(pb)
  
  # Summarize MC results
  mc_summary <- bind_rows(lapply(mc_results, function(x) {
    data.frame(replicate = x$replicate, rmse = x$overall_rmse, 
               mae = x$overall_mae, r2 = x$overall_r2, n = x$n_plots)
  }))
  
  cat("\n  FIA (Fuzzed) MC Results:\n")
  cat("    Mean RMSE:", round(mean(mc_summary$rmse), 2), "±", round(sd(mc_summary$rmse), 2), "Mg/ha\n")
  cat("    Mean R²:", round(mean(mc_summary$r2), 4), "\n")
  cat("    RMSE range:", round(min(mc_summary$rmse), 2), "-", round(max(mc_summary$rmse), 2), "\n")
  
  # ===========================================================================
  # 4. EVALUATE ON NEFIN HOLDOUT
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 4: Evaluating on NEFIN Holdout (True Biomass)\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  X_holdout <- as.matrix(nefin_holdout %>% select(ndvi_modis, tmean, ppt))
  y_holdout <- nefin_holdout$aglb_Mg_per_ha
  
  # NEFIN-trained model predictions
  preds_nefin <- predict(model_nefin, xgb.DMatrix(X_holdout))
  rmse_nefin_holdout <- sqrt(mean((y_holdout - preds_nefin)^2))
  r2_nefin_holdout <- 1 - sum((y_holdout - preds_nefin)^2) / sum((y_holdout - mean(y_holdout))^2)
  
  cat("  NEFIN-trained model on holdout:\n")
  cat("    RMSE:", round(rmse_nefin_holdout, 2), "Mg/ha\n")
  cat("    R²:", round(r2_nefin_holdout, 4), "\n")
  
  # FIA-trained models on holdout (MC)
  fia_holdout_rmse <- c()
  
  for (rep in 1:min(5, n_reps)) {  # Use first 5 reps for speed
    jitter_data <- read_csv(jitter_files[rep], show_col_types = FALSE) %>%
      filter(MEASYEAR >= year_start, MEASYEAR <= year_end) %>%
      inner_join(plot_biomass, by = "CN")
    
    jitter_sf <- st_as_sf(jitter_data, coords = c("lon_jittered", "lat_jittered"), crs = 4326)
    jitter_extracted <- extract_covariates_at_points(jitter_sf, cfg$rasters, window)
    jitter_extracted$aglb_Mg_per_ha <- jitter_data$aglb_Mg_per_ha
    
    jitter_complete <- jitter_extracted %>%
      filter(!is.na(aglb_Mg_per_ha), !is.na(ndvi_modis), !is.na(tmean), !is.na(ppt))
    
    X_fia <- as.matrix(jitter_complete %>% select(ndvi_modis, tmean, ppt))
    y_fia <- jitter_complete$aglb_Mg_per_ha
    
    model_fia <- xgb.train(
      data = xgb.DMatrix(X_fia, label = y_fia),
      nrounds = 50,
      params = list(objective = "reg:squarederror", eta = 0.1, max_depth = 4),
      verbose = 0
    )
    
    preds_fia <- predict(model_fia, xgb.DMatrix(X_holdout))
    fia_holdout_rmse <- c(fia_holdout_rmse, sqrt(mean((y_holdout - preds_fia)^2)))
  }
  
  cat("\n  FIA-trained models on holdout:\n")
  cat("    Mean RMSE:", round(mean(fia_holdout_rmse), 2), "±", round(sd(fia_holdout_rmse), 2), "Mg/ha\n")
  
  delta_rmse <- mean(fia_holdout_rmse) - rmse_nefin_holdout
  cat("\n  ΔRMSE (FIA - NEFIN):", round(delta_rmse, 2), "Mg/ha\n")
  cat("  → True coordinates improve RMSE by", round(100 * delta_rmse / mean(fia_holdout_rmse), 1), "%\n")
  
  # ===========================================================================
  # 5. CREATE FIGURES
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 5: Creating Figures\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # --- Figure 1: CV RMSE Distribution ---
  cv_comparison <- bind_rows(
    nefin_results$fold_metrics %>% mutate(dataset = "NEFIN (true)"),
    bind_rows(lapply(mc_results, function(x) x$fold_metrics)) %>% mutate(dataset = "FIA (fuzzed)")
  )
  
  fig1 <- ggplot(cv_comparison, aes(x = dataset, y = rmse, fill = dataset)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_fill_manual(values = c("FIA (fuzzed)" = "#e74c3c", "NEFIN (true)" = "#27ae60")) +
    labs(
      title = "Cross-Validation RMSE: FIA (Fuzzed) vs NEFIN (True)",
      subtitle = paste0("Spatial CV with ", cfg$n_cv_folds, " folds | ", n_reps, " MC replicates for FIA"),
      x = "", y = "RMSE (Mg/ha)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")
  
  ggsave(fs::path(cfg$output_dir, "figures", "fig1_cv_rmse_comparison.png"), fig1,
         width = 8, height = 6, dpi = 300)
  cat("  ✓ fig1_cv_rmse_comparison.png\n")
  
  # --- Figure 2: MC Variability ---
  fig2 <- ggplot(mc_summary, aes(x = rmse)) +
    geom_histogram(bins = 15, fill = "#e74c3c", alpha = 0.7) +
    geom_vline(xintercept = nefin_results$overall_rmse, 
               linetype = "dashed", color = "#27ae60", linewidth = 1.2) +
    annotate("text", x = nefin_results$overall_rmse, y = Inf, 
             label = paste0("NEFIN: ", round(nefin_results$overall_rmse, 1)),
             hjust = -0.1, vjust = 2, color = "#27ae60", fontface = "bold") +
    labs(
      title = "FIA Model RMSE Variability Across MC Replicates",
      subtitle = "Dashed line = NEFIN (true coordinates) performance",
      x = "RMSE (Mg/ha)", y = "Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(fs::path(cfg$output_dir, "figures", "fig2_mc_variability.png"), fig2,
         width = 10, height = 6, dpi = 300)
  cat("  ✓ fig2_mc_variability.png\n")
  
  # --- Figure 3: Holdout Performance ---
  holdout_df <- data.frame(
    model = c("NEFIN-trained", rep("FIA-trained", length(fia_holdout_rmse))),
    rmse = c(rmse_nefin_holdout, fia_holdout_rmse)
  )
  
  fig3 <- ggplot(holdout_df, aes(x = model, y = rmse, fill = model)) +
    geom_boxplot(alpha = 0.7) +
    geom_point(size = 3) +
    scale_fill_manual(values = c("FIA-trained" = "#e74c3c", "NEFIN-trained" = "#27ae60")) +
    labs(
      title = "Holdout Prediction RMSE",
      subtitle = paste0("Tested on ", nrow(nefin_holdout), " NEFIN plots with true biomass"),
      x = "", y = "RMSE (Mg/ha)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")
  
  ggsave(fs::path(cfg$output_dir, "figures", "fig3_holdout_rmse.png"), fig3,
         width = 8, height = 6, dpi = 300)
  cat("  ✓ fig3_holdout_rmse.png\n")
  
  # --- Figure 4: Observed vs Predicted ---
  scatter_df <- data.frame(
    observed = y_holdout,
    predicted_nefin = preds_nefin,
    residual_nefin = y_holdout - preds_nefin
  )
  
  fig4 <- ggplot(scatter_df, aes(x = observed, y = predicted_nefin)) +
    geom_point(alpha = 0.4, color = "#27ae60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    labs(
      title = "NEFIN Model: Observed vs Predicted Biomass",
      subtitle = paste0("Holdout R² = ", round(r2_nefin_holdout, 3), 
                        " | RMSE = ", round(rmse_nefin_holdout, 1), " Mg/ha"),
      x = "Observed Biomass (Mg/ha)", y = "Predicted Biomass (Mg/ha)"
    ) +
    coord_fixed() +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(fs::path(cfg$output_dir, "figures", "fig4_observed_vs_predicted.png"), fig4,
         width = 8, height = 8, dpi = 300)
  cat("  ✓ fig4_observed_vs_predicted.png\n")
  
  # ===========================================================================
  # 6. CREATE MAPS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 6: Creating Maps\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Load hex grid and states
  hex <- if (fs::file_exists(cfg$hex_path)) st_read(cfg$hex_path, quiet = TRUE) else NULL
  states <- if (fs::file_exists(cfg$states_path)) st_read(cfg$states_path, quiet = TRUE) else NULL
  
  # --- Map 1: Residual Map (NEFIN model) ---
  if (!is.null(hex)) {
    nefin_holdout$residual <- scatter_df$residual_nefin
    nefin_holdout$lon <- nefin_complete$lon[holdout_idx]
    nefin_holdout$lat <- nefin_complete$lat[holdout_idx]
    
    holdout_sf <- st_as_sf(nefin_holdout, coords = c("lon", "lat"), crs = 4326) %>%
      st_transform(st_crs(hex))
    
    map1 <- ggplot() +
      geom_sf(data = hex, fill = "gray95", color = "gray70", linewidth = 0.1) +
      geom_sf(data = holdout_sf, aes(color = residual), size = 2, alpha = 0.7) +
      scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                            midpoint = 0, limits = c(-150, 150), oob = scales::squish) +
      labs(
        title = "NEFIN Model Residuals (Observed - Predicted)",
        subtitle = "Blue = underprediction | Red = overprediction",
        color = "Residual\n(Mg/ha)"
      ) +
      theme_void(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))
    
    if (!is.null(states)) {
      map1 <- map1 + geom_sf(data = states, fill = NA, color = "black", linewidth = 0.5)
    }
    
    ggsave(fs::path(cfg$output_dir, "maps", "map1_residuals_nefin.png"), map1,
           width = 10, height = 8, dpi = 300)
    cat("  ✓ map1_residuals_nefin.png\n")
  }
  
  # --- Map 2: MODIS NDVI with Plot Overlay ---
  ndvi_path <- cfg$rasters[[window]]$ndvi_modis
  if (fs::file_exists(ndvi_path)) {
    r_ndvi <- terra::rast(ndvi_path)
    
    # Convert to data frame for ggplot
    ndvi_df <- as.data.frame(r_ndvi, xy = TRUE)
    names(ndvi_df)[3] <- "ndvi"
    
    # Subsample for plotting
    set.seed(42)
    ndvi_sample <- ndvi_df %>% 
      filter(!is.na(ndvi)) %>%
      sample_n(min(100000, nrow(.)))
    
    map2 <- ggplot() +
      geom_raster(data = ndvi_sample, aes(x = x, y = y, fill = ndvi)) +
      scale_fill_viridis_c(option = "G", direction = -1, na.value = "white") +
      coord_sf(crs = crs(r_ndvi)) +
      labs(
        title = paste0("MODIS NDVI (", window, ") with NEFIN Plots"),
        fill = "NDVI"
      ) +
      theme_void(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))
    
    # Add plots
    holdout_reproj <- st_transform(holdout_sf, crs(r_ndvi))
    map2 <- map2 + geom_sf(data = holdout_reproj, color = "red", size = 0.5, alpha = 0.5)
    
    if (!is.null(states)) {
      states_reproj <- st_transform(states, crs(r_ndvi))
      map2 <- map2 + geom_sf(data = states_reproj, fill = NA, color = "black", linewidth = 0.3)
    }
    
    ggsave(fs::path(cfg$output_dir, "maps", "map2_modis_ndvi_with_plots.png"), map2,
           width = 12, height = 10, dpi = 300)
    cat("  ✓ map2_modis_ndvi_with_plots.png\n")
  }
  
  # ===========================================================================
  # 7. SAVE RESULTS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 7: Saving Results\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Summary table
  summary_df <- data.frame(
    metric = c("NEFIN CV RMSE", "NEFIN Holdout RMSE", "FIA Mean CV RMSE", 
               "FIA Holdout RMSE", "ΔRMSE (FIA-NEFIN)", "Improvement %"),
    value = c(
      nefin_results$overall_rmse,
      rmse_nefin_holdout,
      mean(mc_summary$rmse),
      mean(fia_holdout_rmse),
      delta_rmse,
      100 * delta_rmse / mean(fia_holdout_rmse)
    )
  )
  
  write_csv(summary_df, fs::path(cfg$output_dir, "model_comparison_summary.csv"))
  write_csv(mc_summary, fs::path(cfg$output_dir, "mc_replicate_results.csv"))
  write_csv(cv_comparison, fs::path(cfg$output_dir, "cv_fold_results.csv"))
  
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
  cat("  NEFIN (true coords) CV RMSE:     ", round(nefin_results$overall_rmse, 2), "Mg/ha\n")
  cat("  FIA (fuzzed coords) Mean CV RMSE:", round(mean(mc_summary$rmse), 2), "Mg/ha\n")
  cat("  NEFIN holdout RMSE:              ", round(rmse_nefin_holdout, 2), "Mg/ha\n")
  cat("  FIA holdout RMSE:                ", round(mean(fia_holdout_rmse), 2), "Mg/ha\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  cat("  ΔRMSE (FIA - NEFIN):             ", round(delta_rmse, 2), "Mg/ha\n")
  cat("  True coordinates improve RMSE by:", round(100 * delta_rmse / mean(fia_holdout_rmse), 1), "%\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  
  cat("\nOutput directory:", cfg$output_dir, "\n")
  
  invisible(list(
    nefin_results = nefin_results,
    mc_summary = mc_summary,
    delta_rmse = delta_rmse
  ))
}

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  run_fia_nefin_modeling()
}
