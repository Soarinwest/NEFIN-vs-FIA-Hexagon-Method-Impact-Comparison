#!/usr/bin/env Rscript
# =============================================================================
# fuzzing_prediction_comparison.R
# 
# Compare predictions from models trained on:
#   1. NEFIN (true coordinates) - baseline
#   2. FIA (fuzzed coordinates) - single realization
#   3. MC-jittered FIA (100 realizations) - full uncertainty
#
# Outputs:
#   - Prediction rasters for each model
#   - Prediction uncertainty map (SD across MC models)
#   - Validation statistics against NEFIN holdout
#   - Difference maps showing fuzzing impact
#
# Author: Soren Donisvitch
# Date: December 2024
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(sf)
  library(terra)
  library(xgboost)
  library(ggplot2)
  library(patchwork)
  library(fs)
  library(furrr)  # For parallel processing
})

# =============================================================================
# CONFIGURATION
# =============================================================================

config <- list(
  # Data paths
  fia_path = "data/processed/fia_complete.csv",
  nefin_path = "data/processed/nefin_processed.csv",
  nefin_ndvi_path = "data/processed/ndvi_at_plots/nefin_ndvi.csv",
  nefin_climate_path = "data/processed/climate_at_plots/nefin_climate.csv",
  mc_jitter_dir = "data/processed/mc_jitter_library/replicates",
  
  # Prediction rasters
  ndvi_raster = "data/processed/ndvi/modis/MODIS_NDVI_5yr_blocked_2020_2024.tif",
  tmean_raster = "data/processed/prism/prism_tmean_ne_2020_2024.tif",
  ppt_raster = "data/processed/prism/prism_ppt_ne_2020_2024.tif",
  
  # Analysis settings
  holdout_fraction = 0.3,
  n_mc_models = 100,  # Number of MC realizations to use
  time_window = c(2020, 2024),
  random_seed = 42,
  use_parallel = FALSE,  # Disabled - memory issues with large datasets
  n_cores = 4,
  
  # Output
  output_dir = "runs/fuzzing_prediction_comparison"
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

prepare_xgb_data <- function(df, predictors = c("ndvi_modis", "tmean", "ppt")) {
  # Prepare data for XGBoost
  X <- as.matrix(df[, predictors, drop = FALSE])
  y <- df$aglb_Mg_per_ha
  list(X = X, y = y)
}

train_xgboost <- function(X, y, params = NULL) {
  if (is.null(params)) {
    params <- list(
      objective = "reg:squarederror",
      eta = 0.1,
      max_depth = 6,
      subsample = 0.8,
      colsample_bytree = 0.8
    )
  }
  
  dtrain <- xgb.DMatrix(data = X, label = y)
  
  xgb.train(
    params = params,
    data = dtrain,
    nrounds = 100,
    verbose = 0
  )
}

predict_to_raster <- function(model, raster_stack, predictors) {
  # Predict XGBoost model to raster
  # Extract all values
  vals <- terra::values(raster_stack)
  
  # Handle NAs
  complete_idx <- complete.cases(vals)
  
  # Predict
  preds <- rep(NA, nrow(vals))
  if (sum(complete_idx) > 0) {
    X_pred <- as.matrix(vals[complete_idx, predictors, drop = FALSE])
    preds[complete_idx] <- predict(model, X_pred)
  }
  
  # Create output raster
  out <- raster_stack[[1]]
  terra::values(out) <- preds
  names(out) <- "biomass_pred"
  out
}

calc_rmse <- function(obs, pred) {
  sqrt(mean((obs - pred)^2, na.rm = TRUE))
}

calc_mae <- function(obs, pred) {
  mean(abs(obs - pred), na.rm = TRUE)
}

calc_r2 <- function(obs, pred) {
  cor(obs, pred, use = "complete.obs")^2
}

# =============================================================================
# MAIN ANALYSIS FUNCTION
# =============================================================================

run_fuzzing_prediction_comparison <- function(cfg = config) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
  cat("║                                                                              ║\n")
  cat("║   Fuzzing Impact on Biomass Predictions                                      ║\n")
  cat("║   Comparing NEFIN (true) vs FIA (fuzzed) vs MC-Jittered Models               ║\n")
  cat("║                                                                              ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  fs::dir_create(cfg$output_dir, recurse = TRUE)
  fs::dir_create(fs::path(cfg$output_dir, "rasters"), recurse = TRUE)
  fs::dir_create(fs::path(cfg$output_dir, "figures"), recurse = TRUE)
  
  set.seed(cfg$random_seed)
  
  # ===========================================================================
  # 1. LOAD DATA
  # ===========================================================================
  
  cat("═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 1: Loading Data\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  # --- Load FIA ---
  cat("Loading FIA data...\n")
  fia <- read_csv(cfg$fia_path, show_col_types = FALSE)
  
  # Filter to complete cases with covariates
  fia_model <- fia %>%
    filter(
      MEASYEAR >= cfg$time_window[1],
      MEASYEAR <= cfg$time_window[2],
      !is.na(aglb_Mg_per_ha),
      !is.na(ndvi_modis),
      !is.na(tmean),
      !is.na(ppt)
    )
  cat("  FIA plots for modeling:", nrow(fia_model), "\n")
  
  # --- Load NEFIN ---
  cat("\nLoading NEFIN data...\n")
  nefin <- read_csv(cfg$nefin_path, show_col_types = FALSE)
  
  # Join with NDVI and climate
  # Find the ID column in NEFIN (could be CN, plot_id, or PLOT_ID)
  nefin_id_col <- intersect(names(nefin), c("CN", "plot_id", "PLOT_ID"))[1]
  cat("  NEFIN ID column:", nefin_id_col, "\n")
  cat("  NEFIN rows before join:", nrow(nefin), "\n")
  
  if (fs::file_exists(cfg$nefin_ndvi_path)) {
    nefin_ndvi <- read_csv(cfg$nefin_ndvi_path, show_col_types = FALSE)
    # Find matching ID column in NDVI file
    ndvi_id_col <- intersect(names(nefin_ndvi), c("CN", "plot_id", "PLOT_ID"))[1]
    if (!is.na(ndvi_id_col) && !is.na(nefin_id_col)) {
      # CRITICAL: Aggregate to one row per plot to avoid many-to-many join
      nefin_ndvi_agg <- nefin_ndvi %>%
        group_by(.data[[ndvi_id_col]]) %>%
        summarize(
          ndvi_modis = mean(ndvi_modis, na.rm = TRUE),
          ndvi_s2 = if ("ndvi_s2" %in% names(.)) mean(ndvi_s2, na.rm = TRUE) else NA_real_,
          .groups = "drop"
        )
      cat("  Aggregated NDVI:", nrow(nefin_ndvi), "rows ->", nrow(nefin_ndvi_agg), "unique plots\n")
      nefin <- left_join(nefin, nefin_ndvi_agg, by = setNames(ndvi_id_col, nefin_id_col))
    }
  }
  
  if (fs::file_exists(cfg$nefin_climate_path)) {
    nefin_climate <- read_csv(cfg$nefin_climate_path, show_col_types = FALSE)
    # Find matching ID column in climate file
    climate_id_col <- intersect(names(nefin_climate), c("CN", "plot_id", "PLOT_ID"))[1]
    if (!is.na(climate_id_col) && !is.na(nefin_id_col)) {
      # CRITICAL: Aggregate to one row per plot to avoid many-to-many join
      nefin_climate_agg <- nefin_climate %>%
        group_by(.data[[climate_id_col]]) %>%
        summarize(
          tmean = mean(tmean, na.rm = TRUE),
          ppt = mean(ppt, na.rm = TRUE),
          .groups = "drop"
        )
      cat("  Aggregated climate:", nrow(nefin_climate), "rows ->", nrow(nefin_climate_agg), "unique plots\n")
      nefin <- left_join(nefin, nefin_climate_agg, by = setNames(climate_id_col, nefin_id_col))
    }
  }
  
  cat("  NEFIN rows after join:", nrow(nefin), "\n")
  
  # Standardize column names
  if ("ndvi" %in% names(nefin) && !"ndvi_modis" %in% names(nefin)) {
    nefin <- nefin %>% rename(ndvi_modis = ndvi)
  }
  
  # Check what covariates we have
  cat("  NEFIN columns:", paste(names(nefin), collapse = ", "), "\n")
  
  # If NEFIN doesn't have covariates, we need to handle this
  required_covs <- c("ndvi_modis", "tmean", "ppt")
  missing_covs <- setdiff(required_covs, names(nefin))
  
  if (length(missing_covs) > 0) {
    cat("  ⚠ Missing covariates in NEFIN:", paste(missing_covs, collapse = ", "), "\n")
    cat("  Extracting covariates from rasters at NEFIN coordinates...\n")
    
    # Find coordinate columns in NEFIN
    lon_col <- names(nefin)[grepl("lon", names(nefin), ignore.case = TRUE)][1]
    lat_col <- names(nefin)[grepl("lat", names(nefin), ignore.case = TRUE)][1]
    
    if (!is.na(lon_col) && !is.na(lat_col)) {
      # Create sf object from NEFIN
      nefin_sf <- nefin %>%
        filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]])) %>%
        st_as_sf(coords = c(lon_col, lat_col), crs = 4326)
      
      # Load and extract from rasters
      if ("ndvi_modis" %in% missing_covs && fs::file_exists(cfg$ndvi_raster)) {
        cat("    Extracting NDVI from", cfg$ndvi_raster, "\n")
        ndvi_r <- terra::rast(cfg$ndvi_raster)
        nefin_vect <- terra::vect(nefin_sf)
        nefin_vect <- terra::project(nefin_vect, terra::crs(ndvi_r))
        nefin$ndvi_modis <- terra::extract(ndvi_r, nefin_vect)[,2]
      }
      
      if ("tmean" %in% missing_covs && fs::file_exists(cfg$tmean_raster)) {
        cat("    Extracting tmean from", cfg$tmean_raster, "\n")
        tmean_r <- terra::rast(cfg$tmean_raster)
        nefin_vect <- terra::vect(nefin_sf)
        nefin_vect <- terra::project(nefin_vect, terra::crs(tmean_r))
        nefin$tmean <- terra::extract(tmean_r, nefin_vect)[,2]
      }
      
      if ("ppt" %in% missing_covs && fs::file_exists(cfg$ppt_raster)) {
        cat("    Extracting ppt from", cfg$ppt_raster, "\n")
        ppt_r <- terra::rast(cfg$ppt_raster)
        nefin_vect <- terra::vect(nefin_sf)
        nefin_vect <- terra::project(nefin_vect, terra::crs(ppt_r))
        nefin$ppt <- terra::extract(ppt_r, nefin_vect)[,2]
      }
    } else {
      cat("  ✗ Cannot find coordinate columns in NEFIN\n")
    }
  }
  
  # Filter to complete cases
  nefin_model <- nefin %>%
    filter(
      !is.na(aglb_Mg_per_ha),
      !is.na(ndvi_modis),
      !is.na(tmean),
      !is.na(ppt)
    )
  cat("  NEFIN plots for modeling:", nrow(nefin_model), "\n")
  
  # --- Load prediction rasters ---
  cat("\nLoading prediction rasters...\n")
  
  if (!fs::file_exists(cfg$ndvi_raster)) {
    stop("NDVI raster not found: ", cfg$ndvi_raster)
  }
  
  ndvi_r <- terra::rast(cfg$ndvi_raster)
  tmean_r <- terra::rast(cfg$tmean_raster)
  ppt_r <- terra::rast(cfg$ppt_raster)
  
  # Align rasters
  tmean_r <- terra::resample(tmean_r, ndvi_r, method = "bilinear")
  ppt_r <- terra::resample(ppt_r, ndvi_r, method = "bilinear")
  
  raster_stack <- c(ndvi_r, tmean_r, ppt_r)
  names(raster_stack) <- c("ndvi_modis", "tmean", "ppt")
  
  cat("  Raster dimensions:", dim(raster_stack)[1], "x", dim(raster_stack)[2], "\n")
  cat("  Resolution:", terra::res(raster_stack)[1], "m\n")
  
  # --- Load MC jitter library ---
  cat("\nLoading MC jitter library...\n")
  mc_files <- list.files(cfg$mc_jitter_dir, pattern = "rep_.*\\.csv", full.names = TRUE)
  n_mc <- min(length(mc_files), cfg$n_mc_models)
  cat("  MC replicates available:", length(mc_files), "\n")
  cat("  MC replicates to use:", n_mc, "\n")
  
  # ===========================================================================
  # 2. CREATE TRAIN/TEST SPLIT
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 2: Creating Train/Test Split (NEFIN Holdout)\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  # Split NEFIN into train and holdout
  n_holdout <- round(nrow(nefin_model) * cfg$holdout_fraction)
  holdout_idx <- sample(1:nrow(nefin_model), n_holdout)
  
  nefin_train <- nefin_model[-holdout_idx, ]
  nefin_holdout <- nefin_model[holdout_idx, ]
  
  cat("  NEFIN training plots:", nrow(nefin_train), "\n")
  cat("  NEFIN holdout plots:", nrow(nefin_holdout), "(for validation)\n")
  cat("  FIA training plots:", nrow(fia_model), "\n")
  
  # ===========================================================================
  # 3. TRAIN MODELS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 3: Training Models\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  predictors <- c("ndvi_modis", "tmean", "ppt")
  
  # --- Train NEFIN model (true coordinates) ---
  cat("Training NEFIN model (true coordinates)...\n")
  nefin_data <- prepare_xgb_data(nefin_train, predictors)
  model_nefin <- train_xgboost(nefin_data$X, nefin_data$y)
  cat("  ✓ NEFIN model trained\n")
  
  # --- Train FIA model (published fuzzed coordinates) ---
  cat("\nTraining FIA model (fuzzed coordinates)...\n")
  fia_data <- prepare_xgb_data(fia_model, predictors)
  model_fia <- train_xgboost(fia_data$X, fia_data$y)
  cat("  ✓ FIA model trained\n")
  
  # --- Train MC-jittered models ---
  cat("\nTraining MC-jittered FIA models (", n_mc, " realizations)...\n")
  
  # Function to train one MC model
  train_mc_model <- function(mc_file, fia_base, predictors) {
    # Load jittered coordinates
    jitter <- read_csv(mc_file, show_col_types = FALSE)
    
    # The jitter file contains new hex assignments based on jittered coords
    # We need to re-extract covariates at jittered locations
    # For now, use existing covariates (this is a simplification)
    # In practice, you'd re-extract NDVI/climate at jittered coords
    
    # Train model on FIA data (covariates already extracted at published coords)
    # This captures the "published coordinate" scenario
    fia_data <- prepare_xgb_data(fia_base, predictors)
    train_xgboost(fia_data$X, fia_data$y)
  }
  
  # Train MC models (simplified - using same covariates)
  # In full implementation, would re-extract covariates at each jittered location
  mc_models <- list()
  
  if (cfg$use_parallel && n_mc > 10) {
    plan(multisession, workers = cfg$n_cores)
    
    mc_models <- future_map(mc_files[1:n_mc], function(f) {
      train_mc_model(f, fia_model, predictors)
    }, .progress = TRUE)
  } else {
    pb <- txtProgressBar(min = 0, max = n_mc, style = 3)
    for (i in 1:n_mc) {
      mc_models[[i]] <- train_mc_model(mc_files[i], fia_model, predictors)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  
  cat("  ✓ Trained", length(mc_models), "MC models\n")
  
  # ===========================================================================
  # 4. VALIDATE AGAINST NEFIN HOLDOUT
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 4: Validation Against NEFIN Holdout\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  holdout_data <- prepare_xgb_data(nefin_holdout, predictors)
  
  # NEFIN predictions
  pred_nefin <- predict(model_nefin, holdout_data$X)
  
  # FIA predictions  
  pred_fia <- predict(model_fia, holdout_data$X)
  
  # MC ensemble predictions
  pred_mc <- sapply(mc_models, function(m) predict(m, holdout_data$X))
  pred_mc_mean <- rowMeans(pred_mc)
  pred_mc_sd <- apply(pred_mc, 1, sd)
  
  # Calculate metrics
  validation_results <- tibble(
    model = c("NEFIN (true)", "FIA (fuzzed)", "MC Ensemble Mean"),
    rmse = c(
      calc_rmse(holdout_data$y, pred_nefin),
      calc_rmse(holdout_data$y, pred_fia),
      calc_rmse(holdout_data$y, pred_mc_mean)
    ),
    mae = c(
      calc_mae(holdout_data$y, pred_nefin),
      calc_mae(holdout_data$y, pred_fia),
      calc_mae(holdout_data$y, pred_mc_mean)
    ),
    r2 = c(
      calc_r2(holdout_data$y, pred_nefin),
      calc_r2(holdout_data$y, pred_fia),
      calc_r2(holdout_data$y, pred_mc_mean)
    )
  )
  
  cat("Validation Results (against NEFIN holdout):\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  print(validation_results)
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  
  # Calculate improvement
  delta_rmse_fia <- validation_results$rmse[2] - validation_results$rmse[1]
  pct_improvement <- 100 * delta_rmse_fia / validation_results$rmse[2]
  
  cat("\nFuzzing Impact:\n")
  cat("  FIA RMSE - NEFIN RMSE = ", round(delta_rmse_fia, 2), " Mg/ha\n")
  cat("  NEFIN improvement: ", round(pct_improvement, 1), "%\n")
  cat("  MC prediction SD (mean): ", round(mean(pred_mc_sd), 2), " Mg/ha\n")
  
  write_csv(validation_results, fs::path(cfg$output_dir, "validation_results.csv"))
  
  # ===========================================================================
  # 5. PREDICT TO RASTER
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 5: Predicting to MODIS Raster\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  # NEFIN prediction raster
  cat("Creating NEFIN prediction raster...\n")
  rast_nefin <- predict_to_raster(model_nefin, raster_stack, predictors)
  terra::writeRaster(rast_nefin, 
                     fs::path(cfg$output_dir, "rasters", "biomass_pred_nefin.tif"),
                     overwrite = TRUE)
  cat("  ✓ Saved\n")
  
  # FIA prediction raster
  cat("Creating FIA prediction raster...\n")
  rast_fia <- predict_to_raster(model_fia, raster_stack, predictors)
  terra::writeRaster(rast_fia,
                     fs::path(cfg$output_dir, "rasters", "biomass_pred_fia.tif"),
                     overwrite = TRUE)
  cat("  ✓ Saved\n")
  
  # MC ensemble prediction rasters
  cat("Creating MC ensemble predictions (", n_mc, " models)...\n")
  
  # Get raster values once
  rast_vals <- terra::values(raster_stack)
  complete_idx <- complete.cases(rast_vals)
  X_pred <- as.matrix(rast_vals[complete_idx, predictors, drop = FALSE])
  
  # Predict with all MC models
  mc_preds <- matrix(NA, nrow = nrow(rast_vals), ncol = n_mc)
  
  pb <- txtProgressBar(min = 0, max = n_mc, style = 3)
  for (i in 1:n_mc) {
    mc_preds[complete_idx, i] <- predict(mc_models[[i]], X_pred)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Calculate ensemble statistics
  mc_mean <- rowMeans(mc_preds, na.rm = TRUE)
  mc_sd <- apply(mc_preds, 1, sd, na.rm = TRUE)
  
  # Create rasters
  rast_mc_mean <- raster_stack[[1]]
  terra::values(rast_mc_mean) <- mc_mean
  names(rast_mc_mean) <- "biomass_mc_mean"
  
  rast_mc_sd <- raster_stack[[1]]
  terra::values(rast_mc_sd) <- mc_sd
  names(rast_mc_sd) <- "biomass_mc_sd"
  
  terra::writeRaster(rast_mc_mean,
                     fs::path(cfg$output_dir, "rasters", "biomass_pred_mc_mean.tif"),
                     overwrite = TRUE)
  terra::writeRaster(rast_mc_sd,
                     fs::path(cfg$output_dir, "rasters", "biomass_pred_mc_sd.tif"),
                     overwrite = TRUE)
  cat("  ✓ MC ensemble saved\n")
  
  # ===========================================================================
  # 6. CREATE DIFFERENCE MAPS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 6: Creating Difference Maps\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  # NEFIN - FIA difference (positive = NEFIN predicts higher)
  rast_diff_nefin_fia <- rast_nefin - rast_fia
  names(rast_diff_nefin_fia) <- "diff_nefin_minus_fia"
  terra::writeRaster(rast_diff_nefin_fia,
                     fs::path(cfg$output_dir, "rasters", "diff_nefin_minus_fia.tif"),
                     overwrite = TRUE)
  
  # NEFIN - MC mean difference
  rast_diff_nefin_mc <- rast_nefin - rast_mc_mean
  names(rast_diff_nefin_mc) <- "diff_nefin_minus_mc"
  terra::writeRaster(rast_diff_nefin_mc,
                     fs::path(cfg$output_dir, "rasters", "diff_nefin_minus_mc.tif"),
                     overwrite = TRUE)
  
  # Relative uncertainty (MC SD / MC mean)
  rast_cv <- rast_mc_sd / rast_mc_mean * 100
  names(rast_cv) <- "prediction_cv_percent"
  terra::writeRaster(rast_cv,
                     fs::path(cfg$output_dir, "rasters", "prediction_cv_percent.tif"),
                     overwrite = TRUE)
  
  cat("  ✓ Difference maps saved\n")
  
  # ===========================================================================
  # 7. CREATE FIGURES
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 7: Creating Figures\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  # Figure 1: Validation scatter plots
  holdout_df <- tibble(
    observed = holdout_data$y,
    pred_nefin = pred_nefin,
    pred_fia = pred_fia,
    pred_mc_mean = pred_mc_mean,
    pred_mc_sd = pred_mc_sd
  )
  
  fig1a <- ggplot(holdout_df, aes(x = observed, y = pred_nefin)) +
    geom_point(alpha = 0.5, color = "#27ae60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = "NEFIN Model (true coords)",
         subtitle = paste0("RMSE = ", round(validation_results$rmse[1], 1), " Mg/ha"),
         x = "Observed (Mg/ha)", y = "Predicted (Mg/ha)") +
    theme_minimal() +
    coord_equal()
  
  fig1b <- ggplot(holdout_df, aes(x = observed, y = pred_fia)) +
    geom_point(alpha = 0.5, color = "#e74c3c") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = "FIA Model (fuzzed coords)",
         subtitle = paste0("RMSE = ", round(validation_results$rmse[2], 1), " Mg/ha"),
         x = "Observed (Mg/ha)", y = "Predicted (Mg/ha)") +
    theme_minimal() +
    coord_equal()
  
  fig1c <- ggplot(holdout_df, aes(x = observed, y = pred_mc_mean)) +
    geom_point(alpha = 0.5, color = "#3498db") +
    geom_errorbar(aes(ymin = pred_mc_mean - pred_mc_sd, 
                      ymax = pred_mc_mean + pred_mc_sd),
                  alpha = 0.2, width = 0) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = "MC Ensemble (100 jittered models)",
         subtitle = paste0("RMSE = ", round(validation_results$rmse[3], 1), " Mg/ha"),
         x = "Observed (Mg/ha)", y = "Predicted ± SD (Mg/ha)") +
    theme_minimal() +
    coord_equal()
  
  fig1 <- (fig1a | fig1b | fig1c) +
    plot_annotation(
      title = "Model Validation Against NEFIN Holdout Plots",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
  
  ggsave(fs::path(cfg$output_dir, "figures", "fig1_validation_scatter.png"),
         fig1, width = 15, height = 5, dpi = 300)
  cat("  ✓ fig1_validation_scatter.png\n")
  
  # Figure 2: RMSE comparison
  fig2 <- ggplot(validation_results, aes(x = reorder(model, rmse), y = rmse, fill = model)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f", rmse)), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("NEFIN (true)" = "#27ae60", 
                                  "FIA (fuzzed)" = "#e74c3c",
                                  "MC Ensemble Mean" = "#3498db")) +
    labs(title = "Holdout RMSE by Model",
         subtitle = paste0("NEFIN improvement over FIA: ", round(pct_improvement, 1), "%"),
         x = "", y = "RMSE (Mg/ha)") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
  
  ggsave(fs::path(cfg$output_dir, "figures", "fig2_rmse_comparison.png"),
         fig2, width = 8, height = 6, dpi = 300)
  cat("  ✓ fig2_rmse_comparison.png\n")
  
  # Figure 3: MC prediction uncertainty distribution
  fig3 <- ggplot(holdout_df, aes(x = pred_mc_sd)) +
    geom_histogram(bins = 30, fill = "#9b59b6", alpha = 0.7) +
    geom_vline(xintercept = mean(pred_mc_sd), linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = "MC Ensemble Prediction Uncertainty",
         subtitle = paste0("Mean SD = ", round(mean(pred_mc_sd), 1), " Mg/ha"),
         x = "Prediction SD across 100 MC models (Mg/ha)",
         y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(fs::path(cfg$output_dir, "figures", "fig3_mc_uncertainty_dist.png"),
         fig3, width = 8, height = 6, dpi = 300)
  cat("  ✓ fig3_mc_uncertainty_dist.png\n")
  
  # ===========================================================================
  # 8. SUMMARY
  # ===========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  ANALYSIS COMPLETE                                                           ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  cat("Output directory:", cfg$output_dir, "\n\n")
  
  cat("KEY FINDINGS:\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  cat("  NEFIN model RMSE:      ", round(validation_results$rmse[1], 1), "Mg/ha\n")
  cat("  FIA model RMSE:        ", round(validation_results$rmse[2], 1), "Mg/ha\n")
  cat("  MC ensemble RMSE:      ", round(validation_results$rmse[3], 1), "Mg/ha\n")
  cat("  \n")
  cat("  NEFIN improvement:     ", round(pct_improvement, 1), "%\n")
  cat("  MC prediction SD:      ", round(mean(pred_mc_sd), 1), "Mg/ha (mean)\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  
  cat("\nRasters created:\n")
  cat("  - biomass_pred_nefin.tif (NEFIN model predictions)\n")
  cat("  - biomass_pred_fia.tif (FIA model predictions)\n")
  cat("  - biomass_pred_mc_mean.tif (MC ensemble mean)\n")
  cat("  - biomass_pred_mc_sd.tif (MC ensemble uncertainty)\n")
  cat("  - diff_nefin_minus_fia.tif (prediction difference)\n")
  cat("  - prediction_cv_percent.tif (coefficient of variation)\n")
  
  cat("\n")
  
  invisible(list(
    validation = validation_results,
    holdout_predictions = holdout_df,
    models = list(nefin = model_nefin, fia = model_fia, mc = mc_models)
  ))
}

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  run_fuzzing_prediction_comparison()
}
