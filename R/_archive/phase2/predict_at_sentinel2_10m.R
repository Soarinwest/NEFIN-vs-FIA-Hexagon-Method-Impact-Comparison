#!/usr/bin/env Rscript
# =============================================================================
# predict_at_sentinel2_10m_v3.R
# 
# Compare prediction uncertainty at Sentinel-2 10m vs MODIS 250m
# v3: Fixed column detection
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(readr)
  library(dplyr)
  library(xgboost)
  library(ggplot2)
  library(patchwork)
  library(fs)
})
gc()

# =============================================================================
# CONFIGURATION
# =============================================================================

config <- list(
  fia_path = "data/processed/fia_complete.csv",
  nefin_path = "data/processed/nefin_processed.csv",
  s2_ndvi_path = "data/processed/ndvi/s2/S2_NDVI_10m_2020_2025.tif",
  tmean_path = "data/processed/prism/prism_tmean_ne_2020_2024.tif",
  ppt_path = "data/processed/prism/prism_ppt_ne_2020_2024.tif",
  modis_mc_sd_path = "runs/fuzzing_prediction_comparison/rasters/biomass_pred_mc_sd.tif",
  time_window = c(2020, 2024),
  holdout_fraction = 0.3,
  n_mc_models = 20,
  random_seed = 42,
  max_cells = 5e6,
  output_dir = "runs/s2_prediction_comparison"
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

train_xgboost <- function(X, y, seed = 42) {
  if (nrow(X) == 0) stop("No training data!")
  params <- list(
    objective = "reg:squarederror",
    eta = 0.1,
    max_depth = 6,
    subsample = 0.8,
    colsample_bytree = 0.8
  )
  set.seed(seed)
  dtrain <- xgb.DMatrix(data = X, label = y)
  xgb.train(params = params, data = dtrain, nrounds = 100, verbose = 0)
}

extract_covariates_to_points <- function(pts_sf, s2_ndvi, tmean, ppt) {
  ndvi_vals <- terra::extract(s2_ndvi, pts_sf, ID = FALSE)[[1]]
  tmean_vals <- terra::extract(tmean, pts_sf, ID = FALSE)[[1]]
  ppt_vals <- terra::extract(ppt, pts_sf, ID = FALSE)[[1]]
  data.frame(ndvi_s2 = ndvi_vals, tmean = tmean_vals, ppt = ppt_vals)
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

predict_at_s2_resolution <- function(cfg = config) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
  cat("║   Sentinel-2 10m Prediction Comparison (v3)                                  ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  fs::dir_create(cfg$output_dir, recurse = TRUE)
  fs::dir_create(fs::path(cfg$output_dir, "rasters"), recurse = TRUE)
  fs::dir_create(fs::path(cfg$output_dir, "figures"), recurse = TRUE)
  
  set.seed(cfg$random_seed)
  
  # ===========================================================================
  # 1. LOAD RASTERS
  # ===========================================================================
  
  cat("═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 1: Loading Rasters\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  s2_ndvi_full <- rast(cfg$s2_ndvi_path)
  cat("  Sentinel-2 NDVI:", nrow(s2_ndvi_full), "x", ncol(s2_ndvi_full), "\n")
  
  tmean <- rast(cfg$tmean_path)
  ppt <- rast(cfg$ppt_path)
  cat("  Climate rasters loaded\n")
  
  # ===========================================================================
  # 2. LOAD TRAINING DATA
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 2: Loading Training Data\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  # Load FIA
  fia_raw <- read_csv(cfg$fia_path, show_col_types = FALSE)
  cat("  FIA columns available:\n")
  cat("   ", paste(names(fia_raw)[1:20], collapse = ", "), "...\n")
  
  # EXPLICIT column selection - look for what we know exists
  if ("lat_original" %in% names(fia_raw)) {
    cat("  Found: lat_original, lon_original\n")
    
    # Get year column
    if ("MEASYEAR" %in% names(fia_raw)) {
      year_vals <- fia_raw$MEASYEAR
    } else if ("MEASYEAR.x" %in% names(fia_raw)) {
      year_vals <- fia_raw$MEASYEAR.x
    } else {
      year_vals <- rep(2022, nrow(fia_raw))
    }
    
    fia <- data.frame(
      lat = fia_raw$lat_original,
      lon = fia_raw$lon_original,
      year = year_vals,
      aglb_Mg_per_ha = fia_raw$aglb_Mg_per_ha
    )
  } else if ("LAT" %in% names(fia_raw)) {
    cat("  Found: LAT, LON\n")
    fia <- data.frame(
      lat = fia_raw$LAT,
      lon = fia_raw$LON,
      year = fia_raw$MEASYEAR,
      aglb_Mg_per_ha = fia_raw$aglb_Mg_per_ha
    )
  } else {
    stop("Cannot find coordinate columns. Available: ", paste(names(fia_raw), collapse = ", "))
  }
  
  fia <- fia %>%
    filter(
      year >= cfg$time_window[1],
      year <= cfg$time_window[2],
      !is.na(aglb_Mg_per_ha),
      !is.na(lat),
      !is.na(lon)
    )
  cat("  FIA plots (filtered):", nrow(fia), "\n")
  
  # Load NEFIN  
  nefin_raw <- read_csv(cfg$nefin_path, show_col_types = FALSE)
  cat("\n  NEFIN columns available:\n")
  cat("   ", paste(names(nefin_raw)[1:15], collapse = ", "), "...\n")
  
  if ("lat_public" %in% names(nefin_raw)) {
    cat("  Found: lat_public, lon_public\n")
    nefin <- data.frame(
      lat = nefin_raw$lat_public,
      lon = nefin_raw$lon_public,
      aglb_Mg_per_ha = nefin_raw$aglb_Mg_per_ha
    )
  } else if ("lat_original" %in% names(nefin_raw)) {
    cat("  Found: lat_original, lon_original\n")
    nefin <- data.frame(
      lat = nefin_raw$lat_original,
      lon = nefin_raw$lon_original,
      aglb_Mg_per_ha = nefin_raw$aglb_Mg_per_ha
    )
  } else {
    stop("Cannot find NEFIN coordinate columns")
  }
  
  nefin <- nefin %>%
    filter(!is.na(aglb_Mg_per_ha), !is.na(lat), !is.na(lon))
  cat("  NEFIN plots (filtered):", nrow(nefin), "\n")
  
  # ===========================================================================
  # 3. EXTRACT COVARIATES FROM RASTERS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 3: Extracting Covariates from Rasters\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  fia_sf <- st_as_sf(fia, coords = c("lon", "lat"), crs = 4326)
  nefin_sf <- st_as_sf(nefin, coords = c("lon", "lat"), crs = 4326)
  
  cat("  Extracting for FIA plots...\n")
  fia_covs <- extract_covariates_to_points(fia_sf, s2_ndvi_full, tmean, ppt)
  fia$ndvi_s2 <- fia_covs$ndvi_s2
  fia$tmean <- fia_covs$tmean
  fia$ppt <- fia_covs$ppt
  
  cat("  Extracting for NEFIN plots...\n")
  nefin_covs <- extract_covariates_to_points(nefin_sf, s2_ndvi_full, tmean, ppt)
  nefin$ndvi_s2 <- nefin_covs$ndvi_s2
  nefin$tmean <- nefin_covs$tmean
  nefin$ppt <- nefin_covs$ppt
  
  # Filter complete cases
  fia_model <- fia %>% filter(!is.na(ndvi_s2), !is.na(tmean), !is.na(ppt))
  nefin_model <- nefin %>% filter(!is.na(ndvi_s2), !is.na(tmean), !is.na(ppt))
  
  cat("\n  After extraction:\n")
  cat("    FIA with covariates:", nrow(fia_model), "\n")
  cat("    NEFIN with covariates:", nrow(nefin_model), "\n")
  
  # ===========================================================================
  # 4. CREATE PREDICTION SUBSET
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 4: Creating Prediction Subset\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  # Subset to manageable size
  center_lon <- -72.5
  center_lat <- 44.0
  half_size <- 0.25
  
  subset_ext <- ext(
    center_lon - half_size,
    center_lon + half_size,
    center_lat - half_size,
    center_lat + half_size
  )
  
  s2_ndvi <- crop(s2_ndvi_full, subset_ext)
  cat("  Cropped to:", nrow(s2_ndvi), "x", ncol(s2_ndvi), "cells\n")
  
  # Create prediction stack
  tmean_pred <- crop(tmean, ext(s2_ndvi))
  ppt_pred <- crop(ppt, ext(s2_ndvi))
  tmean_pred <- resample(tmean_pred, s2_ndvi, method = "bilinear")
  ppt_pred <- resample(ppt_pred, s2_ndvi, method = "bilinear")
  
  s2_stack <- c(s2_ndvi, tmean_pred, ppt_pred)
  names(s2_stack) <- c("ndvi_s2", "tmean", "ppt")
  cat("  Prediction stack ready\n")
  
  # ===========================================================================
  # 5. TRAIN MODELS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 5: Training Models\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  predictors <- c("ndvi_s2", "tmean", "ppt")
  
  # Split NEFIN
  n_holdout <- round(nrow(nefin_model) * cfg$holdout_fraction)
  holdout_idx <- sample(1:nrow(nefin_model), n_holdout)
  nefin_train <- nefin_model[-holdout_idx, ]
  nefin_holdout <- nefin_model[holdout_idx, ]
  
  cat("  NEFIN train:", nrow(nefin_train), "| holdout:", nrow(nefin_holdout), "\n")
  cat("  FIA train:", nrow(fia_model), "\n")
  
  X_nefin <- as.matrix(nefin_train[, predictors])
  y_nefin <- nefin_train$aglb_Mg_per_ha
  X_fia <- as.matrix(fia_model[, predictors])
  y_fia <- fia_model$aglb_Mg_per_ha
  
  cat("\n  Training NEFIN model...\n")
  model_nefin <- train_xgboost(X_nefin, y_nefin, seed = cfg$random_seed)
  
  cat("  Training FIA model...\n")
  model_fia <- train_xgboost(X_fia, y_fia, seed = cfg$random_seed)
  
  cat("  Training", cfg$n_mc_models, "MC models...\n")
  mc_models <- list()
  pb <- txtProgressBar(min = 0, max = cfg$n_mc_models, style = 3)
  for (i in 1:cfg$n_mc_models) {
    mc_models[[i]] <- train_xgboost(X_fia, y_fia, seed = cfg$random_seed + i)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # ===========================================================================
  # 6. VALIDATE
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 6: Validation\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  X_holdout <- as.matrix(nefin_holdout[, predictors])
  y_holdout <- nefin_holdout$aglb_Mg_per_ha
  
  pred_nefin <- predict(model_nefin, X_holdout)
  pred_fia <- predict(model_fia, X_holdout)
  
  rmse_nefin <- sqrt(mean((y_holdout - pred_nefin)^2))
  rmse_fia <- sqrt(mean((y_holdout - pred_fia)^2))
  
  cat("  NEFIN RMSE:", round(rmse_nefin, 1), "Mg/ha\n")
  cat("  FIA RMSE:  ", round(rmse_fia, 1), "Mg/ha\n")
  cat("  Improvement:", round(100 * (rmse_fia - rmse_nefin) / rmse_fia, 1), "%\n")
  
  # ===========================================================================
  # 7. PREDICT TO RASTER
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 7: Predicting to S2 Raster\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  rast_vals <- values(s2_stack)
  complete_idx <- complete.cases(rast_vals)
  n_valid <- sum(complete_idx)
  cat("  Valid cells:", format(n_valid, big.mark = ","), "\n")
  
  X_pred <- as.matrix(rast_vals[complete_idx, predictors, drop = FALSE])
  
  # NEFIN prediction
  cat("  Predicting NEFIN...\n")
  preds_nefin <- rep(NA, nrow(rast_vals))
  preds_nefin[complete_idx] <- predict(model_nefin, X_pred)
  rast_nefin <- s2_stack[[1]]
  values(rast_nefin) <- preds_nefin
  writeRaster(rast_nefin, fs::path(cfg$output_dir, "rasters", "biomass_s2_nefin.tif"), overwrite = TRUE)
  
  # FIA prediction
  cat("  Predicting FIA...\n")
  preds_fia <- rep(NA, nrow(rast_vals))
  preds_fia[complete_idx] <- predict(model_fia, X_pred)
  rast_fia <- s2_stack[[1]]
  values(rast_fia) <- preds_fia
  writeRaster(rast_fia, fs::path(cfg$output_dir, "rasters", "biomass_s2_fia.tif"), overwrite = TRUE)
  
  # MC predictions
  cat("  Predicting MC ensemble...\n")
  mc_preds <- matrix(NA, nrow = nrow(rast_vals), ncol = cfg$n_mc_models)
  pb <- txtProgressBar(min = 0, max = cfg$n_mc_models, style = 3)
  for (i in 1:cfg$n_mc_models) {
    mc_preds[complete_idx, i] <- predict(mc_models[[i]], X_pred)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  mc_mean <- rowMeans(mc_preds, na.rm = TRUE)
  mc_sd <- apply(mc_preds, 1, sd, na.rm = TRUE)
  
  rast_mc_sd <- s2_stack[[1]]
  values(rast_mc_sd) <- mc_sd
  writeRaster(rast_mc_sd, fs::path(cfg$output_dir, "rasters", "biomass_s2_mc_sd.tif"), overwrite = TRUE)
  
  # ===========================================================================
  # 8. COMPARE TO MODIS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 8: Comparing to MODIS\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  s2_sd_vals <- mc_sd[!is.na(mc_sd)]
  s2_mean_sd <- mean(s2_sd_vals)
  s2_max_sd <- max(s2_sd_vals)
  
  cat("  Sentinel-2 (10m) uncertainty:\n")
  cat("    Mean SD:", round(s2_mean_sd, 1), "Mg/ha\n")
  cat("    Max SD: ", round(s2_max_sd, 1), "Mg/ha\n")
  
  if (fs::file_exists(cfg$modis_mc_sd_path)) {
    modis_mc_sd <- rast(cfg$modis_mc_sd_path)
    modis_crop <- crop(modis_mc_sd, ext(rast_mc_sd))
    modis_vals <- values(modis_crop)
    modis_vals <- modis_vals[!is.na(modis_vals)]
    
    if (length(modis_vals) > 0) {
      modis_mean_sd <- mean(modis_vals)
      cat("\n  MODIS (250m) uncertainty:\n")
      cat("    Mean SD:", round(modis_mean_sd, 1), "Mg/ha\n")
      cat("\n  ══════════════════════════════════════════════════════════════════\n")
      cat("  S2/MODIS ratio:", round(s2_mean_sd / modis_mean_sd, 2), "×\n")
      cat("  ══════════════════════════════════════════════════════════════════\n")
    }
  }
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  COMPLETE                                                                    ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")
  cat("\nOutput:", cfg$output_dir, "\n")
}

# Run
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  predict_at_s2_resolution()
}
