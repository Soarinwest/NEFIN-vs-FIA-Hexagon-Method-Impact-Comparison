#!/usr/bin/env Rscript
# R/spatial_model_comparison.R
# Phase 2: Does spatial fidelity affect model-based biomass prediction?
# Compare FIA (fuzzed) vs NEFIN (true) coordinates
# Author: Soren Donisvitch
# Date: December 2024

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sf)
  library(terra)
  library(mgcv)      # GAMs
  library(xgboost)   # XGBoost
  library(patchwork)
  library(viridis)
  library(fs)
})

# =============================================================================
# MAIN FUNCTION
# =============================================================================

run_spatial_model_comparison <- function(
    fia_covariates_path = "data/processed/fia_covariates.csv",
    fia_biomass_path = "data/processed/plot_hex_assignments.csv",
    nefin_path = "data/processed/nefin_processed.csv",
    nefin_climate_path = "data/processed/climate_at_plots/nefin_climate.csv",
    nefin_ndvi_path = "data/processed/ndvi_at_plots/nefin_ndvi.csv",
    output_dir = "runs/spatial_model_comparison",
    time_window = "2020_2024"  # Which 5-year window to use
) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════╗\n")
  cat("║  Phase 2: Spatial Model Comparison                                    ║\n")
  cat("║  Does spatial fidelity affect model-based biomass prediction?         ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  fs::dir_create(output_dir, recurse = TRUE)
  
  # =========================================================================
  # 1. BUILD MODELING DATASETS
  # =========================================================================
  
  cat("═══════════════════════════════════════════════════════════════════════\n")
  cat("STEP 1: Building Modeling Datasets\n")
  cat("═══════════════════════════════════════════════════════════════════════\n\n")
  
  # --- FIA Dataset (fuzzed coordinates) ---
  cat("Loading FIA data (fuzzed coordinates)...\n")
  
  fia_cov <- read_csv(fia_covariates_path, show_col_types = FALSE)
  fia_plots <- read_csv(fia_biomass_path, show_col_types = FALSE)
  
  # Get biomass from tree data (need to join with interim data)
  # For now, use a placeholder - we need to compute plot-level biomass
  # This should come from the FIA tree table aggregated to plot level
  
  # Check what columns we have
  cat("  FIA covariate columns:", paste(names(fia_cov), collapse = ", "), "\n")
  cat("  FIA plot columns:", paste(names(fia_plots)[1:10], collapse = ", "), "...\n")
  
  # Filter to time window based on MEASYEAR
  year_start <- as.integer(substr(time_window, 1, 4))
  year_end <- as.integer(substr(time_window, 6, 9))
  
  fia_cov <- fia_cov %>%
    filter(MEASYEAR >= year_start, MEASYEAR <= year_end)
  
  cat("  FIA plots in", time_window, "window:", nrow(fia_cov), "\n")
  
  # --- NEFIN Dataset (true coordinates) ---
  cat("\nLoading NEFIN data (true coordinates)...\n")
  
  nefin <- read_csv(nefin_path, show_col_types = FALSE)
  
  cat("  NEFIN plots:", nrow(nefin), "\n")
  cat("  NEFIN columns:", paste(names(nefin), collapse = ", "), "\n")
  
  # Load NEFIN covariates - these should already be at plot level
  # Don't join if it would create duplicates
  nefin_climate <- if (fs::file_exists(nefin_climate_path)) {
    clim <- read_csv(nefin_climate_path, show_col_types = FALSE)
    # Deduplicate by taking first row per plot
    if ("CN" %in% names(clim)) {
      clim <- clim %>% distinct(CN, .keep_all = TRUE)
    }
    clim
  } else NULL
  
  nefin_ndvi <- if (fs::file_exists(nefin_ndvi_path)) {
    ndvi <- read_csv(nefin_ndvi_path, show_col_types = FALSE)
    if ("CN" %in% names(ndvi)) {
      ndvi <- ndvi %>% distinct(CN, .keep_all = TRUE)
    }
    ndvi
  } else NULL
  
  # Join covariates carefully
  if (!is.null(nefin_climate) && "CN" %in% names(nefin) && "CN" %in% names(nefin_climate)) {
    # Only keep climate columns we need
    clim_cols <- intersect(names(nefin_climate), c("CN", "tmean", "ppt"))
    if (length(clim_cols) > 1) {
      nefin <- nefin %>% 
        left_join(nefin_climate %>% select(all_of(clim_cols)), by = "CN")
    }
  }
  
  if (!is.null(nefin_ndvi) && "CN" %in% names(nefin) && "CN" %in% names(nefin_ndvi)) {
    ndvi_cols <- intersect(names(nefin_ndvi), c("CN", "ndvi_modis", "ndvi_s2"))
    if (length(ndvi_cols) > 1) {
      nefin <- nefin %>% 
        left_join(nefin_ndvi %>% select(all_of(ndvi_cols)), by = "CN")
    }
  }
  
  cat("  NEFIN after covariate join:", nrow(nefin), "rows\n")
  
  # =========================================================================
  # 2. PREPARE ANALYSIS DATASETS WITH HOLDOUT
  # =========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════\n")
  cat("STEP 2: Preparing Analysis Datasets (with NEFIN holdout)\n")
  cat("═══════════════════════════════════════════════════════════════════════\n\n")
  
  # --- FIA modeling dataset ---
  # Need: biomass, ndvi, tmean, ppt, coordinates
  
  fia_model <- fia_cov %>%
    select(CN, STATECD, MEASYEAR, 
           any_of(c("tmean", "ppt", "ndvi_modis", "ndvi_s2",
                    "lon_original", "lat_original"))) %>%
    filter(!is.na(tmean), !is.na(ndvi_modis))
  
  cat("  FIA modeling dataset (with covariates):", nrow(fia_model), "plots\n")
  
  # --- NEFIN modeling dataset ---
  # Find coordinate columns
  lon_col <- names(nefin)[grepl("lon", names(nefin), ignore.case = TRUE)][1]
  lat_col <- names(nefin)[grepl("lat", names(nefin), ignore.case = TRUE)][1]
  
  nefin_model <- nefin %>%
    select(any_of(c("CN", "PLOT_ID", "aglb_Mg_per_ha", 
                    "tmean", "ppt", "ndvi_modis", "ndvi_s2",
                    lon_col, lat_col))) %>%
    filter(!is.na(aglb_Mg_per_ha))
  
  # Rename coordinate columns for consistency
  if (!is.null(lon_col) && lon_col %in% names(nefin_model)) {
    nefin_model <- nefin_model %>% rename(lon = !!lon_col, lat = !!lat_col)
  }
  
  cat("  NEFIN modeling dataset:", nrow(nefin_model), "plots\n")
  
  # --- SPLIT NEFIN: 70% train, 30% holdout for prediction testing ---
  cat("\n  Splitting NEFIN into train/holdout...\n")
  
  set.seed(42)
  nefin_with_covars <- nefin_model %>%
    filter(!is.na(tmean), !is.na(ppt), !is.na(ndvi_modis))
  
  n_nefin <- nrow(nefin_with_covars)
  holdout_idx <- sample(1:n_nefin, size = round(0.3 * n_nefin))
  
  nefin_train <- nefin_with_covars[-holdout_idx, ]
  nefin_holdout <- nefin_with_covars[holdout_idx, ]
  
  cat("    NEFIN train:", nrow(nefin_train), "plots\n")
  cat("    NEFIN holdout:", nrow(nefin_holdout), "plots (for prediction accuracy)\n")
  
  # Check covariate availability
  cat("\n  Covariate availability:\n")
  cat("    FIA        - tmean:", sum(!is.na(fia_model$tmean)), 
      ", ppt:", sum(!is.na(fia_model$ppt)),
      ", ndvi:", sum(!is.na(fia_model$ndvi_modis)), "\n")
  cat("    NEFIN train - tmean:", sum(!is.na(nefin_train$tmean)), 
      ", ppt:", sum(!is.na(nefin_train$ppt)),
      ", ndvi:", sum(!is.na(nefin_train$ndvi_modis)), "\n")
  cat("    NEFIN holdout - tmean:", sum(!is.na(nefin_holdout$tmean)), 
      ", ppt:", sum(!is.na(nefin_holdout$ppt)),
      ", ndvi:", sum(!is.na(nefin_holdout$ndvi_modis)), "\n")
  
  # =========================================================================
  # 3. COMPUTE PLOT-LEVEL BIOMASS FOR FIA
  # =========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════\n")
  cat("STEP 3: Computing FIA Plot-Level Biomass\n")
  cat("═══════════════════════════════════════════════════════════════════════\n\n")
  
  # Load tree data to compute plot-level biomass
  tree_files <- list.files("data/interim/fia/states", pattern = "tree.csv", 
                           recursive = TRUE, full.names = TRUE)
  
  if (length(tree_files) > 0) {
    cat("  Loading tree data from", length(tree_files), "state files...\n")
    
    tree_list <- lapply(tree_files, function(f) {
      read_csv(f, show_col_types = FALSE) %>%
        select(any_of(c("CN", "PLT_CN", "STATECD", "COUNTYCD", "PLOT", 
                        "INVYR", "CARBON_AG", "DRYBIO_AG", "TPA_UNADJ")))
    })
    
    trees <- bind_rows(tree_list)
    cat("  Total tree records:", nrow(trees), "\n")
    
    # Compute plot-level aboveground biomass (Mg/ha)
    # DRYBIO_AG is in pounds, TPA_UNADJ is trees per acre
    # Convert to Mg/ha: pounds * TPA * 0.0004536 * 2.471
    
    plot_biomass <- trees %>%
      filter(!is.na(DRYBIO_AG), !is.na(TPA_UNADJ)) %>%
      group_by(PLT_CN) %>%
      summarise(
        aglb_Mg_per_ha = sum(DRYBIO_AG * TPA_UNADJ * 0.0004536 * 2.471, na.rm = TRUE),
        n_trees = n(),
        .groups = "drop"
      ) %>%
      rename(CN = PLT_CN)
    
    cat("  Computed biomass for", nrow(plot_biomass), "plots\n")
    cat("  Mean biomass:", round(mean(plot_biomass$aglb_Mg_per_ha), 1), "Mg/ha\n")
    
    # Join with FIA covariates
    fia_model <- fia_model %>%
      left_join(plot_biomass, by = "CN") %>%
      filter(!is.na(aglb_Mg_per_ha))
    
    cat("  FIA plots with biomass + covariates:", nrow(fia_model), "\n")
    
  } else {
    cat("  ⚠ No tree data found. Using regional tree data...\n")
    
    tree_file <- "data/interim/fia_region/tree.csv"
    if (fs::file_exists(tree_file)) {
      trees <- read_csv(tree_file, show_col_types = FALSE)
      
      plot_biomass <- trees %>%
        filter(!is.na(DRYBIO_AG), !is.na(TPA_UNADJ)) %>%
        group_by(PLT_CN) %>%
        summarise(
          aglb_Mg_per_ha = sum(DRYBIO_AG * TPA_UNADJ * 0.0004536 * 2.471, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        rename(CN = PLT_CN)
      
      fia_model <- fia_model %>%
        left_join(plot_biomass, by = "CN") %>%
        filter(!is.na(aglb_Mg_per_ha))
      
      cat("  FIA plots with biomass + covariates:", nrow(fia_model), "\n")
    }
  }
  
  # =========================================================================
  # 4. FIT MODELS
  # =========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════\n")
  cat("STEP 4: Fitting Models\n")
  cat("═══════════════════════════════════════════════════════════════════════\n\n")
  
  results <- list()
  
  # --- 4a. Linear Models ---
  cat("4a. Linear Models (Biomass ~ NDVI + Climate)\n")
  cat("─────────────────────────────────────────────\n")
  
  # FIA model
  if (nrow(fia_model) > 50 && "aglb_Mg_per_ha" %in% names(fia_model)) {
    fia_lm <- tryCatch({
      lm(aglb_Mg_per_ha ~ ndvi_modis + tmean + ppt, data = fia_model)
    }, error = function(e) NULL)
    
    if (!is.null(fia_lm)) {
      fia_lm_summary <- summary(fia_lm)
      cat("\n  FIA Linear Model:\n")
      cat("    R² =", round(fia_lm_summary$r.squared, 4), "\n")
      cat("    RMSE =", round(sqrt(mean(fia_lm$residuals^2)), 2), "Mg/ha\n")
      cat("    NDVI coef =", round(coef(fia_lm)["ndvi_modis"], 4), "\n")
      
      results$fia_lm <- list(
        r2 = fia_lm_summary$r.squared,
        rmse = sqrt(mean(fia_lm$residuals^2)),
        coefs = coef(fia_lm),
        n = nrow(fia_model)
      )
    }
  } else {
    cat("  ⚠ Insufficient FIA data for modeling\n")
  }
  
  # NEFIN model (trained on nefin_train only, NOT holdout)
  if (nrow(nefin_train) > 50 && sum(!is.na(nefin_train$tmean)) > 30) {
    nefin_lm <- tryCatch({
      lm(aglb_Mg_per_ha ~ ndvi_modis + tmean + ppt, data = nefin_train)
    }, error = function(e) NULL)
    
    if (!is.null(nefin_lm)) {
      nefin_lm_summary <- summary(nefin_lm)
      cat("\n  NEFIN Linear Model (trained on", nrow(nefin_train), "plots):\n")
      cat("    R² =", round(nefin_lm_summary$r.squared, 4), "\n")
      cat("    RMSE =", round(sqrt(mean(nefin_lm$residuals^2)), 2), "Mg/ha\n")
      cat("    NDVI coef =", round(coef(nefin_lm)["ndvi_modis"], 4), "\n")
      
      results$nefin_lm <- list(
        r2 = nefin_lm_summary$r.squared,
        rmse = sqrt(mean(nefin_lm$residuals^2)),
        coefs = coef(nefin_lm),
        n = nrow(nefin_train %>% filter(!is.na(tmean)))
      )
    }
  } else {
    cat("  ⚠ Insufficient NEFIN data with covariates for modeling\n")
  }
  
  # --- 4b. XGBoost Models ---
  cat("\n4b. XGBoost Models\n")
  cat("─────────────────────────────────────────────\n")
  
  # FIA XGBoost
  if (exists("fia_lm") && !is.null(fia_lm)) {
    fia_xgb_data <- fia_model %>%
      filter(!is.na(aglb_Mg_per_ha), !is.na(ndvi_modis), !is.na(tmean), !is.na(ppt))
    
    if (nrow(fia_xgb_data) > 100) {
      # Prepare data for xgboost
      fia_X <- as.matrix(fia_xgb_data %>% select(ndvi_modis, tmean, ppt))
      fia_y <- fia_xgb_data$aglb_Mg_per_ha
      
      fia_dtrain <- xgb.DMatrix(data = fia_X, label = fia_y)
      
      # Cross-validation to get OOB-like estimate
      set.seed(42)
      fia_cv <- tryCatch({
        xgb.cv(
          data = fia_dtrain,
          nrounds = 200,
          nfold = 5,
          params = list(
            objective = "reg:squarederror",
            eta = 0.05,
            max_depth = 4,
            subsample = 0.8,
            colsample_bytree = 0.8,
            min_child_weight = 5
          ),
          early_stopping_rounds = 20,
          verbose = 0
        )
      }, error = function(e) NULL)
      
      if (!is.null(fia_cv) && !is.null(fia_cv$best_iteration) && fia_cv$best_iteration > 0) {
        # Fit final model
        best_nrounds <- max(fia_cv$best_iteration, 10)
        fia_xgb <- xgb.train(
          data = fia_dtrain,
          nrounds = best_nrounds,
          params = list(
            objective = "reg:squarederror",
            eta = 0.05,
            max_depth = 4,
            subsample = 0.8,
            colsample_bytree = 0.8,
            min_child_weight = 5
          ),
          verbose = 0
        )
        
        # Get predictions and metrics
        fia_preds <- predict(fia_xgb, fia_dtrain)
        fia_xgb_rmse <- sqrt(mean((fia_y - fia_preds)^2))
        fia_xgb_r2 <- 1 - sum((fia_y - fia_preds)^2) / sum((fia_y - mean(fia_y))^2)
        
        # CV RMSE (more honest estimate)
        fia_cv_rmse <- min(fia_cv$evaluation_log$test_rmse_mean)
        
        # Variable importance
        fia_importance <- xgb.importance(model = fia_xgb)
        
        cat("\n  FIA XGBoost:\n")
        cat("    R² (training) =", round(fia_xgb_r2, 4), "\n")
        cat("    RMSE (CV) =", round(fia_cv_rmse, 2), "Mg/ha\n")
        cat("    Best rounds =", best_nrounds, "\n")
        cat("    Variable importance:\n")
        for (i in 1:nrow(fia_importance)) {
          cat("      ", fia_importance$Feature[i], ":", round(fia_importance$Gain[i], 3), "\n")
        }
        
        results$fia_xgb <- list(
          r2 = fia_xgb_r2,
          rmse = fia_cv_rmse,
          importance = setNames(fia_importance$Gain, fia_importance$Feature),
          n = nrow(fia_xgb_data),
          best_rounds = best_nrounds
        )
      } else {
        cat("\n  ⚠ FIA XGBoost CV failed - using simple train/test split\n")
        
        # Fallback: simple 80/20 split
        set.seed(42)
        train_idx <- sample(1:nrow(fia_xgb_data), 0.8 * nrow(fia_xgb_data))
        
        train_X <- fia_X[train_idx, ]
        train_y <- fia_y[train_idx]
        test_X <- fia_X[-train_idx, ]
        test_y <- fia_y[-train_idx]
        
        dtrain <- xgb.DMatrix(data = train_X, label = train_y)
        dtest <- xgb.DMatrix(data = test_X, label = test_y)
        
        fia_xgb <- xgb.train(
          data = dtrain,
          nrounds = 50,
          params = list(
            objective = "reg:squarederror",
            eta = 0.1,
            max_depth = 4
          ),
          verbose = 0
        )
        
        test_preds <- predict(fia_xgb, dtest)
        fia_cv_rmse <- sqrt(mean((test_y - test_preds)^2))
        fia_xgb_r2 <- 1 - sum((test_y - test_preds)^2) / sum((test_y - mean(test_y))^2)
        
        fia_importance <- xgb.importance(model = fia_xgb)
        
        cat("    R² (test) =", round(fia_xgb_r2, 4), "\n")
        cat("    RMSE (test) =", round(fia_cv_rmse, 2), "Mg/ha\n")
        
        results$fia_xgb <- list(
          r2 = fia_xgb_r2,
          rmse = fia_cv_rmse,
          importance = setNames(fia_importance$Gain, fia_importance$Feature),
          n = nrow(fia_xgb_data),
          best_rounds = 50
        )
      }
    }
  }
  
  # NEFIN XGBoost (trained on nefin_train only)
  nefin_xgb_data <- nefin_train %>%
    filter(!is.na(aglb_Mg_per_ha), !is.na(ndvi_modis), !is.na(tmean), !is.na(ppt))
  
  cat("\n  NEFIN plots for XGBoost (train set):", nrow(nefin_xgb_data), "\n")
  
  if (nrow(nefin_xgb_data) > 100) {
    # Prepare data for xgboost
    nefin_X <- as.matrix(nefin_xgb_data %>% select(ndvi_modis, tmean, ppt))
    nefin_y <- nefin_xgb_data$aglb_Mg_per_ha
    
    nefin_dtrain <- xgb.DMatrix(data = nefin_X, label = nefin_y)
    
    # Cross-validation
    set.seed(42)
    nefin_cv <- tryCatch({
      xgb.cv(
        data = nefin_dtrain,
        nrounds = 200,
        nfold = 5,
        params = list(
          objective = "reg:squarederror",
          eta = 0.05,
          max_depth = 4,
          subsample = 0.8,
          colsample_bytree = 0.8,
          min_child_weight = 5
        ),
        early_stopping_rounds = 20,
        verbose = 0
      )
    }, error = function(e) NULL)
    
    if (!is.null(nefin_cv) && !is.null(nefin_cv$best_iteration) && nefin_cv$best_iteration > 0) {
      # Fit final model
      best_nrounds <- max(nefin_cv$best_iteration, 10)
      nefin_xgb <- xgb.train(
        data = nefin_dtrain,
        nrounds = best_nrounds,
        params = list(
          objective = "reg:squarederror",
          eta = 0.05,
          max_depth = 4,
          subsample = 0.8,
          colsample_bytree = 0.8,
          min_child_weight = 5
        ),
        verbose = 0
      )
      
      # Get predictions and metrics
      nefin_preds <- predict(nefin_xgb, nefin_dtrain)
      nefin_xgb_rmse <- sqrt(mean((nefin_y - nefin_preds)^2))
      nefin_xgb_r2 <- 1 - sum((nefin_y - nefin_preds)^2) / sum((nefin_y - mean(nefin_y))^2)
      
      # CV RMSE
      nefin_cv_rmse <- min(nefin_cv$evaluation_log$test_rmse_mean)
      
      # Variable importance
      nefin_importance <- xgb.importance(model = nefin_xgb)
      
      cat("\n  NEFIN XGBoost:\n")
      cat("    R² (training) =", round(nefin_xgb_r2, 4), "\n")
      cat("    RMSE (CV) =", round(nefin_cv_rmse, 2), "Mg/ha\n")
      cat("    Best rounds =", best_nrounds, "\n")
      cat("    Variable importance:\n")
      for (i in 1:nrow(nefin_importance)) {
        cat("      ", nefin_importance$Feature[i], ":", round(nefin_importance$Gain[i], 3), "\n")
      }
      
      results$nefin_xgb <- list(
        r2 = nefin_xgb_r2,
        rmse = nefin_cv_rmse,
        importance = setNames(nefin_importance$Gain, nefin_importance$Feature),
        n = nrow(nefin_xgb_data),
        best_rounds = best_nrounds
      )
    } else {
      cat("\n  ⚠ NEFIN XGBoost CV failed - using simple train/test split\n")
      
      # Fallback: simple 80/20 split
      set.seed(42)
      train_idx <- sample(1:nrow(nefin_xgb_data), 0.8 * nrow(nefin_xgb_data))
      
      train_X <- nefin_X[train_idx, ]
      train_y <- nefin_y[train_idx]
      test_X <- nefin_X[-train_idx, ]
      test_y <- nefin_y[-train_idx]
      
      dtrain <- xgb.DMatrix(data = train_X, label = train_y)
      dtest <- xgb.DMatrix(data = test_X, label = test_y)
      
      nefin_xgb <- xgb.train(
        data = dtrain,
        nrounds = 50,
        params = list(
          objective = "reg:squarederror",
          eta = 0.1,
          max_depth = 4
        ),
        verbose = 0
      )
      
      test_preds <- predict(nefin_xgb, dtest)
      nefin_cv_rmse <- sqrt(mean((test_y - test_preds)^2))
      nefin_xgb_r2 <- 1 - sum((test_y - test_preds)^2) / sum((test_y - mean(test_y))^2)
      
      nefin_importance <- xgb.importance(model = nefin_xgb)
      
      cat("    R² (test) =", round(nefin_xgb_r2, 4), "\n")
      cat("    RMSE (test) =", round(nefin_cv_rmse, 2), "Mg/ha\n")
      
      results$nefin_xgb <- list(
        r2 = nefin_xgb_r2,
        rmse = nefin_cv_rmse,
        importance = setNames(nefin_importance$Gain, nefin_importance$Feature),
        n = nrow(nefin_xgb_data),
        best_rounds = 50
      )
    }
  } else {
    cat("  ⚠ Insufficient NEFIN data for XGBoost\n")
  }
  
  # =========================================================================
  # 5. EVALUATE ON NEFIN HOLDOUT (TRUE BIOMASS)
  # =========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════\n")
  cat("STEP 5: Prediction Accuracy on NEFIN Holdout\n")
  cat("═══════════════════════════════════════════════════════════════════════\n\n")
  
  cat("Testing models on", nrow(nefin_holdout), "NEFIN plots with true coordinates & biomass\n\n")
  
  # Prepare holdout data
  holdout_X <- as.matrix(nefin_holdout %>% select(ndvi_modis, tmean, ppt))
  holdout_y <- nefin_holdout$aglb_Mg_per_ha
  
  holdout_results <- list()
  
  # --- Linear Model Predictions ---
  cat("5a. Linear Model Predictions on Holdout\n")
  cat("─────────────────────────────────────────────\n")
  
  if (!is.null(results$fia_lm)) {
    # FIA-trained model
    fia_lm_preds <- predict(fia_lm, newdata = nefin_holdout)
    fia_lm_holdout_rmse <- sqrt(mean((holdout_y - fia_lm_preds)^2, na.rm = TRUE))
    fia_lm_holdout_r2 <- 1 - sum((holdout_y - fia_lm_preds)^2, na.rm = TRUE) / 
                             sum((holdout_y - mean(holdout_y, na.rm = TRUE))^2, na.rm = TRUE)
    fia_lm_holdout_bias <- mean(fia_lm_preds - holdout_y, na.rm = TRUE)
    
    cat("\n  FIA-trained Linear Model → NEFIN holdout:\n")
    cat("    RMSE =", round(fia_lm_holdout_rmse, 2), "Mg/ha\n")
    cat("    R² =", round(fia_lm_holdout_r2, 4), "\n")
    cat("    Bias =", round(fia_lm_holdout_bias, 2), "Mg/ha\n")
    
    holdout_results$fia_lm <- list(
      rmse = fia_lm_holdout_rmse,
      r2 = fia_lm_holdout_r2,
      bias = fia_lm_holdout_bias
    )
  }
  
  if (!is.null(results$nefin_lm)) {
    # NEFIN-trained model (trained on nefin_train, test on nefin_holdout)
    nefin_lm_preds <- predict(nefin_lm, newdata = nefin_holdout)
    nefin_lm_holdout_rmse <- sqrt(mean((holdout_y - nefin_lm_preds)^2, na.rm = TRUE))
    nefin_lm_holdout_r2 <- 1 - sum((holdout_y - nefin_lm_preds)^2, na.rm = TRUE) / 
                               sum((holdout_y - mean(holdout_y, na.rm = TRUE))^2, na.rm = TRUE)
    nefin_lm_holdout_bias <- mean(nefin_lm_preds - holdout_y, na.rm = TRUE)
    
    cat("\n  NEFIN-trained Linear Model → NEFIN holdout:\n")
    cat("    RMSE =", round(nefin_lm_holdout_rmse, 2), "Mg/ha\n")
    cat("    R² =", round(nefin_lm_holdout_r2, 4), "\n")
    cat("    Bias =", round(nefin_lm_holdout_bias, 2), "Mg/ha\n")
    
    holdout_results$nefin_lm <- list(
      rmse = nefin_lm_holdout_rmse,
      r2 = nefin_lm_holdout_r2,
      bias = nefin_lm_holdout_bias
    )
  }
  
  # --- XGBoost Predictions ---
  cat("\n5b. XGBoost Predictions on Holdout\n")
  cat("─────────────────────────────────────────────\n")
  
  holdout_dmat <- xgb.DMatrix(data = holdout_X)
  
  if (!is.null(results$fia_xgb) && exists("fia_xgb")) {
    fia_xgb_preds <- predict(fia_xgb, holdout_dmat)
    fia_xgb_holdout_rmse <- sqrt(mean((holdout_y - fia_xgb_preds)^2, na.rm = TRUE))
    fia_xgb_holdout_r2 <- 1 - sum((holdout_y - fia_xgb_preds)^2, na.rm = TRUE) / 
                              sum((holdout_y - mean(holdout_y, na.rm = TRUE))^2, na.rm = TRUE)
    fia_xgb_holdout_bias <- mean(fia_xgb_preds - holdout_y, na.rm = TRUE)
    
    cat("\n  FIA-trained XGBoost → NEFIN holdout:\n")
    cat("    RMSE =", round(fia_xgb_holdout_rmse, 2), "Mg/ha\n")
    cat("    R² =", round(fia_xgb_holdout_r2, 4), "\n")
    cat("    Bias =", round(fia_xgb_holdout_bias, 2), "Mg/ha\n")
    
    holdout_results$fia_xgb <- list(
      rmse = fia_xgb_holdout_rmse,
      r2 = fia_xgb_holdout_r2,
      bias = fia_xgb_holdout_bias
    )
  }
  
  if (!is.null(results$nefin_xgb) && exists("nefin_xgb")) {
    nefin_xgb_preds <- predict(nefin_xgb, holdout_dmat)
    nefin_xgb_holdout_rmse <- sqrt(mean((holdout_y - nefin_xgb_preds)^2, na.rm = TRUE))
    nefin_xgb_holdout_r2 <- 1 - sum((holdout_y - nefin_xgb_preds)^2, na.rm = TRUE) / 
                                sum((holdout_y - mean(holdout_y, na.rm = TRUE))^2, na.rm = TRUE)
    nefin_xgb_holdout_bias <- mean(nefin_xgb_preds - holdout_y, na.rm = TRUE)
    
    cat("\n  NEFIN-trained XGBoost → NEFIN holdout:\n")
    cat("    RMSE =", round(nefin_xgb_holdout_rmse, 2), "Mg/ha\n")
    cat("    R² =", round(nefin_xgb_holdout_r2, 4), "\n")
    cat("    Bias =", round(nefin_xgb_holdout_bias, 2), "Mg/ha\n")
    
    holdout_results$nefin_xgb <- list(
      rmse = nefin_xgb_holdout_rmse,
      r2 = nefin_xgb_holdout_r2,
      bias = nefin_xgb_holdout_bias
    )
  }
  
  # --- Summary comparison ---
  cat("\n─────────────────────────────────────────────\n")
  cat("HOLDOUT PREDICTION SUMMARY\n")
  cat("─────────────────────────────────────────────\n\n")
  
  if (length(holdout_results) >= 2) {
    holdout_df <- data.frame(
      model = c("Linear", "Linear", "XGBoost", "XGBoost"),
      trained_on = c("FIA (fuzzed)", "NEFIN (true)", "FIA (fuzzed)", "NEFIN (true)"),
      holdout_rmse = c(
        holdout_results$fia_lm$rmse %||% NA,
        holdout_results$nefin_lm$rmse %||% NA,
        holdout_results$fia_xgb$rmse %||% NA,
        holdout_results$nefin_xgb$rmse %||% NA
      ),
      holdout_r2 = c(
        holdout_results$fia_lm$r2 %||% NA,
        holdout_results$nefin_lm$r2 %||% NA,
        holdout_results$fia_xgb$r2 %||% NA,
        holdout_results$nefin_xgb$r2 %||% NA
      ),
      holdout_bias = c(
        holdout_results$fia_lm$bias %||% NA,
        holdout_results$nefin_lm$bias %||% NA,
        holdout_results$fia_xgb$bias %||% NA,
        holdout_results$nefin_xgb$bias %||% NA
      )
    ) %>% filter(!is.na(holdout_rmse))
    
    print(holdout_df, row.names = FALSE)
    
    # Calculate improvement
    if (!is.na(holdout_results$fia_lm$rmse) && !is.na(holdout_results$nefin_lm$rmse)) {
      lm_rmse_improvement <- holdout_results$fia_lm$rmse - holdout_results$nefin_lm$rmse
      lm_rmse_pct <- 100 * lm_rmse_improvement / holdout_results$fia_lm$rmse
      cat("\n  Linear Model: NEFIN reduces holdout RMSE by", 
          round(lm_rmse_improvement, 2), "Mg/ha (", round(lm_rmse_pct, 1), "%)\n")
    }
    
    if (!is.na(holdout_results$fia_xgb$rmse) && !is.na(holdout_results$nefin_xgb$rmse)) {
      xgb_rmse_improvement <- holdout_results$fia_xgb$rmse - holdout_results$nefin_xgb$rmse
      xgb_rmse_pct <- 100 * xgb_rmse_improvement / holdout_results$fia_xgb$rmse
      cat("  XGBoost: NEFIN reduces holdout RMSE by", 
          round(xgb_rmse_improvement, 2), "Mg/ha (", round(xgb_rmse_pct, 1), "%)\n")
    }
    
    # Save holdout results
    write_csv(holdout_df, fs::path(output_dir, "holdout_prediction_results.csv"))
    cat("\n  ✓ Saved holdout_prediction_results.csv\n")
  }
  
  # Store for plotting
  results$holdout <- holdout_results
  
  # =========================================================================
  # 6. CREATE COMPARISON FIGURES
  # =========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════\n")
  cat("STEP 6: Creating Comparison Figures\n")
  cat("═══════════════════════════════════════════════════════════════════════\n\n")
  
  # --- Figure 1: Model Performance Comparison ---
  if (length(results) >= 2) {
    
    # Build comparison dataframe
    perf_data <- data.frame(
      model = c("Linear Model", "Linear Model", "XGBoost", "XGBoost"),
      dataset = c("FIA (fuzzed)", "NEFIN (true)", "FIA (fuzzed)", "NEFIN (true)"),
      r2 = c(
        results$fia_lm$r2 %||% NA,
        results$nefin_lm$r2 %||% NA,
        results$fia_xgb$r2 %||% NA,
        results$nefin_xgb$r2 %||% NA
      ),
      rmse = c(
        results$fia_lm$rmse %||% NA,
        results$nefin_lm$rmse %||% NA,
        results$fia_xgb$rmse %||% NA,
        results$nefin_xgb$rmse %||% NA
      )
    ) %>% filter(!is.na(r2))
    
    if (nrow(perf_data) > 0) {
      # R² comparison
      fig_r2 <- ggplot(perf_data, aes(x = model, y = r2, fill = dataset)) +
        geom_col(position = "dodge") +
        geom_text(aes(label = sprintf("%.3f", r2)), 
                  position = position_dodge(0.9), vjust = -0.5, size = 4) +
        scale_fill_manual(values = c("FIA (fuzzed)" = "#e74c3c", "NEFIN (true)" = "#27ae60")) +
        scale_y_continuous(limits = c(0, max(perf_data$r2, na.rm = TRUE) * 1.2)) +
        labs(
          title = "A) Model R² Comparison",
          subtitle = "Does spatial fidelity improve model fit?",
          x = "", y = "R²", fill = "Coordinates"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold"),
          legend.position = "bottom"
        )
      
      # RMSE comparison
      fig_rmse <- ggplot(perf_data, aes(x = model, y = rmse, fill = dataset)) +
        geom_col(position = "dodge") +
        geom_text(aes(label = sprintf("%.1f", rmse)), 
                  position = position_dodge(0.9), vjust = -0.5, size = 4) +
        scale_fill_manual(values = c("FIA (fuzzed)" = "#e74c3c", "NEFIN (true)" = "#27ae60")) +
        labs(
          title = "B) Model RMSE Comparison",
          subtitle = "Lower is better",
          x = "", y = "RMSE (Mg/ha)", fill = "Coordinates"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold"),
          legend.position = "bottom"
        )
      
      fig1 <- fig_r2 + fig_rmse +
        plot_annotation(
          title = "Phase 2: Spatial Fidelity Impact on Model Performance",
          subtitle = "Comparing FIA fuzzed coordinates vs NEFIN true coordinates",
          theme = theme(
            plot.title = element_text(face = "bold", size = 16),
            plot.subtitle = element_text(size = 12)
          )
        )
      
      ggsave(fs::path(output_dir, "fig1_model_performance.png"), fig1,
             width = 12, height = 6, dpi = 300)
      cat("  ✓ Saved fig1_model_performance.png\n")
    }
    
    # --- Figure 2: Coefficient Comparison (Linear Models) ---
    if (!is.null(results$fia_lm) && !is.null(results$nefin_lm)) {
      coef_data <- data.frame(
        variable = rep(c("Intercept", "NDVI", "Tmean", "PPT"), 2),
        dataset = rep(c("FIA (fuzzed)", "NEFIN (true)"), each = 4),
        coefficient = c(
          results$fia_lm$coefs,
          results$nefin_lm$coefs
        )
      ) %>%
        filter(variable != "Intercept")
      
      fig2 <- ggplot(coef_data, aes(x = variable, y = coefficient, fill = dataset)) +
        geom_col(position = "dodge") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_fill_manual(values = c("FIA (fuzzed)" = "#e74c3c", "NEFIN (true)" = "#27ae60")) +
        labs(
          title = "Linear Model Coefficients: FIA vs NEFIN",
          subtitle = "Attenuation bias would show smaller FIA coefficients",
          x = "Predictor", y = "Coefficient", fill = "Coordinates"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold"),
          legend.position = "bottom"
        )
      
      ggsave(fs::path(output_dir, "fig2_coefficient_comparison.png"), fig2,
             width = 10, height = 6, dpi = 300)
      cat("  ✓ Saved fig2_coefficient_comparison.png\n")
    }
    
    # --- Figure 3: Variable Importance (XGBoost) ---
    if (!is.null(results$fia_xgb) && !is.null(results$nefin_xgb)) {
      imp_data <- data.frame(
        variable = rep(names(results$fia_xgb$importance), 2),
        dataset = rep(c("FIA (fuzzed)", "NEFIN (true)"), each = length(results$fia_xgb$importance)),
        importance = c(
          as.numeric(results$fia_xgb$importance),
          as.numeric(results$nefin_xgb$importance)
        )
      )
      
      fig3 <- ggplot(imp_data, aes(x = reorder(variable, importance), y = importance, fill = dataset)) +
        geom_col(position = "dodge") +
        scale_fill_manual(values = c("FIA (fuzzed)" = "#e74c3c", "NEFIN (true)" = "#27ae60")) +
        coord_flip() +
        labs(
          title = "XGBoost Variable Importance (Gain)",
          subtitle = "Comparing feature importance with fuzzed vs true coordinates",
          x = "", y = "Importance (Gain)", fill = "Coordinates"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold"),
          legend.position = "bottom"
        )
      
      ggsave(fs::path(output_dir, "fig3_variable_importance.png"), fig3,
             width = 10, height = 6, dpi = 300)
      cat("  ✓ Saved fig3_variable_importance.png\n")
    }
    
    # --- Figure 4: Holdout Prediction Comparison ---
    if (length(holdout_results) >= 2) {
      holdout_plot_data <- data.frame(
        model = c("Linear", "Linear", "XGBoost", "XGBoost"),
        trained_on = c("FIA (fuzzed)", "NEFIN (true)", "FIA (fuzzed)", "NEFIN (true)"),
        rmse = c(
          holdout_results$fia_lm$rmse %||% NA,
          holdout_results$nefin_lm$rmse %||% NA,
          holdout_results$fia_xgb$rmse %||% NA,
          holdout_results$nefin_xgb$rmse %||% NA
        ),
        r2 = c(
          holdout_results$fia_lm$r2 %||% NA,
          holdout_results$nefin_lm$r2 %||% NA,
          holdout_results$fia_xgb$r2 %||% NA,
          holdout_results$nefin_xgb$r2 %||% NA
        )
      ) %>% filter(!is.na(rmse))
      
      fig4a <- ggplot(holdout_plot_data, aes(x = model, y = rmse, fill = trained_on)) +
        geom_col(position = "dodge") +
        geom_text(aes(label = sprintf("%.1f", rmse)), 
                  position = position_dodge(0.9), vjust = -0.5, size = 4) +
        scale_fill_manual(values = c("FIA (fuzzed)" = "#e74c3c", "NEFIN (true)" = "#27ae60")) +
        labs(
          title = "A) Holdout RMSE (lower is better)",
          subtitle = paste0("Tested on ", nrow(nefin_holdout), " NEFIN plots with true biomass"),
          x = "", y = "RMSE (Mg/ha)", fill = "Trained On"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold"),
          legend.position = "bottom"
        )
      
      fig4b <- ggplot(holdout_plot_data, aes(x = model, y = r2, fill = trained_on)) +
        geom_col(position = "dodge") +
        geom_text(aes(label = sprintf("%.3f", r2)), 
                  position = position_dodge(0.9), vjust = -0.5, size = 4) +
        scale_fill_manual(values = c("FIA (fuzzed)" = "#e74c3c", "NEFIN (true)" = "#27ae60")) +
        labs(
          title = "B) Holdout R² (higher is better)",
          x = "", y = "R²", fill = "Trained On"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold"),
          legend.position = "bottom"
        )
      
      fig4 <- fig4a + fig4b +
        plot_annotation(
          title = "Prediction Accuracy on Held-Out NEFIN Plots",
          subtitle = "Which training data produces better predictions at true plot locations?",
          theme = theme(
            plot.title = element_text(face = "bold", size = 16),
            plot.subtitle = element_text(size = 12)
          )
        )
      
      ggsave(fs::path(output_dir, "fig4_holdout_prediction.png"), fig4,
             width = 12, height = 6, dpi = 300)
      cat("  ✓ Saved fig4_holdout_prediction.png\n")
      
      # --- Figure 5: Observed vs Predicted scatter ---
      if (exists("nefin_xgb") && exists("fia_xgb")) {
        scatter_data <- data.frame(
          observed = rep(holdout_y, 2),
          predicted = c(
            predict(fia_xgb, holdout_dmat),
            predict(nefin_xgb, holdout_dmat)
          ),
          model = rep(c("FIA-trained", "NEFIN-trained"), each = length(holdout_y))
        )
        
        fig5 <- ggplot(scatter_data, aes(x = observed, y = predicted)) +
          geom_point(alpha = 0.3, size = 1) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
          geom_smooth(method = "lm", se = FALSE, color = "blue") +
          facet_wrap(~model) +
          labs(
            title = "XGBoost Predictions vs Observed Biomass (NEFIN Holdout)",
            subtitle = "Red dashed = 1:1 line | Blue = linear fit",
            x = "Observed Biomass (Mg/ha)", 
            y = "Predicted Biomass (Mg/ha)"
          ) +
          coord_fixed() +
          theme_minimal(base_size = 12) +
          theme(
            plot.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold", size = 12)
          )
        
        ggsave(fs::path(output_dir, "fig5_observed_vs_predicted.png"), fig5,
               width = 12, height = 6, dpi = 300)
        cat("  ✓ Saved fig5_observed_vs_predicted.png\n")
      }
    }
  }
  
  # =========================================================================
  # 7. SAVE RESULTS
  # =========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════\n")
  cat("STEP 7: Saving Results\n")
  cat("═══════════════════════════════════════════════════════════════════════\n\n")
  
  # Save summary table
  if (length(results) > 0) {
    summary_df <- data.frame(
      dataset = c("FIA", "FIA", "NEFIN", "NEFIN"),
      model = c("Linear", "XGBoost", "Linear", "XGBoost"),
      n_plots = c(
        results$fia_lm$n %||% NA,
        results$fia_xgb$n %||% NA,
        results$nefin_lm$n %||% NA,
        results$nefin_xgb$n %||% NA
      ),
      r2 = c(
        results$fia_lm$r2 %||% NA,
        results$fia_xgb$r2 %||% NA,
        results$nefin_lm$r2 %||% NA,
        results$nefin_xgb$r2 %||% NA
      ),
      rmse = c(
        results$fia_lm$rmse %||% NA,
        results$fia_xgb$rmse %||% NA,
        results$nefin_lm$rmse %||% NA,
        results$nefin_xgb$rmse %||% NA
      )
    )
    
    write_csv(summary_df, fs::path(output_dir, "model_comparison_summary.csv"))
    cat("  ✓ Saved model_comparison_summary.csv\n")
    
    # Print summary
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════════════════╗\n")
    cat("║  PHASE 2 RESULTS SUMMARY                                              ║\n")
    cat("╚══════════════════════════════════════════════════════════════════════╝\n")
    cat("\n")
    print(summary_df, row.names = FALSE)
    
    # Calculate improvements
    if (!is.na(results$fia_lm$r2) && !is.na(results$nefin_lm$r2)) {
      r2_diff_lm <- results$nefin_lm$r2 - results$fia_lm$r2
      cat("\n  Linear Model R² difference (NEFIN - FIA):", round(r2_diff_lm, 4), "\n")
    }
    
    if (!is.null(results$fia_xgb) && !is.null(results$nefin_xgb)) {
      r2_diff_xgb <- results$nefin_xgb$r2 - results$fia_xgb$r2
      cat("  XGBoost R² difference (NEFIN - FIA):", round(r2_diff_xgb, 4), "\n")
    }
  }
  
  cat("\n  Output directory:", output_dir, "\n")
  
  invisible(results)
}

# Helper function
`%||%` <- function(a, b) if (!is.null(a)) a else b

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  run_spatial_model_comparison()
}
