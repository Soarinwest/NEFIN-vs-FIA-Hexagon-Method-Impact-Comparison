#!/usr/bin/env Rscript
# R/spatial_model_comparison_v2.R
# Phase 2: Spatial Model Comparison - UPDATED to use fia_complete.csv
# Does spatial fidelity affect model-based biomass prediction?
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
  library(fs)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

config <- list(
  # Data paths - UPDATED to use fia_complete.csv
  fia_path = "data/processed/fia_complete.csv",
  nefin_path = "data/processed/nefin_processed.csv",
  nefin_climate_path = "data/processed/climate_at_plots/nefin_climate.csv",
  nefin_ndvi_path = "data/processed/ndvi_at_plots/nefin_ndvi.csv",
  
  # Analysis settings
  time_window = c(2020, 2024),
  holdout_fraction = 0.3,
  cv_folds = 5,
  
  # Output
  output_dir = "runs/spatial_model_comparison"
)

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

run_model_comparison <- function(cfg = config) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  Phase 2: Spatial Model Comparison (v2 - with FIA biomass)               ║\n")
  cat("║  Does spatial fidelity affect model-based biomass prediction?            ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  fs::dir_create(cfg$output_dir, recurse = TRUE)
  
  # ===========================================================================
  # 1. LOAD DATA
  # ===========================================================================
  
  cat("═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 1: Loading Data\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # --- Load FIA (fuzzed coordinates) ---
  cat("Loading FIA data (fuzzed coordinates)...\n")
  
  if (!fs::file_exists(cfg$fia_path)) {
    stop("FIA complete dataset not found: ", cfg$fia_path, 
         "\n  Run: Rscript R/build_fia_biomass.R first")
  }
  
  fia <- read_csv(cfg$fia_path, show_col_types = FALSE)
  cat("  FIA plots loaded:", nrow(fia), "\n")
  
  # Check for required columns
  fia_cols <- c("CN", "aglb_Mg_per_ha", "ndvi_modis", "tmean", "ppt", "MEASYEAR")
  missing <- setdiff(fia_cols, names(fia))
  if (length(missing) > 0) {
    cat("  ⚠ Missing columns:", paste(missing, collapse = ", "), "\n")
  }
  
  # Filter to time window and complete cases
  fia_model <- fia %>%
    filter(
      MEASYEAR >= cfg$time_window[1], 
      MEASYEAR <= cfg$time_window[2],
      !is.na(aglb_Mg_per_ha),
      !is.na(ndvi_modis),
      !is.na(tmean),
      !is.na(ppt)
    )
  
  cat("  FIA plots for modeling (", cfg$time_window[1], "-", cfg$time_window[2], 
      " with biomass + covariates):", nrow(fia_model), "\n")
  
  # --- Load NEFIN (true coordinates) ---
  cat("\nLoading NEFIN data (true coordinates)...\n")
  
  nefin <- read_csv(cfg$nefin_path, show_col_types = FALSE)
  cat("  NEFIN plots loaded:", nrow(nefin), "\n")
  
  # Load NEFIN covariates
  if (fs::file_exists(cfg$nefin_climate_path)) {
    nefin_climate <- read_csv(cfg$nefin_climate_path, show_col_types = FALSE) %>%
      distinct(CN, .keep_all = TRUE)
    nefin <- nefin %>% left_join(nefin_climate %>% select(CN, tmean, ppt), by = "CN")
  }
  
  if (fs::file_exists(cfg$nefin_ndvi_path)) {
    nefin_ndvi <- read_csv(cfg$nefin_ndvi_path, show_col_types = FALSE) %>%
      distinct(CN, .keep_all = TRUE)
    nefin <- nefin %>% left_join(nefin_ndvi %>% select(CN, ndvi_modis, ndvi_s2), by = "CN")
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
  
  # ===========================================================================
  # 2. PREPARE TRAIN/HOLDOUT SPLITS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 2: Preparing Train/Holdout Splits\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  set.seed(42)
  
  # NEFIN split (holdout for final evaluation)
  n_nefin <- nrow(nefin_model)
  holdout_idx <- sample(1:n_nefin, size = round(cfg$holdout_fraction * n_nefin))
  
  nefin_train <- nefin_model[-holdout_idx, ]
  nefin_holdout <- nefin_model[holdout_idx, ]
  
  cat("  NEFIN train:", nrow(nefin_train), "\n")
  cat("  NEFIN holdout:", nrow(nefin_holdout), "(for final evaluation)\n")
  cat("  FIA train:", nrow(fia_model), "\n")
  
  # ===========================================================================
  # 3. FIT MODELS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 3: Fitting Models\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # --- 3a. Linear Models ---
  cat("3a. Linear Models (Biomass ~ NDVI + Climate)\n")
  cat("─────────────────────────────────────────────\n\n")
  
  # FIA Linear Model
  lm_fia <- lm(aglb_Mg_per_ha ~ ndvi_modis + tmean + ppt, data = fia_model)
  fia_lm_r2 <- summary(lm_fia)$r.squared
  fia_lm_rmse <- sqrt(mean(residuals(lm_fia)^2))
  
  cat("  FIA Linear Model:\n")
  cat("    R² =", round(fia_lm_r2, 4), "\n")
  cat("    RMSE =", round(fia_lm_rmse, 2), "Mg/ha\n")
  cat("    NDVI coef =", round(coef(lm_fia)["ndvi_modis"], 2), "\n\n")
  
  # NEFIN Linear Model
  lm_nefin <- lm(aglb_Mg_per_ha ~ ndvi_modis + tmean + ppt, data = nefin_train)
  nefin_lm_r2 <- summary(lm_nefin)$r.squared
  nefin_lm_rmse <- sqrt(mean(residuals(lm_nefin)^2))
  
  cat("  NEFIN Linear Model:\n")
  cat("    R² =", round(nefin_lm_r2, 4), "\n")
  cat("    RMSE =", round(nefin_lm_rmse, 2), "Mg/ha\n")
  cat("    NDVI coef =", round(coef(lm_nefin)["ndvi_modis"], 2), "\n\n")
  
  # --- 3b. XGBoost Models ---
  cat("3b. XGBoost Models\n")
  cat("─────────────────────────────────────────────\n\n")
  
  xgb_params <- list(
    objective = "reg:squarederror",
    eta = 0.1,
    max_depth = 4,
    subsample = 0.8,
    colsample_bytree = 0.8
  )
  
  # FIA XGBoost
  X_fia <- as.matrix(fia_model %>% select(ndvi_modis, tmean, ppt))
  y_fia <- fia_model$aglb_Mg_per_ha
  
  dtrain_fia <- xgb.DMatrix(data = X_fia, label = y_fia)
  xgb_fia <- xgb.train(data = dtrain_fia, nrounds = 50, params = xgb_params, verbose = 0)
  
  pred_fia_train <- predict(xgb_fia, dtrain_fia)
  fia_xgb_rmse <- sqrt(mean((y_fia - pred_fia_train)^2))
  fia_xgb_r2 <- 1 - sum((y_fia - pred_fia_train)^2) / sum((y_fia - mean(y_fia))^2)
  
  cat("  FIA XGBoost:\n")
  cat("    Train R² =", round(fia_xgb_r2, 4), "\n")
  cat("    Train RMSE =", round(fia_xgb_rmse, 2), "Mg/ha\n\n")
  
  # NEFIN XGBoost
  X_nefin <- as.matrix(nefin_train %>% select(ndvi_modis, tmean, ppt))
  y_nefin <- nefin_train$aglb_Mg_per_ha
  
  dtrain_nefin <- xgb.DMatrix(data = X_nefin, label = y_nefin)
  xgb_nefin <- xgb.train(data = dtrain_nefin, nrounds = 50, params = xgb_params, verbose = 0)
  
  pred_nefin_train <- predict(xgb_nefin, dtrain_nefin)
  nefin_xgb_rmse <- sqrt(mean((y_nefin - pred_nefin_train)^2))
  nefin_xgb_r2 <- 1 - sum((y_nefin - pred_nefin_train)^2) / sum((y_nefin - mean(y_nefin))^2)
  
  cat("  NEFIN XGBoost:\n")
  cat("    Train R² =", round(nefin_xgb_r2, 4), "\n")
  cat("    Train RMSE =", round(nefin_xgb_rmse, 2), "Mg/ha\n")
  
  # ===========================================================================
  # 4. EVALUATE ON HOLDOUT
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 4: Evaluating on NEFIN Holdout (True Coordinates & Biomass)\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  X_holdout <- as.matrix(nefin_holdout %>% select(ndvi_modis, tmean, ppt))
  y_holdout <- nefin_holdout$aglb_Mg_per_ha
  
  # Linear models on holdout
  pred_lm_fia_holdout <- predict(lm_fia, newdata = nefin_holdout)
  pred_lm_nefin_holdout <- predict(lm_nefin, newdata = nefin_holdout)
  
  rmse_lm_fia <- sqrt(mean((y_holdout - pred_lm_fia_holdout)^2))
  rmse_lm_nefin <- sqrt(mean((y_holdout - pred_lm_nefin_holdout)^2))
  
  r2_lm_fia <- 1 - sum((y_holdout - pred_lm_fia_holdout)^2) / sum((y_holdout - mean(y_holdout))^2)
  r2_lm_nefin <- 1 - sum((y_holdout - pred_lm_nefin_holdout)^2) / sum((y_holdout - mean(y_holdout))^2)
  
  cat("  Linear Models on Holdout:\n")
  cat("    FIA-trained:   RMSE =", round(rmse_lm_fia, 2), "| R² =", round(r2_lm_fia, 4), "\n")
  cat("    NEFIN-trained: RMSE =", round(rmse_lm_nefin, 2), "| R² =", round(r2_lm_nefin, 4), "\n\n")
  
  # XGBoost on holdout
  dholdout <- xgb.DMatrix(data = X_holdout)
  pred_xgb_fia_holdout <- predict(xgb_fia, dholdout)
  pred_xgb_nefin_holdout <- predict(xgb_nefin, dholdout)
  
  rmse_xgb_fia <- sqrt(mean((y_holdout - pred_xgb_fia_holdout)^2))
  rmse_xgb_nefin <- sqrt(mean((y_holdout - pred_xgb_nefin_holdout)^2))
  
  r2_xgb_fia <- 1 - sum((y_holdout - pred_xgb_fia_holdout)^2) / sum((y_holdout - mean(y_holdout))^2)
  r2_xgb_nefin <- 1 - sum((y_holdout - pred_xgb_nefin_holdout)^2) / sum((y_holdout - mean(y_holdout))^2)
  
  cat("  XGBoost Models on Holdout:\n")
  cat("    FIA-trained:   RMSE =", round(rmse_xgb_fia, 2), "| R² =", round(r2_xgb_fia, 4), "\n")
  cat("    NEFIN-trained: RMSE =", round(rmse_xgb_nefin, 2), "| R² =", round(r2_xgb_nefin, 4), "\n")
  
  # Delta RMSE
  delta_lm <- rmse_lm_fia - rmse_lm_nefin
  delta_xgb <- rmse_xgb_fia - rmse_xgb_nefin
  
  cat("\n─────────────────────────────────────────────\n")
  cat("HOLDOUT SUMMARY\n")
  cat("─────────────────────────────────────────────\n")
  cat("  Linear Model ΔRMSE (FIA - NEFIN):", round(delta_lm, 2), "Mg/ha",
      "(", round(100 * delta_lm / rmse_lm_fia, 1), "% improvement)\n")
  cat("  XGBoost ΔRMSE (FIA - NEFIN):     ", round(delta_xgb, 2), "Mg/ha",
      "(", round(100 * delta_xgb / rmse_xgb_fia, 1), "% improvement)\n")
  
  # ===========================================================================
  # 5. CREATE FIGURES
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 5: Creating Figures\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Holdout results dataframe
  holdout_results <- data.frame(
    model = c("Linear", "Linear", "XGBoost", "XGBoost"),
    trained_on = c("FIA (fuzzed)", "NEFIN (true)", "FIA (fuzzed)", "NEFIN (true)"),
    holdout_rmse = c(rmse_lm_fia, rmse_lm_nefin, rmse_xgb_fia, rmse_xgb_nefin),
    holdout_r2 = c(r2_lm_fia, r2_lm_nefin, r2_xgb_fia, r2_xgb_nefin)
  )
  
  # Figure 1: Holdout RMSE comparison
  fig1 <- ggplot(holdout_results, aes(x = model, y = holdout_rmse, fill = trained_on)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f", holdout_rmse)),
              position = position_dodge(0.9), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("FIA (fuzzed)" = "#e74c3c", "NEFIN (true)" = "#27ae60")) +
    labs(
      title = "Holdout Prediction RMSE: FIA vs NEFIN Training",
      subtitle = paste0("Tested on ", nrow(nefin_holdout), " NEFIN plots with true coordinates"),
      x = "Model Type", y = "RMSE (Mg/ha)", fill = "Trained On"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(fs::path(cfg$output_dir, "fig1_holdout_rmse.png"), fig1,
         width = 10, height = 6, dpi = 300)
  cat("  ✓ fig1_holdout_rmse.png\n")
  
  # Figure 2: Coefficient comparison
  coef_df <- data.frame(
    variable = c("NDVI", "Temperature", "Precipitation"),
    FIA = c(coef(lm_fia)["ndvi_modis"], coef(lm_fia)["tmean"], coef(lm_fia)["ppt"]),
    NEFIN = c(coef(lm_nefin)["ndvi_modis"], coef(lm_nefin)["tmean"], coef(lm_nefin)["ppt"])
  ) %>%
    pivot_longer(cols = c(FIA, NEFIN), names_to = "dataset", values_to = "coefficient")
  
  fig2 <- ggplot(coef_df, aes(x = variable, y = coefficient, fill = dataset)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("FIA" = "#e74c3c", "NEFIN" = "#27ae60")) +
    labs(
      title = "Linear Model Coefficients: FIA vs NEFIN",
      subtitle = "Attenuation in NDVI coefficient suggests measurement error from fuzzing",
      x = "", y = "Coefficient", fill = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(fs::path(cfg$output_dir, "fig2_coefficients.png"), fig2,
         width = 10, height = 6, dpi = 300)
  cat("  ✓ fig2_coefficients.png\n")
  
  # Figure 3: Observed vs Predicted (XGBoost)
  scatter_df <- data.frame(
    observed = y_holdout,
    pred_fia = pred_xgb_fia_holdout,
    pred_nefin = pred_xgb_nefin_holdout
  )
  
  fig3a <- ggplot(scatter_df, aes(x = observed, y = pred_fia)) +
    geom_point(alpha = 0.3, color = "#e74c3c") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    labs(title = "FIA-trained", 
         subtitle = paste0("RMSE = ", round(rmse_xgb_fia, 1), " | R² = ", round(r2_xgb_fia, 3)),
         x = "Observed (Mg/ha)", y = "Predicted (Mg/ha)") +
    coord_fixed(xlim = c(0, 400), ylim = c(0, 400)) +
    theme_minimal()
  
  fig3b <- ggplot(scatter_df, aes(x = observed, y = pred_nefin)) +
    geom_point(alpha = 0.3, color = "#27ae60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    labs(title = "NEFIN-trained",
         subtitle = paste0("RMSE = ", round(rmse_xgb_nefin, 1), " | R² = ", round(r2_xgb_nefin, 3)),
         x = "Observed (Mg/ha)", y = "Predicted (Mg/ha)") +
    coord_fixed(xlim = c(0, 400), ylim = c(0, 400)) +
    theme_minimal()
  
  fig3 <- fig3a + fig3b +
    plot_annotation(
      title = "XGBoost: Observed vs Predicted Biomass on Holdout",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
  
  ggsave(fs::path(cfg$output_dir, "fig3_observed_vs_predicted.png"), fig3,
         width = 12, height = 6, dpi = 300)
  cat("  ✓ fig3_observed_vs_predicted.png\n")
  
  # ===========================================================================
  # 6. SAVE RESULTS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 6: Saving Results\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  write_csv(holdout_results, fs::path(cfg$output_dir, "holdout_prediction_results.csv"))
  cat("  ✓ holdout_prediction_results.csv\n")
  
  # Summary
  summary_df <- data.frame(
    metric = c("FIA_n_plots", "NEFIN_train_n", "NEFIN_holdout_n",
               "LM_FIA_RMSE", "LM_NEFIN_RMSE", "LM_delta_RMSE",
               "XGB_FIA_RMSE", "XGB_NEFIN_RMSE", "XGB_delta_RMSE",
               "NDVI_coef_FIA", "NDVI_coef_NEFIN", "NDVI_attenuation_pct"),
    value = c(nrow(fia_model), nrow(nefin_train), nrow(nefin_holdout),
              rmse_lm_fia, rmse_lm_nefin, delta_lm,
              rmse_xgb_fia, rmse_xgb_nefin, delta_xgb,
              coef(lm_fia)["ndvi_modis"], coef(lm_nefin)["ndvi_modis"],
              100 * (coef(lm_nefin)["ndvi_modis"] - coef(lm_fia)["ndvi_modis"]) / coef(lm_nefin)["ndvi_modis"])
  )
  
  write_csv(summary_df, fs::path(cfg$output_dir, "model_comparison_summary.csv"))
  cat("  ✓ model_comparison_summary.csv\n")
  
  # ===========================================================================
  # FINAL SUMMARY
  # ===========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  ANALYSIS COMPLETE                                                        ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  cat("KEY RESULTS:\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  cat("  FIA plots (fuzzed):    ", nrow(fia_model), "\n")
  cat("  NEFIN plots (true):    ", nrow(nefin_train), "train +", nrow(nefin_holdout), "holdout\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  cat("  Linear Model ΔRMSE:    ", round(delta_lm, 2), "Mg/ha (", 
      round(100 * delta_lm / rmse_lm_fia, 1), "% improvement)\n")
  cat("  XGBoost ΔRMSE:         ", round(delta_xgb, 2), "Mg/ha (",
      round(100 * delta_xgb / rmse_xgb_fia, 1), "% improvement)\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  cat("  NDVI coefficient (FIA):   ", round(coef(lm_fia)["ndvi_modis"], 1), "\n")
  cat("  NDVI coefficient (NEFIN): ", round(coef(lm_nefin)["ndvi_modis"], 1), "\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  
  cat("\nOutput:", cfg$output_dir, "\n")
  
  invisible(list(
    holdout_results = holdout_results,
    lm_fia = lm_fia,
    lm_nefin = lm_nefin,
    xgb_fia = xgb_fia,
    xgb_nefin = xgb_nefin
  ))
}

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  run_model_comparison()
}
