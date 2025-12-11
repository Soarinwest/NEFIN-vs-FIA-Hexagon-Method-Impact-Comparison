#!/usr/bin/env Rscript
# ==============================================================================
# compare_spatial_models.R
# ==============================================================================
# Purpose:
#   Fit and compare spatial predictive models across different datasets:
#   - FIA fuzzed (single replicate)
#   - FIA mean (aggregated across replicates)
#   - NEFIN only (true coordinates)
#   - FIA + NEFIN combined
#
#   Evaluates whether NEFIN's spatial fidelity improves model accuracy for
#   predicting AGLB from NDVI and climate covariates.
#
# Inputs:
#   - data/processed/modeling/dataset_fia_fuzzed.csv
#   - data/processed/modeling/dataset_fia_mean.csv
#   - data/processed/modeling/dataset_nefin.csv
#   - data/processed/modeling/dataset_combined.csv
#
# Outputs:
#   - runs/model_comparison/cv_results.csv
#   - runs/model_comparison/model_summary.csv
#   - runs/model_comparison/figures/*.png
#   - runs/model_comparison/report.md
#
# Usage:
#   Rscript R/compare_spatial_models.R [--models=lm,gam,rf] [--cv-folds=5]
#                                       [--hex-scale=1kha]
#
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(fs)
  library(yaml)
  library(ggplot2)
  library(glue)
})

# Check optional packages
has_mgcv <- requireNamespace("mgcv", quietly = TRUE)
has_ranger <- requireNamespace("ranger", quietly = TRUE)

if (!has_mgcv) message("Note: mgcv not available - GAM models will be skipped")
if (!has_ranger) message("Note: ranger not available - RF models will be skipped")

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==============================================================================
# Configuration
# ==============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Defaults
  config <- list(
    models = c("lm", "gam", "rf"),
    cv_folds = 5,
    hex_scale = "1kha",
    output_dir = fs::path("runs", "model_comparison")
  )
  
  # Parse --models
  models_arg <- args[grepl("^--models=", args)]
  if (length(models_arg) > 0) {
    config$models <- strsplit(gsub("^--models=", "", models_arg[1]), ",")[[1]]
  }
  
  # Parse --cv-folds
  folds_arg <- args[grepl("^--cv-folds=", args)]
  if (length(folds_arg) > 0) {
    config$cv_folds <- as.integer(gsub("^--cv-folds=", "", folds_arg[1]))
  }
  
  # Parse --hex-scale
  scale_arg <- args[grepl("^--hex-scale=", args)]
  if (length(scale_arg) > 0) {
    config$hex_scale <- gsub("^--hex-scale=", "", scale_arg[1])
  }
  
  config
}

# ==============================================================================
# Data Loading
# ==============================================================================

load_datasets <- function() {
  
  data_dir <- fs::path("data", "processed", "modeling")
  
  datasets <- list()
  
  # FIA fuzzed
  path <- fs::path(data_dir, "dataset_fia_fuzzed.csv")
  if (fs::file_exists(path)) {
    datasets$fia_fuzzed <- read_csv(path, show_col_types = FALSE)
    message("Loaded FIA fuzzed: ", nrow(datasets$fia_fuzzed), " plots")
  }
  
  # FIA mean
  path <- fs::path(data_dir, "dataset_fia_mean.csv")
  if (fs::file_exists(path)) {
    datasets$fia_mean <- read_csv(path, show_col_types = FALSE)
    message("Loaded FIA mean: ", nrow(datasets$fia_mean), " plots")
  }
  
  # NEFIN
  path <- fs::path(data_dir, "dataset_nefin.csv")
  if (fs::file_exists(path)) {
    datasets$nefin <- read_csv(path, show_col_types = FALSE)
    message("Loaded NEFIN: ", nrow(datasets$nefin), " plots")
  }
  
  # Combined
  path <- fs::path(data_dir, "dataset_combined.csv")
  if (fs::file_exists(path)) {
    datasets$combined <- read_csv(path, show_col_types = FALSE)
    message("Loaded Combined: ", nrow(datasets$combined), " plots")
  }
  
  if (length(datasets) == 0) {
    stop("No datasets found in ", data_dir)
  }
  
  datasets
}

# ==============================================================================
# Data Preparation
# ==============================================================================

#' Prepare dataset for modeling
#' 
#' @param data Raw dataset
#' @param hex_scale Hex scale for spatial CV
#' @param response Response variable name
#' @return Clean dataset ready for modeling
prepare_for_modeling <- function(data, hex_scale = "1kha", 
                                  response = "aglb_Mg_per_ha") {
  
  hex_col <- paste0("hex_id_", hex_scale)
  
  # Required columns
  required <- c(response, "ndvi_modis")
  
  if (!all(required %in% names(data))) {
    missing <- setdiff(required, names(data))
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # Select modeling columns
  model_cols <- c("CN", "STATECD", "MEASYEAR", response,
                  "ndvi_modis", "ndvi_s2", "tmean", "ppt", "source")
  model_cols <- intersect(model_cols, names(data))
  
  # Add hex column if exists
  if (hex_col %in% names(data)) {
    model_cols <- c(model_cols, hex_col)
  }
  
  # Filter and clean
  clean <- data %>%
    select(all_of(model_cols)) %>%
    filter(!is.na(.data[[response]])) %>%
    filter(!is.na(ndvi_modis) | !is.na(ndvi_s2))
  
  # Rename hex column for consistency
  if (hex_col %in% names(clean)) {
    clean <- clean %>% rename(hex_id = !!sym(hex_col))
  }
  
  # Rename STATECD to state for consistency
  if ("STATECD" %in% names(clean)) {
    clean <- clean %>% rename(state = STATECD)
  }
  
  clean
}

# ==============================================================================
# Spatial Cross-Validation
# ==============================================================================

#' Create spatial CV folds based on hex IDs
#' 
#' @param data Dataset with hex_id column
#' @param k Number of folds
#' @return Dataset with fold column added
create_spatial_cv_folds <- function(data, k = 5) {
  
  if (!"hex_id" %in% names(data)) {
    # Fall back to random CV
    message("No hex_id found - using random CV")
    set.seed(42)
    return(data %>% mutate(fold = sample(rep(1:k, length.out = n()))))
  }
  
  # Get unique hexes
  hexes <- unique(data$hex_id)
  n_hex <- length(hexes)
  
  if (n_hex < k) {
    message("Only ", n_hex, " hexes - reducing to ", n_hex, " folds")
    k <- n_hex
  }
  
  # Randomly assign hexes to folds
  set.seed(42)
  hex_folds <- data.frame(
    hex_id = hexes,
    fold = sample(rep(1:k, length.out = n_hex)),
    stringsAsFactors = FALSE
  )
  
  data %>%
    left_join(hex_folds, by = "hex_id")
}

# ==============================================================================
# Model Fitting
# ==============================================================================

#' Build model formula
#' 
#' @param response Response variable
#' @param predictors Vector of predictor names
#' @param model_type "lm", "gam", or "rf"
#' @return Formula object
build_formula <- function(response, predictors, model_type = "lm") {
  
  # Filter to non-NA predictors
  if (model_type == "gam") {
    terms <- paste0("s(", predictors, ")")
    formula_str <- paste(response, "~", paste(terms, collapse = " + "))
  } else {
    formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
  }
  
  as.formula(formula_str)
}

#' Fit model and get predictions
#' 
#' @param train Training data
#' @param test Test data
#' @param predictors Predictor columns
#' @param response Response column
#' @param model_type Model type
#' @return Vector of predictions
fit_and_predict <- function(train, test, predictors, response, model_type) {
  
  # Check which predictors are available and non-NA
  available_preds <- predictors[
    sapply(predictors, function(p) {
      p %in% names(train) && sum(!is.na(train[[p]])) > 10
    })
  ]
  
  if (length(available_preds) == 0) {
    warning("No valid predictors available")
    return(rep(NA, nrow(test)))
  }
  
  # Build formula
  formula <- build_formula(response, available_preds, model_type)
  
  # Remove rows with NA in predictors
  train_clean <- train %>%
    filter(if_all(all_of(available_preds), ~ !is.na(.)))
  
  if (nrow(train_clean) < 20) {
    warning("Not enough training data after removing NAs")
    return(rep(NA, nrow(test)))
  }
  
  # Fit model
  model <- tryCatch({
    if (model_type == "lm") {
      lm(formula, data = train_clean)
    } else if (model_type == "gam" && has_mgcv) {
      mgcv::gam(formula, data = train_clean, method = "REML")
    } else if (model_type == "rf" && has_ranger) {
      ranger::ranger(formula, data = train_clean, num.trees = 500, seed = 42)
    } else {
      NULL
    }
  }, error = function(e) {
    warning("Model fitting failed: ", e$message)
    NULL
  })
  
  if (is.null(model)) {
    return(rep(NA, nrow(test)))
  }
  
  # Predict
  preds <- tryCatch({
    if (model_type == "rf") {
      predict(model, data = test)$predictions
    } else {
      predict(model, newdata = test)
    }
  }, error = function(e) {
    rep(NA, nrow(test))
  })
  
  preds
}

#' Calculate performance metrics
#' 
#' @param observed Observed values
#' @param predicted Predicted values
#' @return Named list of metrics
calc_metrics <- function(observed, predicted) {
  
  valid <- !is.na(observed) & !is.na(predicted)
  obs <- observed[valid]
  pred <- predicted[valid]
  
  if (length(obs) < 2) {
    return(list(rmse = NA, mae = NA, r2 = NA, bias = NA, n = 0))
  }
  
  residuals <- obs - pred
  
  list(
    rmse = sqrt(mean(residuals^2)),
    mae = mean(abs(residuals)),
    r2 = cor(obs, pred)^2,
    bias = mean(residuals),
    n = length(obs)
  )
}

# ==============================================================================
# Run CV for Single Dataset
# ==============================================================================

#' Run cross-validation for a single dataset
#' 
#' @param data Prepared dataset with fold column
#' @param dataset_name Name for output
#' @param model_types Vector of model types
#' @param response Response variable
#' @return data.frame of CV results
run_cv <- function(data, dataset_name, model_types, response = "aglb_Mg_per_ha") {
  
  folds <- sort(unique(data$fold))
  predictors <- c("ndvi_modis", "ndvi_s2", "tmean", "ppt")
  predictors <- intersect(predictors, names(data))
  
  message("  Running CV for ", dataset_name, " (", length(folds), " folds, ",
          length(predictors), " predictors)")
  
  results <- list()
  
  for (fold_i in folds) {
    train <- data %>% filter(fold != fold_i)
    test <- data %>% filter(fold == fold_i)
    
    for (model_type in model_types) {
      
      # Skip unavailable models
      if (model_type == "gam" && !has_mgcv) next
      if (model_type == "rf" && !has_ranger) next
      
      # Fit and predict
      preds <- fit_and_predict(train, test, predictors, response, model_type)
      
      # Calculate metrics
      metrics <- calc_metrics(test[[response]], preds)
      
      results[[length(results) + 1]] <- tibble(
        dataset = dataset_name,
        model = model_type,
        fold = fold_i,
        rmse = metrics$rmse,
        mae = metrics$mae,
        r2 = metrics$r2,
        bias = metrics$bias,
        n_test = metrics$n,
        n_train = nrow(train)
      )
    }
  }
  
  bind_rows(results)
}

# ==============================================================================
# Generate Report
# ==============================================================================

generate_report <- function(summary, delta, output_dir) {
  
  lines <- c(
    "# Spatial Model Comparison Report",
    "",
    glue("**Generated:** {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}"),
    "",
    "## Overview",
    "",
    "This report compares predictive models for Above-Ground Live Biomass (AGLB)",
    "using NDVI and climate covariates across datasets with different spatial fidelity:",
    "",
    "- **FIA fuzzed**: Single jitter replicate (±1.6 km uncertainty)",
    "- **FIA mean**: Mean across 80 jitter replicates",
    "- **NEFIN**: True plot coordinates",
    "- **Combined**: FIA (mean) + NEFIN together",
    "",
    "## Model Performance",
    "",
    "| Dataset | Model | RMSE (Mg/ha) | R² | MAE |",
    "|---------|-------|--------------|-----|-----|"
  )
  
  for (i in 1:nrow(summary)) {
    r <- summary[i, ]
    lines <- c(lines,
      glue("| {r$dataset} | {r$model} | {round(r$rmse_mean, 2)} ± {round(r$rmse_sd, 2)} | {round(r$r2_mean * 100, 1)}% | {round(r$mae_mean, 2)} |")
    )
  }
  
  lines <- c(lines,
    "",
    "## Key Findings",
    ""
  )
  
  # Find best dataset
  best <- summary %>%
    filter(model == "lm") %>%
    slice_min(rmse_mean)
  
  lines <- c(lines,
    glue("**Best performing dataset (LM):** {best$dataset}"),
    glue("- RMSE: {round(best$rmse_mean, 2)} Mg/ha"),
    glue("- R²: {round(best$r2_mean * 100, 1)}%"),
    ""
  )
  
  # NEFIN impact
  if (!is.null(delta) && "delta_rmse" %in% names(delta)) {
    combined_delta <- delta %>% filter(comparison == "combined_vs_fia_mean")
    
    if (nrow(combined_delta) > 0) {
      d <- combined_delta[1, ]
      if (!is.na(d$delta_rmse) && d$delta_rmse < 0) {
        lines <- c(lines,
          "**NEFIN Impact:** Adding NEFIN to FIA **improves** predictions",
          glue("- ΔRMSE: {round(d$delta_rmse, 2)} Mg/ha (negative = improvement)"),
          glue("- ΔR²: {round(d$delta_r2 * 100, 2)} percentage points"),
          ""
        )
      } else {
        lines <- c(lines,
          "**NEFIN Impact:** Adding NEFIN did not improve predictions",
          glue("- ΔRMSE: {round(d$delta_rmse, 2)} Mg/ha"),
          ""
        )
      }
    }
  }
  
  lines <- c(lines,
    "## Recommendations",
    "",
    "Based on this analysis:",
    ""
  )
  
  # Add recommendations based on results
  if (!is.null(delta) && any(delta$delta_rmse < 0, na.rm = TRUE)) {
    lines <- c(lines,
      "1. **Include NEFIN data** in spatial modeling for improved accuracy",
      "2. The true coordinates provide value beyond FIA's fuzzed locations",
      ""
    )
  } else {
    lines <- c(lines,
      "1. **FIA alone may be sufficient** for NDVI-based AGLB predictions",
      "2. Coordinate fuzzing has limited impact at this spatial scale",
      ""
    )
  }
  
  lines <- c(lines,
    "---",
    "",
    glue("*Report generated by compare_spatial_models.R*")
  )
  
  report_path <- fs::path(output_dir, "report.md")
  writeLines(lines, report_path)
  message("Saved report: ", report_path)
}

# ==============================================================================
# Main Entry Point
# ==============================================================================

main <- function() {
  
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  SPATIAL MODEL COMPARISON\n")
  cat(strrep("=", 70), "\n\n")
  
  # Parse configuration
  config <- parse_args()
  
  message("Configuration:")
  message("  Models: ", paste(config$models, collapse = ", "))
  message("  CV folds: ", config$cv_folds)
  message("  Hex scale: ", config$hex_scale)
  
  # Create output directory
  fs::dir_create(config$output_dir, recurse = TRUE)
  figures_dir <- fs::path(config$output_dir, "figures")
  fs::dir_create(figures_dir)
  
  # Load datasets
  message("\nLoading datasets...")
  datasets <- load_datasets()
  
  # Prepare datasets
  message("\nPreparing datasets for modeling...")
  prepared <- list()
  
  for (name in names(datasets)) {
    prepared[[name]] <- datasets[[name]] %>%
      prepare_for_modeling(hex_scale = config$hex_scale) %>%
      create_spatial_cv_folds(k = config$cv_folds)
    
    message("  ", name, ": ", nrow(prepared[[name]]), " plots, ",
            length(unique(prepared[[name]]$fold)), " folds")
  }
  
  # Run CV for each dataset
  message("\nRunning cross-validation...")
  all_results <- list()
  
  for (name in names(prepared)) {
    results <- run_cv(prepared[[name]], name, config$models)
    all_results[[name]] <- results
  }
  
  cv_results <- bind_rows(all_results)
  
  # Save CV results
  cv_path <- fs::path(config$output_dir, "cv_results.csv")
  write_csv(cv_results, cv_path)
  message("\nSaved CV results: ", cv_path)
  
  # Summarize results
  summary <- cv_results %>%
    group_by(dataset, model) %>%
    summarise(
      rmse_mean = mean(rmse, na.rm = TRUE),
      rmse_sd = sd(rmse, na.rm = TRUE),
      mae_mean = mean(mae, na.rm = TRUE),
      r2_mean = mean(r2, na.rm = TRUE),
      r2_sd = sd(r2, na.rm = TRUE),
      n_folds = n(),
      .groups = "drop"
    )
  
  summary_path <- fs::path(config$output_dir, "model_summary.csv")
  write_csv(summary, summary_path)
  message("Saved summary: ", summary_path)
  
  # Calculate delta metrics
  delta <- NULL
  if ("combined" %in% names(prepared) && "fia_mean" %in% names(prepared)) {
    delta <- summary %>%
      select(dataset, model, rmse_mean, r2_mean) %>%
      pivot_wider(names_from = dataset, values_from = c(rmse_mean, r2_mean)) %>%
      mutate(
        comparison = "combined_vs_fia_mean",
        delta_rmse = rmse_mean_combined - rmse_mean_fia_mean,
        delta_r2 = r2_mean_combined - r2_mean_fia_mean
      )
    
    delta_path <- fs::path(config$output_dir, "delta_metrics.csv")
    write_csv(delta, delta_path)
    message("Saved delta metrics: ", delta_path)
  }
  
  # Generate plots
  message("\nGenerating plots...")
  
  # RMSE comparison
  p1 <- ggplot(summary, aes(x = model, y = rmse_mean, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(
      aes(ymin = rmse_mean - rmse_sd, ymax = rmse_mean + rmse_sd),
      position = position_dodge(width = 0.8),
      width = 0.2
    ) +
    scale_fill_brewer(palette = "Set2", name = "Dataset") +
    labs(
      title = "Model RMSE by Dataset",
      x = "Model Type",
      y = "RMSE (Mg/ha)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(fs::path(figures_dir, "rmse_comparison.png"), p1, width = 10, height = 6, dpi = 150)
  
  # R² comparison
  p2 <- ggplot(summary, aes(x = model, y = r2_mean, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_brewer(palette = "Set2", name = "Dataset") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
      title = "Model R² by Dataset",
      x = "Model Type",
      y = "R²"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(fs::path(figures_dir, "r2_comparison.png"), p2, width = 10, height = 6, dpi = 150)
  
  message("Saved plots to: ", figures_dir)
  
  # Generate report
  generate_report(summary, delta, config$output_dir)
  
  # Print summary
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  RESULTS SUMMARY\n")
  cat(strrep("=", 70), "\n\n")
  
  print(summary %>% arrange(model, rmse_mean))
  
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  COMPARISON COMPLETE\n")
  cat("  Outputs: ", config$output_dir, "\n")
  cat(strrep("=", 70), "\n\n")
}

# Run if called directly
if (!interactive()) {
  main()
}
