#!/usr/bin/env Rscript
# ==============================================================================
# ndvi_model_comparison.R
# ==============================================================================
# Purpose:
#   Fit and compare spatial predictive models for AGLB using NDVI and climate
#   covariates. Compares FIA-only (fuzzed coordinates) vs FIA+NEFIN (combined
#   with true coordinates) to quantify the value of high-fidelity plot locations
#   for spatial prediction.
#
# Inputs:
#   - data/processed/ndvi/modeling_dataset_fia.csv (FIA plots with NDVI)
#   - data/processed/ndvi/modeling_dataset_nefin.csv (NEFIN plots with NDVI)
#   - data/processed/ndvi/modeling_dataset_combined.csv (both datasets)
#   - configs/process.yml (modeling configuration)
#
# Outputs:
#   - runs/consolidated_.../ndvi_model_comparison.csv (performance metrics)
#   - runs/consolidated_.../ndvi_model_figures/*.png (comparison plots)
#   - runs/consolidated_.../ndvi_model_summary.yml (summary statistics)
#
# Usage:
#   Rscript R/ndvi_model_comparison.R [--overwrite] [--models=lm,gam,rf]
#                                      [--cv-folds=5] [--output-dir=path]
#
# Pipeline integration:
#   Run after build_ndvi_plot_dataset.R, before ndvi_model_report.R
#
# Author: Claude (Anthropic)
# Date: 2025-12
# ==============================================================================

# --- Setup -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(yaml)
  library(fs)
  library(sf)
  library(glue)
  library(ggplot2)
  library(cli)
})

# Check for optional packages
has_mgcv <- requireNamespace("mgcv", quietly = TRUE)
has_ranger <- requireNamespace("ranger", quietly = TRUE)

if (!has_mgcv) {
  cli_alert_warning("Package 'mgcv' not available - GAM models will be skipped")
}
if (!has_ranger) {
  cli_alert_warning("Package 'ranger' not available - Random Forest models will be skipped")
}

# --- Configuration -----------------------------------------------------------

# Load configuration
config_path <- fs::path("configs", "process.yml")
if (!fs::file_exists(config_path)) {
  cli_abort("Configuration file not found: {config_path}")
}
config <- yaml::read_yaml(config_path)

# Extract modeling config with defaults
model_config <- config$ndvi_modeling %||% list()
cv_method <- model_config$cv_method %||% "spatial_hex"
cv_folds <- model_config$cv_folds %||% 5
models_to_fit <- model_config$models %||% c("lm", "gam", "rf")
response_var <- model_config$response %||% "aglb"
include_climate <- model_config$include_climate %||% TRUE
min_plots_per_hex <- model_config$min_plots_per_hex %||% 3
fuzz_pressure_quantiles <- model_config$fuzz_pressure_quantiles %||% c(0.25, 0.5, 0.75)

# Paths
ndvi_dir <- fs::path("data", "processed", "ndvi")
fia_data_path <- fs::path(ndvi_dir, "modeling_dataset_fia.csv")
nefin_data_path <- fs::path(ndvi_dir, "modeling_dataset_nefin.csv")
combined_data_path <- fs::path(ndvi_dir, "modeling_dataset_combined.csv")

# --- Parse CLI Arguments -----------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

overwrite <- "--overwrite" %in% args

# Parse --models flag
models_arg <- args[grepl("^--models=", args)]
if (length(models_arg) > 0) {
  models_to_fit <- strsplit(gsub("^--models=", "", models_arg[1]), ",")[[1]]
}

# Parse --cv-folds flag
cv_folds_arg <- args[grepl("^--cv-folds=", args)]
if (length(cv_folds_arg) > 0) {
  cv_folds <- as.integer(gsub("^--cv-folds=", "", cv_folds_arg[1]))
}

# Parse --output-dir flag
output_dir_arg <- args[grepl("^--output-dir=", args)]
if (length(output_dir_arg) > 0) {
  output_dir <- fs::path(gsub("^--output-dir=", "", output_dir_arg[1]))
} else {
  # Find most recent consolidated run directory
  runs_dir <- fs::path("runs")
  if (fs::dir_exists(runs_dir)) {
    run_dirs <- fs::dir_ls(runs_dir, type = "directory", regexp = "consolidated_")
    if (length(run_dirs) > 0) {
      output_dir <- run_dirs[length(run_dirs)]  # Most recent
    } else {
      output_dir <- fs::path(runs_dir, glue("consolidated_{format(Sys.time(), '%Y%m%d_%H%M%S')}"))
    }
  } else {
    output_dir <- fs::path(runs_dir, glue("consolidated_{format(Sys.time(), '%Y%m%d_%H%M%S')}"))
  }
}

# Create output directories
fs::dir_create(output_dir, recurse = TRUE)
figures_dir <- fs::path(output_dir, "ndvi_model_figures")
fs::dir_create(figures_dir)

cli_h1("NDVI Spatial Model Comparison")
cli_alert_info("Output directory: {output_dir}")

# --- Load Data ---------------------------------------------------------------

cli_h2("Loading datasets")

if (!fs::file_exists(fia_data_path)) {
  cli_abort("FIA modeling dataset not found: {fia_data_path}")
}
if (!fs::file_exists(combined_data_path)) {
  cli_abort("Combined modeling dataset not found: {combined_data_path}")
}

fia_data <- read_csv(fia_data_path, show_col_types = FALSE)
combined_data <- read_csv(combined_data_path, show_col_types = FALSE)

# Also load NEFIN-only if available
nefin_data <- NULL
if (fs::file_exists(nefin_data_path)) {
  nefin_data <- read_csv(nefin_data_path, show_col_types = FALSE)
}

cli_alert_success("FIA dataset: {nrow(fia_data)} plots")
cli_alert_success("Combined dataset: {nrow(combined_data)} plots")
if (!is.null(nefin_data)) {
  cli_alert_success("NEFIN dataset: {nrow(nefin_data)} plots")
}

# --- Data Preparation --------------------------------------------------------

cli_h2("Preparing data for modeling")

#' Prepare dataset for modeling
#' 
#' @param data Data frame with plot-level data
#' @param response Response variable name
#' @param include_climate Whether to include climate covariates
#' @param ndvi_vars Which NDVI variables to include
#' @return Cleaned data frame ready for modeling
prepare_model_data <- function(data, response = "aglb", include_climate = TRUE,
                                ndvi_vars = c("ndvi_modis", "ndvi_s2")) {
  
  # Start with response and identifiers
  keep_cols <- c("plot_id", "hex_id", "state", response)
  
  # Add NDVI variables (use available ones)
  available_ndvi <- intersect(ndvi_vars, names(data))
  if (length(available_ndvi) == 0) {
    cli_abort("No NDVI variables found in dataset")
  }
  keep_cols <- c(keep_cols, available_ndvi)
  
  # Add uncertainty columns if available
  ndvi_sd_cols <- paste0(available_ndvi, "_sd")
  available_sd <- intersect(ndvi_sd_cols, names(data))
  keep_cols <- c(keep_cols, available_sd)
  
  # Add climate if requested
  if (include_climate) {
    climate_cols <- c("tmean_mean", "ppt_mean")
    available_climate <- intersect(climate_cols, names(data))
    keep_cols <- c(keep_cols, available_climate)
  }
  
  # Add source indicator if present
  if ("source" %in% names(data)) {
    keep_cols <- c(keep_cols, "source")
  }
  
  # Select and filter
  data_clean <- data %>%
    select(any_of(keep_cols)) %>%
    filter(!is.na(.data[[response]])) %>%
    filter(if_any(all_of(available_ndvi), ~ !is.na(.)))
  
  return(data_clean)
}

# Prepare both datasets
fia_model_data <- prepare_model_data(fia_data, response_var, include_climate)
combined_model_data <- prepare_model_data(combined_data, response_var, include_climate)

cli_alert_info("FIA model-ready plots: {nrow(fia_model_data)}")
cli_alert_info("Combined model-ready plots: {nrow(combined_model_data)}")

# --- Spatial Cross-Validation Setup ------------------------------------------

cli_h2("Setting up spatial cross-validation")

#' Create spatial CV folds based on hexagons
#' 
#' @param data Data frame with hex_id column
#' @param k Number of folds
#' @param min_plots_per_fold Minimum plots per fold
#' @return Data frame with fold assignments
create_hex_cv_folds <- function(data, k = 5, min_plots_per_fold = 10) {
  
  # Get unique hexes
  hex_ids <- unique(data$hex_id)
  n_hex <- length(hex_ids)
  
  if (n_hex < k) {
    cli_alert_warning("Only {n_hex} hexes available, reducing folds to {n_hex}")
    k <- n_hex
  }
  
  # Randomly assign hexes to folds
  set.seed(42)  # Reproducibility
  hex_folds <- data.frame(
    hex_id = hex_ids,
    fold = sample(rep(1:k, length.out = n_hex))
  )
  
  # Join back to data
  data_with_folds <- data %>%
    left_join(hex_folds, by = "hex_id")
  
  # Check fold sizes
  fold_sizes <- data_with_folds %>%
    count(fold, name = "n_plots")
  
  cli_alert_info("Fold sizes: {paste(fold_sizes$n_plots, collapse = ', ')}")
  
  return(data_with_folds)
}

#' Create spatial CV folds based on state
#' 
#' @param data Data frame with state column
#' @param k Number of folds (may be adjusted)
#' @return Data frame with fold assignments
create_state_cv_folds <- function(data, k = 5) {
  
  # Get unique states
  states <- unique(data$state)
  n_states <- length(states)
  
  if (n_states < k) {
    cli_alert_warning("Only {n_states} states available, reducing folds to {n_states}")
    k <- n_states
  }
  
  # Assign states to folds (roughly equal plot counts)
  state_counts <- data %>%
    count(state, name = "n_plots") %>%
    arrange(desc(n_plots))
  
  # Greedy assignment to balance folds
  fold_totals <- rep(0, k)
  state_folds <- data.frame(state = character(), fold = integer())
  
  for (i in 1:nrow(state_counts)) {
    # Assign to smallest fold
    target_fold <- which.min(fold_totals)
    state_folds <- bind_rows(state_folds,
                              data.frame(state = state_counts$state[i],
                                        fold = target_fold))
    fold_totals[target_fold] <- fold_totals[target_fold] + state_counts$n_plots[i]
  }
  
  # Join back
  data_with_folds <- data %>%
    left_join(state_folds, by = "state")
  
  return(data_with_folds)
}

#' Create random k-fold CV (non-spatial)
#' 
#' @param data Data frame
#' @param k Number of folds
#' @return Data frame with fold assignments
create_random_cv_folds <- function(data, k = 5) {
  set.seed(42)
  data %>%
    mutate(fold = sample(rep(1:k, length.out = n())))
}

# Create folds based on method
if (cv_method == "spatial_hex") {
  fia_cv <- create_hex_cv_folds(fia_model_data, k = cv_folds)
  combined_cv <- create_hex_cv_folds(combined_model_data, k = cv_folds)
} else if (cv_method == "spatial_state") {
  fia_cv <- create_state_cv_folds(fia_model_data, k = cv_folds)
  combined_cv <- create_state_cv_folds(combined_model_data, k = cv_folds)
} else {
  fia_cv <- create_random_cv_folds(fia_model_data, k = cv_folds)
  combined_cv <- create_random_cv_folds(combined_model_data, k = cv_folds)
}

cli_alert_success("CV method: {cv_method} with {cv_folds} folds")

# --- Model Fitting Functions -------------------------------------------------

cli_h2("Model fitting functions")

#' Build model formula
#' 
#' @param response Response variable name
#' @param ndvi_vars NDVI predictors
#' @param climate_vars Climate predictors
#' @param model_type Type of model (lm, gam, rf)
#' @return Formula object
build_formula <- function(response, ndvi_vars, climate_vars = NULL, model_type = "lm") {
  
  predictors <- ndvi_vars
  if (!is.null(climate_vars) && length(climate_vars) > 0) {
    predictors <- c(predictors, climate_vars)
  }
  
  if (model_type == "gam") {
    # Use smooth terms for GAM
    smooth_terms <- paste0("s(", predictors, ")")
    formula_str <- paste(response, "~", paste(smooth_terms, collapse = " + "))
  } else {
    # Linear terms for lm and rf
    formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
  }
  
  as.formula(formula_str)
}

#' Fit linear model
#' 
#' @param train_data Training data
#' @param formula Model formula
#' @return Fitted model object
fit_lm <- function(train_data, formula) {
  lm(formula, data = train_data)
}

#' Fit GAM model
#' 
#' @param train_data Training data
#' @param formula Model formula
#' @return Fitted model object
fit_gam <- function(train_data, formula) {
  if (!has_mgcv) {
    return(NULL)
  }
  mgcv::gam(formula, data = train_data, method = "REML")
}

#' Fit Random Forest model
#' 
#' @param train_data Training data
#' @param formula Model formula
#' @return Fitted model object
fit_rf <- function(train_data, formula) {
  if (!has_ranger) {
    return(NULL)
  }
  ranger::ranger(formula, data = train_data, num.trees = 500, seed = 42)
}

#' Predict from model
#' 
#' @param model Fitted model
#' @param newdata Data to predict on
#' @param model_type Type of model
#' @return Vector of predictions
predict_model <- function(model, newdata, model_type) {
  if (is.null(model)) return(rep(NA, nrow(newdata)))
  
  if (model_type == "rf") {
    predict(model, data = newdata)$predictions
  } else {
    predict(model, newdata = newdata)
  }
}

#' Calculate performance metrics
#' 
#' @param observed Observed values
#' @param predicted Predicted values
#' @return Named list of metrics
calc_metrics <- function(observed, predicted) {
  
  # Remove NAs
  valid <- !is.na(observed) & !is.na(predicted)
  obs <- observed[valid]
  pred <- predicted[valid]
  
  if (length(obs) < 2) {
    return(list(rmse = NA, mae = NA, r2 = NA, bias = NA, cor = NA, n = 0))
  }
  
  residuals <- obs - pred
  
  list(
    rmse = sqrt(mean(residuals^2)),
    mae = mean(abs(residuals)),
    r2 = cor(obs, pred)^2,
    bias = mean(residuals),
    cor = cor(obs, pred),
    n = length(obs)
  )
}

# --- Cross-Validation Loop ---------------------------------------------------

cli_h2("Running cross-validation")

#' Run cross-validation for a single dataset
#' 
#' @param data Data with fold assignments
#' @param response Response variable
#' @param model_types Vector of model types to fit
#' @param include_climate Include climate predictors
#' @param dataset_name Name for output
#' @return Data frame of CV results
run_cv <- function(data, response, model_types, include_climate, dataset_name) {
  
  folds <- sort(unique(data$fold))
  
  # Identify available predictors
  ndvi_vars <- intersect(c("ndvi_modis", "ndvi_s2"), names(data))
  climate_vars <- if (include_climate) {
    intersect(c("tmean_mean", "ppt_mean"), names(data))
  } else {
    NULL
  }
  
  # Storage for results
  results <- list()
  
  # Progress
  pb <- cli_progress_bar(
    format = "{dataset_name} | Fold {fold}/{n_folds} | {model} | {cli::pb_percent}",
    total = length(folds) * length(model_types)
  )
  
  for (fold_i in folds) {
    # Split data
    train_data <- data %>% filter(fold != fold_i)
    test_data <- data %>% filter(fold == fold_i)
    
    for (model_type in model_types) {
      
      cli_progress_update(pb, set = which(folds == fold_i) * length(model_types) - 
                            (length(model_types) - which(model_types == model_type)))
      
      # Skip unavailable models
      if (model_type == "gam" && !has_mgcv) next
      if (model_type == "rf" && !has_ranger) next
      
      # Build formula
      formula <- build_formula(response, ndvi_vars, climate_vars, model_type)
      
      # Fit model
      model <- tryCatch({
        if (model_type == "lm") fit_lm(train_data, formula)
        else if (model_type == "gam") fit_gam(train_data, formula)
        else if (model_type == "rf") fit_rf(train_data, formula)
      }, error = function(e) {
        cli_alert_warning("Model {model_type} failed in fold {fold_i}: {e$message}")
        NULL
      })
      
      if (is.null(model)) next
      
      # Predict
      predictions <- predict_model(model, test_data, model_type)
      
      # Calculate metrics
      metrics <- calc_metrics(test_data[[response]], predictions)
      
      # Store
      results[[length(results) + 1]] <- tibble(
        dataset = dataset_name,
        model = model_type,
        fold = fold_i,
        rmse = metrics$rmse,
        mae = metrics$mae,
        r2 = metrics$r2,
        bias = metrics$bias,
        cor = metrics$cor,
        n_test = metrics$n,
        n_train = nrow(train_data)
      )
    }
  }
  
  cli_progress_done(pb)
  
  bind_rows(results)
}

# Run CV for FIA-only
cli_alert("Running CV for FIA-only dataset...")
fia_results <- run_cv(fia_cv, response_var, models_to_fit, include_climate, "FIA-only")

# Run CV for Combined
cli_alert("Running CV for Combined dataset...")
combined_results <- run_cv(combined_cv, response_var, models_to_fit, include_climate, "FIA+NEFIN")

# Combine results
all_results <- bind_rows(fia_results, combined_results)

cli_alert_success("CV complete: {nrow(all_results)} model-fold combinations")

# --- Calculate Summary Statistics --------------------------------------------

cli_h2("Calculating summary statistics")

#' Summarize CV results by model and dataset
#' 
#' @param results CV results data frame
#' @return Summary data frame
summarize_cv_results <- function(results) {
  results %>%
    group_by(dataset, model) %>%
    summarise(
      rmse_mean = mean(rmse, na.rm = TRUE),
      rmse_sd = sd(rmse, na.rm = TRUE),
      mae_mean = mean(mae, na.rm = TRUE),
      mae_sd = sd(mae, na.rm = TRUE),
      r2_mean = mean(r2, na.rm = TRUE),
      r2_sd = sd(r2, na.rm = TRUE),
      bias_mean = mean(bias, na.rm = TRUE),
      bias_sd = sd(bias, na.rm = TRUE),
      n_folds = n(),
      .groups = "drop"
    )
}

cv_summary <- summarize_cv_results(all_results)

# Calculate delta metrics (Combined - FIA)
delta_metrics <- cv_summary %>%
  select(dataset, model, rmse_mean, r2_mean, mae_mean) %>%
  pivot_wider(
    names_from = dataset,
    values_from = c(rmse_mean, r2_mean, mae_mean)
  ) %>%
  mutate(
    delta_rmse = `rmse_mean_FIA+NEFIN` - `rmse_mean_FIA-only`,
    delta_r2 = `r2_mean_FIA+NEFIN` - `r2_mean_FIA-only`,
    delta_mae = `mae_mean_FIA+NEFIN` - `mae_mean_FIA-only`,
    rmse_pct_change = (delta_rmse / `rmse_mean_FIA-only`) * 100,
    r2_pct_change = (delta_r2 / `r2_mean_FIA-only`) * 100
  )

cli_alert_info("Delta metrics calculated")

# Print summary
cli_h3("Performance Summary")
print(cv_summary)

cli_h3("Delta Metrics (Combined vs FIA-only)")
print(delta_metrics %>% select(model, delta_rmse, delta_r2, rmse_pct_change, r2_pct_change))

# --- Visualization -----------------------------------------------------------

cli_h2("Generating visualizations")

# Set theme
theme_set(theme_minimal(base_size = 12))

#' Plot RMSE comparison
#' 
#' @param summary CV summary data frame
#' @param output_path Path to save plot
plot_rmse_comparison <- function(summary, output_path) {
  
  p <- ggplot(summary, aes(x = model, y = rmse_mean, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(
      aes(ymin = rmse_mean - rmse_sd, ymax = rmse_mean + rmse_sd),
      position = position_dodge(width = 0.8),
      width = 0.2
    ) +
    scale_fill_manual(
      values = c("FIA-only" = "#E69F00", "FIA+NEFIN" = "#56B4E9"),
      name = "Dataset"
    ) +
    labs(
      title = "Model RMSE Comparison: FIA-only vs FIA+NEFIN",
      subtitle = glue("CV method: {cv_method}, {cv_folds} folds"),
      x = "Model Type",
      y = "RMSE (Mg/ha)",
      caption = "Error bars show ± 1 SD across CV folds"
    ) +
    theme(legend.position = "bottom")
  
  ggsave(output_path, p, width = 8, height = 6, dpi = 150)
  cli_alert_success("Saved: {output_path}")
}

#' Plot R² comparison
#' 
#' @param summary CV summary data frame
#' @param output_path Path to save plot
plot_r2_comparison <- function(summary, output_path) {
  
  p <- ggplot(summary, aes(x = model, y = r2_mean, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(
      aes(ymin = pmax(0, r2_mean - r2_sd), ymax = pmin(1, r2_mean + r2_sd)),
      position = position_dodge(width = 0.8),
      width = 0.2
    ) +
    scale_fill_manual(
      values = c("FIA-only" = "#E69F00", "FIA+NEFIN" = "#56B4E9"),
      name = "Dataset"
    ) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title = "Model R² Comparison: FIA-only vs FIA+NEFIN",
      subtitle = glue("CV method: {cv_method}, {cv_folds} folds"),
      x = "Model Type",
      y = "R² (explained variance)",
      caption = "Error bars show ± 1 SD across CV folds"
    ) +
    theme(legend.position = "bottom")
  
  ggsave(output_path, p, width = 8, height = 6, dpi = 150)
  cli_alert_success("Saved: {output_path}")
}

#' Plot delta metrics
#' 
#' @param delta Delta metrics data frame
#' @param output_path Path to save plot
plot_delta_metrics <- function(delta, output_path) {
  
  # Reshape for plotting
  delta_long <- delta %>%
    select(model, delta_rmse, delta_r2) %>%
    pivot_longer(
      cols = c(delta_rmse, delta_r2),
      names_to = "metric",
      values_to = "delta"
    ) %>%
    mutate(
      metric_label = case_when(
        metric == "delta_rmse" ~ "ΔRMSE (Mg/ha)",
        metric == "delta_r2" ~ "ΔR²"
      ),
      improvement = if_else(
        (metric == "delta_rmse" & delta < 0) | (metric == "delta_r2" & delta > 0),
        "Improved", "Degraded"
      )
    )
  
  p <- ggplot(delta_long, aes(x = model, y = delta, fill = improvement)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    facet_wrap(~metric_label, scales = "free_y") +
    scale_fill_manual(
      values = c("Improved" = "#2E7D32", "Degraded" = "#C62828"),
      name = ""
    ) +
    labs(
      title = "Change in Model Performance: FIA+NEFIN vs FIA-only",
      subtitle = "Negative ΔRMSE and positive ΔR² indicate improvement",
      x = "Model Type",
      y = "Change",
      caption = "Delta = Combined - FIA-only"
    ) +
    theme(legend.position = "bottom")
  
  ggsave(output_path, p, width = 10, height = 5, dpi = 150)
  cli_alert_success("Saved: {output_path}")
}

#' Plot fold-level variability
#' 
#' @param results Full CV results
#' @param output_path Path to save plot
plot_fold_variability <- function(results, output_path) {
  
  p <- ggplot(results, aes(x = factor(fold), y = rmse, color = dataset)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(aes(group = dataset), alpha = 0.5) +
    facet_wrap(~model, scales = "free_y") +
    scale_color_manual(
      values = c("FIA-only" = "#E69F00", "FIA+NEFIN" = "#56B4E9"),
      name = "Dataset"
    ) +
    labs(
      title = "CV Fold Variability by Model and Dataset",
      x = "CV Fold",
      y = "RMSE (Mg/ha)"
    ) +
    theme(legend.position = "bottom")
  
  ggsave(output_path, p, width = 10, height = 6, dpi = 150)
  cli_alert_success("Saved: {output_path}")
}

# Generate plots
plot_rmse_comparison(cv_summary, fs::path(figures_dir, "rmse_comparison.png"))
plot_r2_comparison(cv_summary, fs::path(figures_dir, "r2_comparison.png"))
plot_delta_metrics(delta_metrics, fs::path(figures_dir, "delta_metrics.png"))
plot_fold_variability(all_results, fs::path(figures_dir, "fold_variability.png"))

# --- Scale-Stratified Analysis -----------------------------------------------

cli_h2("Scale-stratified analysis")

# Check if we have hex scale information
if ("hex_scale_ha" %in% names(fia_cv) && "hex_scale_ha" %in% names(combined_cv)) {
  
  #' Run CV stratified by hex scale
  #' 
  #' @param data Data with fold and hex_scale_ha
  #' @param response Response variable
  #' @param model_types Models to fit
  #' @param dataset_name Dataset identifier
  #' @return Stratified results
  run_cv_by_scale <- function(data, response, model_types, include_climate, dataset_name) {
    
    scales <- sort(unique(data$hex_scale_ha))
    results_list <- list()
    
    for (scale in scales) {
      scale_data <- data %>% filter(hex_scale_ha == scale)
      
      if (nrow(scale_data) < 50) {
        cli_alert_warning("Skipping scale {scale}ha: only {nrow(scale_data)} plots")
        next
      }
      
      # Recreate CV folds for this scale
      scale_cv <- create_hex_cv_folds(scale_data, k = min(cv_folds, 3))
      
      scale_results <- run_cv(scale_cv, response, model_types, include_climate, dataset_name) %>%
        mutate(hex_scale_ha = scale)
      
      results_list[[length(results_list) + 1]] <- scale_results
    }
    
    bind_rows(results_list)
  }
  
  cli_alert("Running scale-stratified CV (this may take a while)...")
  
  # Run for both datasets
  fia_scale_results <- run_cv_by_scale(fia_cv, response_var, "lm", include_climate, "FIA-only")
  combined_scale_results <- run_cv_by_scale(combined_cv, response_var, "lm", include_climate, "FIA+NEFIN")
  
  scale_results <- bind_rows(fia_scale_results, combined_scale_results)
  
  if (nrow(scale_results) > 0) {
    # Summarize by scale
    scale_summary <- scale_results %>%
      group_by(dataset, hex_scale_ha) %>%
      summarise(
        rmse_mean = mean(rmse, na.rm = TRUE),
        r2_mean = mean(r2, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Plot scale analysis
    p_scale <- ggplot(scale_summary, aes(x = factor(hex_scale_ha), y = rmse_mean, 
                                          fill = dataset)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      scale_fill_manual(
        values = c("FIA-only" = "#E69F00", "FIA+NEFIN" = "#56B4E9"),
        name = "Dataset"
      ) +
      labs(
        title = "RMSE by Hex Scale: FIA-only vs FIA+NEFIN",
        x = "Hex Scale (ha)",
        y = "RMSE (Mg/ha)"
      ) +
      theme(legend.position = "bottom")
    
    ggsave(fs::path(figures_dir, "rmse_by_scale.png"), p_scale, width = 10, height = 6, dpi = 150)
    cli_alert_success("Saved scale analysis plot")
    
    # Save scale results
    write_csv(scale_results, fs::path(output_dir, "ndvi_model_scale_results.csv"))
  }
  
} else {
  cli_alert_info("Hex scale column not found - skipping scale-stratified analysis")
}

# --- Save Results ------------------------------------------------------------

cli_h2("Saving results")

# Save full CV results
write_csv(all_results, fs::path(output_dir, "ndvi_model_cv_results.csv"))
cli_alert_success("Saved: {fs::path(output_dir, 'ndvi_model_cv_results.csv')}")

# Save summary
write_csv(cv_summary, fs::path(output_dir, "ndvi_model_summary.csv"))
cli_alert_success("Saved: {fs::path(output_dir, 'ndvi_model_summary.csv')}")

# Save delta metrics
write_csv(delta_metrics, fs::path(output_dir, "ndvi_model_delta_metrics.csv"))
cli_alert_success("Saved: {fs::path(output_dir, 'ndvi_model_delta_metrics.csv')}")

# Save metadata YAML
metadata <- list(
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  cv_method = cv_method,
  cv_folds = cv_folds,
  models_fit = models_to_fit,
  response_variable = response_var,
  include_climate = include_climate,
  n_fia_plots = nrow(fia_model_data),
  n_combined_plots = nrow(combined_model_data),
  summary = list(
    fia_best_model = cv_summary %>% 
      filter(dataset == "FIA-only") %>% 
      slice_min(rmse_mean) %>% 
      pull(model),
    combined_best_model = cv_summary %>% 
      filter(dataset == "FIA+NEFIN") %>% 
      slice_min(rmse_mean) %>% 
      pull(model),
    nefin_improves_rmse = any(delta_metrics$delta_rmse < 0),
    nefin_improves_r2 = any(delta_metrics$delta_r2 > 0)
  )
)

write_yaml(metadata, fs::path(output_dir, "ndvi_model_metadata.yml"))
cli_alert_success("Saved: {fs::path(output_dir, 'ndvi_model_metadata.yml')}")

# --- Final Summary -----------------------------------------------------------

cli_h1("NDVI Model Comparison Complete")

cli_alert_success("Outputs saved to: {output_dir}")
cli_alert_info("CV Results: ndvi_model_cv_results.csv")
cli_alert_info("Summary: ndvi_model_summary.csv")
cli_alert_info("Delta Metrics: ndvi_model_delta_metrics.csv")
cli_alert_info("Figures: {figures_dir}")

# Print key findings
cli_h3("Key Findings")

best_model <- cv_summary %>%
  slice_min(rmse_mean) %>%
  slice(1)

cli_alert("Best overall model: {best_model$model} ({best_model$dataset})")
cli_alert("  RMSE: {round(best_model$rmse_mean, 2)} ± {round(best_model$rmse_sd, 2)} Mg/ha")
cli_alert("  R²: {round(best_model$r2_mean * 100, 1)}%")

if (any(delta_metrics$delta_rmse < 0, na.rm = TRUE)) {
  improved_models <- delta_metrics %>% filter(delta_rmse < 0) %>% pull(model)
  cli_alert_success("NEFIN improved RMSE for: {paste(improved_models, collapse = ', ')}")
} else {
  cli_alert_warning("NEFIN did not improve RMSE for any model")
}

if (any(delta_metrics$delta_r2 > 0, na.rm = TRUE)) {
  improved_models <- delta_metrics %>% filter(delta_r2 > 0) %>% pull(model)
  cli_alert_success("NEFIN improved R² for: {paste(improved_models, collapse = ', ')}")
} else {
  cli_alert_warning("NEFIN did not improve R² for any model")
}

cli_alert_info("Run ndvi_model_report.R for detailed interpretation")
