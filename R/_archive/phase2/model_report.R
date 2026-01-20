#!/usr/bin/env Rscript
# ==============================================================================
# ndvi_model_report.R
# ==============================================================================
# Purpose:
#   Generate interpretable text summaries and scale guidance from NDVI model
#   comparison results. Creates markdown reports similar to the existing
#   fuzz_vs_uncertainty scale guidance analysis.
#
# Inputs:
#   - runs/consolidated_.../ndvi_model_cv_results.csv
#   - runs/consolidated_.../ndvi_model_summary.csv
#   - runs/consolidated_.../ndvi_model_delta_metrics.csv
#   - runs/consolidated_.../ndvi_model_scale_results.csv (if available)
#   - runs/consolidated_.../ndvi_model_metadata.yml
#
# Outputs:
#   - runs/consolidated_.../ndvi_model_report.md (comprehensive markdown report)
#   - runs/consolidated_.../ndvi_scale_guidance.md (scale selection guidance)
#
# Usage:
#   Rscript R/ndvi_model_report.R [--output-dir=path]
#
# Pipeline integration:
#   Run after ndvi_model_comparison.R as final reporting step
#
# Author: Claude (Anthropic)
# Date: 2025-12
# ==============================================================================

# --- Setup -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(yaml)
  library(fs)
  library(glue)
  library(cli)
})

# --- Parse CLI Arguments -----------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

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
      output_dir <- run_dirs[length(run_dirs)]
    } else {
      cli_abort("No consolidated run directories found")
    }
  } else {
    cli_abort("Runs directory not found")
  }
}

cli_h1("NDVI Model Report Generator")
cli_alert_info("Reading results from: {output_dir}")

# --- Load Results ------------------------------------------------------------

# Required files
cv_results_path <- fs::path(output_dir, "ndvi_model_cv_results.csv")
summary_path <- fs::path(output_dir, "ndvi_model_summary.csv")
delta_path <- fs::path(output_dir, "ndvi_model_delta_metrics.csv")
metadata_path <- fs::path(output_dir, "ndvi_model_metadata.yml")

# Optional files
scale_results_path <- fs::path(output_dir, "ndvi_model_scale_results.csv")

# Check required files exist
required_files <- c(cv_results_path, summary_path, delta_path, metadata_path)
missing <- required_files[!fs::file_exists(required_files)]
if (length(missing) > 0) {
  cli_abort("Missing required files: {paste(missing, collapse = ', ')}")
}

# Load data
cv_results <- read_csv(cv_results_path, show_col_types = FALSE)
summary_df <- read_csv(summary_path, show_col_types = FALSE)
delta_df <- read_csv(delta_path, show_col_types = FALSE)
metadata <- read_yaml(metadata_path)

# Load optional scale results
scale_results <- NULL
if (fs::file_exists(scale_results_path)) {
  scale_results <- read_csv(scale_results_path, show_col_types = FALSE)
}

cli_alert_success("Loaded all result files")

# --- Helper Functions --------------------------------------------------------

#' Format number with appropriate precision
fmt <- function(x, digits = 2) {
  format(round(x, digits), nsmall = digits)
}

#' Format percentage
fmt_pct <- function(x, digits = 1) {
  paste0(format(round(x * 100, digits), nsmall = digits), "%")
}

#' Generate interpretation of RMSE change
interpret_rmse_change <- function(delta, baseline, threshold_meaningful = 0.05) {
  pct_change <- delta / baseline
  
  if (is.na(delta)) return("Unable to calculate change")
  
  if (abs(pct_change) < threshold_meaningful) {
    "negligible change (< 5%)"
  } else if (delta < 0) {
    glue("improvement of {fmt(abs(delta))} Mg/ha ({fmt_pct(abs(pct_change))} reduction)")
  } else {
    glue("degradation of {fmt(delta)} Mg/ha ({fmt_pct(pct_change)} increase)")
  }
}

#' Generate interpretation of R² change
interpret_r2_change <- function(delta, baseline, threshold_meaningful = 0.02) {
  
  if (is.na(delta)) return("Unable to calculate change")
  
  if (abs(delta) < threshold_meaningful) {
    "negligible change (< 2 percentage points)"
  } else if (delta > 0) {
    glue("improvement of {fmt_pct(delta)} explained variance")
  } else {
    glue("reduction of {fmt_pct(abs(delta))} explained variance")
  }
}

# --- Generate Main Report ----------------------------------------------------

cli_h2("Generating main report")

report_lines <- c(
  "# NDVI Spatial Model Comparison Report",
  "",
  glue("**Generated:** {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}"),
  "",
  glue("**Analysis Directory:** `{output_dir}`"),
  "",
  "---",
  "",
  "## Executive Summary",
  "",
  "This report compares spatial predictive models for Above-Ground Live Biomass (AGLB)",
  "using two datasets:",
  "",
  "1. **FIA-only**: Forest Inventory and Analysis plots with privacy-fuzzed coordinates (±1.6 km jitter)",
  "2. **FIA+NEFIN**: Combined dataset including high-fidelity NEFIN plot locations (true coordinates)",
  "",
  "The analysis quantifies whether adding NEFIN's precise coordinates improves NDVI-based",
  "biomass predictions compared to FIA alone.",
  "",
  "---",
  "",
  "## Methods",
  "",
  glue("- **Cross-validation method:** {metadata$cv_method}"),
  glue("- **Number of folds:** {metadata$cv_folds}"),
  glue("- **Models evaluated:** {paste(metadata$models_fit, collapse = ', ')}"),
  glue("- **Response variable:** {metadata$response_variable}"),
  glue("- **Climate covariates included:** {metadata$include_climate}"),
  "",
  glue("- **FIA plots:** {format(metadata$n_fia_plots, big.mark = ',')}"),
  glue("- **Combined plots (FIA+NEFIN):** {format(metadata$n_combined_plots, big.mark = ',')}"),
  "",
  "---",
  "",
  "## Results",
  "",
  "### Model Performance Summary",
  "",
  "| Dataset | Model | RMSE (Mg/ha) | R² | MAE (Mg/ha) |",
  "|---------|-------|--------------|-----|-------------|"
)

# Add summary rows
for (i in 1:nrow(summary_df)) {
  row <- summary_df[i, ]
  report_lines <- c(report_lines,
    glue("| {row$dataset} | {row$model} | {fmt(row$rmse_mean)} ± {fmt(row$rmse_sd)} | {fmt_pct(row$r2_mean)} | {fmt(row$mae_mean)} |")
  )
}

report_lines <- c(report_lines,
  "",
  "### Impact of Adding NEFIN (Delta Metrics)",
  "",
  "Positive values for ΔRMSE indicate FIA+NEFIN performed worse; negative values indicate improvement.",
  "Positive values for ΔR² indicate improvement.",
  "",
  "| Model | ΔRMSE (Mg/ha) | ΔR² | RMSE % Change | Interpretation |",
  "|-------|---------------|-----|---------------|----------------|"
)

# Add delta rows with interpretation
for (i in 1:nrow(delta_df)) {
  row <- delta_df[i, ]
  
  # Get baseline RMSE for interpretation
  baseline_rmse <- summary_df %>%
    filter(dataset == "FIA-only", model == row$model) %>%
    pull(rmse_mean)
  
  interp <- if (!is.na(row$delta_rmse) && row$delta_rmse < 0) {
    "✓ NEFIN helps"
  } else if (!is.na(row$delta_rmse) && abs(row$rmse_pct_change) < 5) {
    "≈ No significant change"
  } else {
    "✗ NEFIN may introduce bias"
  }
  
  report_lines <- c(report_lines,
    glue("| {row$model} | {fmt(row$delta_rmse)} | {fmt(row$delta_r2, 3)} | {fmt(row$rmse_pct_change)}% | {interp} |")
  )
}

# Add interpretation section
report_lines <- c(report_lines,
  "",
  "### Interpretation",
  ""
)

# Determine overall finding
nefin_helps_rmse <- any(delta_df$delta_rmse < 0, na.rm = TRUE)
nefin_helps_r2 <- any(delta_df$delta_r2 > 0, na.rm = TRUE)

if (nefin_helps_rmse && nefin_helps_r2) {
  report_lines <- c(report_lines,
    "**Overall Finding:** Adding NEFIN plots with true coordinates **improves** spatial",
    "prediction accuracy for at least some model types. This suggests that the high-fidelity",
    "coordinates provide valuable information for NDVI-based biomass estimation.",
    ""
  )
} else if (!nefin_helps_rmse && !nefin_helps_r2) {
  report_lines <- c(report_lines,
    "**Overall Finding:** Adding NEFIN plots did **not improve** prediction accuracy",
    "across the tested models. This may indicate:",
    "",
    "1. FIA fuzzing has minimal impact on NDVI-AGLB relationships at this scale",
    "2. NEFIN's ownership-biased sampling introduces systematic differences",
    "3. The NDVI products (MODIS 250m, Sentinel-2 10m) are robust to coordinate uncertainty",
    ""
  )
} else {
  report_lines <- c(report_lines,
    "**Overall Finding:** Results are **mixed** - NEFIN improves some metrics but not others.",
    "This nuanced outcome suggests the value of high-fidelity coordinates depends on the",
    "specific model and evaluation metric used.",
    ""
  )
}

# Best performing model
best_fia <- summary_df %>% filter(dataset == "FIA-only") %>% slice_min(rmse_mean)
best_combined <- summary_df %>% filter(dataset == "FIA+NEFIN") %>% slice_min(rmse_mean)

report_lines <- c(report_lines,
  "**Best Performing Models:**",
  "",
  glue("- FIA-only: **{best_fia$model}** (RMSE = {fmt(best_fia$rmse_mean)} Mg/ha, R² = {fmt_pct(best_fia$r2_mean)})"),
  glue("- FIA+NEFIN: **{best_combined$model}** (RMSE = {fmt(best_combined$rmse_mean)} Mg/ha, R² = {fmt_pct(best_combined$r2_mean)})"),
  ""
)

# Add scale analysis if available
if (!is.null(scale_results) && nrow(scale_results) > 0) {
  report_lines <- c(report_lines,
    "---",
    "",
    "## Scale-Stratified Analysis",
    "",
    "Performance varies by hexagon scale. Smaller hexes capture finer spatial variation",
    "but may have fewer plots per unit. Larger hexes provide more stable estimates but",
    "smooth over local heterogeneity.",
    "",
    "| Scale (ha) | Dataset | RMSE (Mg/ha) | R² |",
    "|------------|---------|--------------|-----|"
  )
  
  # Summarize by scale
  scale_summary <- scale_results %>%
    group_by(dataset, hex_scale_ha) %>%
    summarise(
      rmse_mean = mean(rmse, na.rm = TRUE),
      r2_mean = mean(r2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(hex_scale_ha, dataset)
  
  for (i in 1:nrow(scale_summary)) {
    row <- scale_summary[i, ]
    report_lines <- c(report_lines,
      glue("| {format(row$hex_scale_ha, big.mark = ',')} | {row$dataset} | {fmt(row$rmse_mean)} | {fmt_pct(row$r2_mean)} |")
    )
  }
  
  # Find optimal scale
  optimal_scale <- scale_summary %>%
    filter(dataset == "FIA+NEFIN") %>%
    slice_min(rmse_mean) %>%
    pull(hex_scale_ha)
  
  report_lines <- c(report_lines,
    "",
    glue("**Optimal scale for FIA+NEFIN:** {format(optimal_scale, big.mark = ',')} ha"),
    ""
  )
}

# Add recommendations
report_lines <- c(report_lines,
  "---",
  "",
  "## Recommendations",
  "",
  "Based on this analysis:",
  ""
)

if (nefin_helps_rmse) {
  report_lines <- c(report_lines,
    "1. **Include NEFIN data** in spatial modeling workflows where high accuracy is critical",
    "2. The added value is most apparent for [specific model types identified above]",
    ""
  )
} else {
  report_lines <- c(report_lines,
    "1. **FIA alone may be sufficient** for NDVI-based biomass predictions at regional scales",
    "2. Coordinate fuzzing appears to have limited impact on model performance",
    ""
  )
}

report_lines <- c(report_lines,
  "3. Consider model selection based on the specific use case:",
  "   - Simple interpretable models (LM) for transparent, reproducible analysis",
  "   - GAMs for capturing non-linear NDVI-biomass relationships",
  "   - Random Forest for maximum predictive accuracy (if interpretability is less critical)",
  "",
  "4. For operational mapping, validate with independent ground-truth data",
  "",
  "---",
  "",
  "## Technical Notes",
  "",
  "- NDVI values extracted from MODIS (250m) and/or Sentinel-2 (10m) composites",
  "- FIA coordinate uncertainty represented by 80 Monte Carlo jitter replicates",
  "- Spatial cross-validation prevents overfitting due to spatial autocorrelation",
  "- PRISM climate covariates (temperature, precipitation) included where configured",
  "",
  "---",
  "",
  glue("*Report generated by ndvi_model_report.R on {format(Sys.time(), '%Y-%m-%d')}*")
)

# Write main report
report_path <- fs::path(output_dir, "ndvi_model_report.md")
writeLines(report_lines, report_path)
cli_alert_success("Saved main report: {report_path}")

# --- Generate Scale Guidance Document ----------------------------------------

cli_h2("Generating scale guidance document")

guidance_lines <- c(
  "# NDVI Model Scale Guidance",
  "",
  glue("**Generated:** {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}"),
  "",
  "---",
  "",
  "## Purpose",
  "",
  "This document provides guidance on selecting appropriate hexagon scales for",
  "NDVI-based biomass modeling, based on cross-validation performance metrics.",
  "",
  "## Key Findings",
  ""
)

if (!is.null(scale_results) && nrow(scale_results) > 0) {
  
  # Calculate scale-specific deltas
  scale_pivot <- scale_results %>%
    group_by(hex_scale_ha, dataset) %>%
    summarise(rmse = mean(rmse, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = dataset, values_from = rmse)
  
  if (all(c("FIA-only", "FIA+NEFIN") %in% names(scale_pivot))) {
    scale_pivot <- scale_pivot %>%
      mutate(
        delta = `FIA+NEFIN` - `FIA-only`,
        nefin_helps = delta < 0
      )
    
    # Find scales where NEFIN helps
    helpful_scales <- scale_pivot %>% filter(nefin_helps) %>% pull(hex_scale_ha)
    unhelpful_scales <- scale_pivot %>% filter(!nefin_helps) %>% pull(hex_scale_ha)
    
    if (length(helpful_scales) > 0) {
      guidance_lines <- c(guidance_lines,
        glue("**NEFIN improves predictions at:** {paste(format(helpful_scales, big.mark = ','), collapse = ', ')} ha"),
        ""
      )
    }
    
    if (length(unhelpful_scales) > 0) {
      guidance_lines <- c(guidance_lines,
        glue("**NEFIN provides no benefit at:** {paste(format(unhelpful_scales, big.mark = ','), collapse = ', ')} ha"),
        ""
      )
    }
    
    guidance_lines <- c(guidance_lines,
      "## Scale Selection Recommendations",
      "",
      "### For Maximum Accuracy",
      ""
    )
    
    best_scale <- scale_pivot %>% slice_min(`FIA+NEFIN`) %>% pull(hex_scale_ha)
    guidance_lines <- c(guidance_lines,
      glue("Use **{format(best_scale, big.mark = ',')} ha** hexes with the FIA+NEFIN combined dataset."),
      ""
    )
    
    guidance_lines <- c(guidance_lines,
      "### For FIA-only Analysis",
      ""
    )
    
    best_fia_scale <- scale_pivot %>% slice_min(`FIA-only`) %>% pull(hex_scale_ha)
    guidance_lines <- c(guidance_lines,
      glue("Use **{format(best_fia_scale, big.mark = ',')} ha** hexes when only FIA data are available."),
      ""
    )
  }
  
} else {
  guidance_lines <- c(guidance_lines,
    "Scale-stratified results not available. Run ndvi_model_comparison.R with",
    "datasets that include hex_scale_ha to generate scale-specific guidance.",
    ""
  )
}

guidance_lines <- c(guidance_lines,
  "---",
  "",
  "## General Guidance",
  "",
  "| Scale Range | Typical Use Case | NEFIN Value |",
  "|-------------|------------------|-------------|",
  "| 100-1,000 ha | Local forest management | Highest - precise locations matter |",
  "| 1,000-10,000 ha | Landscape analysis | Moderate - depends on forest heterogeneity |",
  "| 10,000-100,000 ha | Regional summaries | Lower - spatial averaging reduces fuzz impact |",
  "",
  "---",
  "",
  glue("*Guidance generated by ndvi_model_report.R on {format(Sys.time(), '%Y-%m-%d')}*")
)

# Write scale guidance
guidance_path <- fs::path(output_dir, "ndvi_scale_guidance.md")
writeLines(guidance_lines, guidance_path)
cli_alert_success("Saved scale guidance: {guidance_path}")

# --- Final Summary -----------------------------------------------------------

cli_h1("Report Generation Complete")
cli_alert_success("Main report: {report_path}")
cli_alert_success("Scale guidance: {guidance_path}")
