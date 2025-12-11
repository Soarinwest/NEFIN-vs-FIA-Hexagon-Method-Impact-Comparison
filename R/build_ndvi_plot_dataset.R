# R/build_ndvi_plot_dataset.R
# ============================================================================
# Purpose: Assemble unified plot-level dataset combining:
#          - FIA biomass (AGLB) from plot-level computations
#          - NDVI values (MODIS and/or S2) from extraction
#          - PRISM climate variables (if available)
#          - Hex assignments at multiple scales
#
# Creates datasets for:
#   1. FIA-only (fuzzed) - uses jitter replicate means or all replicates
#   2. NEFIN-only (true coords)
#   3. FIA+NEFIN combined
#
# Inputs:
#   - data/processed/ndvi/fia_ndvi_extracted.csv
#   - data/processed/ndvi/nefin_ndvi_extracted.csv
#   - data/interim/fia_region/tree.csv, plot.csv (for AGLB computation)
#   - data/processed/nefin_processed.csv (for NEFIN AGLB)
#   - data/processed/hex_prism_values.csv (optional climate)
#   - configs/process.yml
#
# Outputs:
#   - data/processed/ndvi/modeling_dataset_fia.csv
#   - data/processed/ndvi/modeling_dataset_nefin.csv
#   - data/processed/ndvi/modeling_dataset_combined.csv
#   - data/processed/ndvi/dataset_summary.yml
#
# Usage:
#   Rscript R/build_ndvi_plot_dataset.R [--overwrite] [--aggregate-reps]
#
# Pipeline integration:
#   Run after: R/extract_ndvi_to_plots.R
#   Run before: R/ndvi_model_comparison.R
# ============================================================================

suppressPackageStartupMessages({

  library(dplyr)
  library(readr)
  library(fs)
  library(yaml)
  library(tidyr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Try to source utils_metrics.R, provide fallback if not found
if (fs::file_exists("R/utils_metrics.R")) {
  source("R/utils_metrics.R")
} else {
  message("Note: R/utils_metrics.R not found, using inline AGLB computation")
  
  # Inline fallback for build_plot_aglb
  build_plot_aglb <- function(tree_path, plot_path) {
    # Load data
    tree <- readr::read_csv(tree_path, show_col_types = FALSE)
    plot <- readr::read_csv(plot_path, show_col_types = FALSE)
    
    # Compute plot-level AGLB from tree data
    # DRYBIO_AG is in lbs, TPA_UNADJ is trees per acre
    # Convert to Mg/ha: lbs * 0.000453592 * trees_per_acre * 2.471054
    tree_aglb <- tree |>
      dplyr::filter(!is.na(DRYBIO_AG), !is.na(TPA_UNADJ), STATUSCD == 1) |>
      dplyr::mutate(
        # Biomass contribution per tree in Mg/ha
        aglb_contrib = DRYBIO_AG * 0.000453592 * TPA_UNADJ * 2.471054
      ) |>
      dplyr::group_by(PLT_CN) |>
      dplyr::summarise(
        aglb_Mg_per_ha = sum(aglb_contrib, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Join with plot metadata
    plot_aglb <- plot |>
      dplyr::select(CN, STATECD, MEASYEAR) |>
      dplyr::left_join(tree_aglb, by = c("CN" = "PLT_CN"))
    
    plot_aglb
  }
}

# ============================================================================
# Build FIA AGLB at plot level
# ============================================================================

#' Compute plot-level AGLB from FIA tree data
#' 
#' @param project_dir Project root
#' @return data.frame with CN, STATECD, MEASYEAR, aglb_Mg_per_ha
compute_fia_plot_aglb <- function(project_dir = ".") {
  
  fia_root <- if (fs::dir_exists(fs::path(project_dir, "data", "interim", "fia_region"))) {
    fs::path(project_dir, "data", "interim", "fia_region")
  } else {
    fs::path(project_dir, "data", "interim", "fia")
  }
  
  tree_csv <- fs::path(fia_root, "tree.csv")
  plot_csv <- fs::path(fia_root, "plot.csv")
  
  if (!fs::file_exists(tree_csv) || !fs::file_exists(plot_csv)) {
    stop("FIA data not found. Expected:\n  ", tree_csv, "\n  ", plot_csv)
  }
  
  # Use existing utility function
  plot_aglb <- build_plot_aglb(tree_csv, plot_csv)
  
  plot_aglb |>
    dplyr::select(CN, STATECD, MEASYEAR, aglb_Mg_per_ha) |>
    dplyr::filter(!is.na(aglb_Mg_per_ha), is.finite(aglb_Mg_per_ha))
}

# ============================================================================
# Join NDVI with AGLB
# ============================================================================

#' Build FIA modeling dataset
#' 
#' @param fia_ndvi_path Path to FIA NDVI extraction
#' @param fia_aglb data.frame with FIA plot-level AGLB
#' @param aggregate_reps If TRUE, average NDVI across jitter replicates
#' @return data.frame ready for modeling
build_fia_dataset <- function(fia_ndvi_path, fia_aglb, aggregate_reps = TRUE) {
  
  message("→ Building FIA modeling dataset...")
  
  if (!fs::file_exists(fia_ndvi_path)) {
    stop("FIA NDVI extraction not found: ", fia_ndvi_path)
  }
  
  fia_ndvi <- readr::read_csv(fia_ndvi_path, show_col_types = FALSE)
  
  message("  NDVI records: ", format(nrow(fia_ndvi), big.mark = ","))
  
  # Join with AGLB
  fia_joined <- fia_ndvi |>
    dplyr::left_join(fia_aglb, by = c("CN", "STATECD", "MEASYEAR"))
  
  message("  After AGLB join: ", format(nrow(fia_joined), big.mark = ","))
  message("  With valid AGLB: ", sum(!is.na(fia_joined$aglb_Mg_per_ha)))
  
  if (aggregate_reps) {
    # Average NDVI across jitter replicates for each plot
    message("  Aggregating across jitter replicates...")
    
    # Get hex columns
    hex_cols <- names(fia_joined)[grepl("^hex_id_", names(fia_joined))]
    
    # For hex assignments, take the mode (most common assignment)
    get_mode <- function(x) {
      x <- x[!is.na(x)]
      if (length(x) == 0) return(NA_character_)
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    
    fia_agg <- fia_joined |>
      dplyr::group_by(CN, STATECD, MEASYEAR, aglb_Mg_per_ha) |>
      dplyr::summarise(
        # NDVI: mean and SD across replicates
        ndvi_modis_mean = mean(ndvi_modis, na.rm = TRUE),
        ndvi_modis_sd = sd(ndvi_modis, na.rm = TRUE),
        ndvi_s2_mean = mean(ndvi_s2, na.rm = TRUE),
        ndvi_s2_sd = sd(ndvi_s2, na.rm = TRUE),
        
        # Count valid replicates
        n_reps_modis = sum(!is.na(ndvi_modis)),
        n_reps_s2 = sum(!is.na(ndvi_s2)),
        n_reps_total = dplyr::n(),
        
        # Window (should be same for all reps)
        window_name = dplyr::first(window_name),
        
        # Hex assignments (mode across replicates)
        dplyr::across(dplyr::all_of(hex_cols), get_mode),
        
        .groups = "drop"
      )
    
    # Rename NDVI columns for consistency
    fia_agg <- fia_agg |>
      dplyr::rename(
        ndvi_modis = ndvi_modis_mean,
        ndvi_s2 = ndvi_s2_mean
      )
    
    message("  Aggregated to ", nrow(fia_agg), " unique plots")
    
    return(fia_agg)
    
  } else {
    # Return all replicates (for uncertainty analysis)
    message("  Keeping all ", nrow(fia_joined), " plot-replicate records")
    return(fia_joined)
  }
}

#' Build NEFIN modeling dataset
#' 
#' @param nefin_ndvi_path Path to NEFIN NDVI extraction
#' @param nefin_hex_path Path to NEFIN hex assignments
#' @return data.frame ready for modeling
build_nefin_dataset <- function(nefin_ndvi_path, 
                                 nefin_hex_path = "data/processed/nefin_hex_assignments.csv") {
  
  message("→ Building NEFIN modeling dataset...")
  
  if (!fs::file_exists(nefin_ndvi_path)) {
    message("  ⚠ NEFIN NDVI extraction not found")
    return(NULL)
  }
  
  nefin_ndvi <- readr::read_csv(nefin_ndvi_path, show_col_types = FALSE)
  
  message("  NDVI records: ", nrow(nefin_ndvi))
  
  # NEFIN AGLB should already be in the NDVI extraction (from nefin_processed.csv)
  if (!"aglb_Mg_per_ha" %in% names(nefin_ndvi)) {
    # Try to load from processed file
    nefin_processed <- fs::path("data", "processed", "nefin_processed.csv")
    if (fs::file_exists(nefin_processed)) {
      nefin_aglb <- readr::read_csv(nefin_processed, show_col_types = FALSE) |>
        dplyr::select(CN, aglb_Mg_per_ha)
      
      nefin_ndvi <- nefin_ndvi |>
        dplyr::left_join(nefin_aglb, by = "CN")
    }
  }
  
  # Load hex assignments if not in NDVI data
  hex_cols <- names(nefin_ndvi)[grepl("^hex_id_", names(nefin_ndvi))]
  
  if (length(hex_cols) == 0 && fs::file_exists(nefin_hex_path)) {
    message("  Loading hex assignments from: ", nefin_hex_path)
    nefin_hex <- readr::read_csv(nefin_hex_path, show_col_types = FALSE)
    
    # Join on CN
    nefin_ndvi <- nefin_ndvi |>
      dplyr::left_join(
        nefin_hex |> dplyr::select(CN, dplyr::starts_with("hex_id_")),
        by = "CN"
      )
  }
  
  # Select and rename columns for consistency
  nefin_out <- nefin_ndvi |>
    dplyr::select(
      CN, STATECD, MEASYEAR,
      aglb_Mg_per_ha,
      ndvi_modis, ndvi_s2,
      window_name,
      dplyr::starts_with("hex_id_")
    ) |>
    dplyr::mutate(
      source = "NEFIN",
      # NEFIN has no replicate uncertainty for NDVI
      ndvi_modis_sd = NA_real_,
      ndvi_s2_sd = NA_real_,
      n_reps_modis = 1L,
      n_reps_s2 = 1L,
      n_reps_total = 1L
    )
  
  message("  NEFIN dataset: ", nrow(nefin_out), " plots")
  message("  With valid AGLB: ", sum(!is.na(nefin_out$aglb_Mg_per_ha)))
  
  nefin_out
}

#' Add PRISM climate covariates to dataset
#' 
#' @param df data.frame with hex_id columns
#' @param prism_path Path to hex-level PRISM data
#' @param hex_scale Which hex scale to use for climate lookup
#' @return data.frame with climate columns added
add_climate_covariates <- function(df, 
                                    prism_path = "data/processed/hex_prism_values.csv",
                                    hex_scale = NULL) {
  
  if (!fs::file_exists(prism_path)) {
    message("  ℹ PRISM climate data not found, skipping")
    return(df)
  }
  
  message("  Adding PRISM climate covariates...")
  
  prism_data <- readr::read_csv(prism_path, show_col_types = FALSE)
  
  # Determine which hex scale to use
  hex_cols <- names(df)[grepl("^hex_id_", names(df))]
  
  if (is.null(hex_scale) && length(hex_cols) > 0) {
    # Use the finest scale available
    hex_scale <- sub("^hex_id_", "", hex_cols[1])
  }
  
  if (is.null(hex_scale)) {
    message("  ⚠ No hex assignments found, skipping climate")
    return(df)
  }
  
  hex_col <- paste0("hex_id_", hex_scale)
  
  if (!(hex_col %in% names(df))) {
    message("  ⚠ Hex column not found: ", hex_col)
    return(df)
  }
  
  # Get climate columns for this scale
  climate_cols <- names(prism_data)[grepl(paste0("_", hex_scale, "$"), names(prism_data))]
  
  if (length(climate_cols) == 0) {
    message("  ⚠ No PRISM data for scale: ", hex_scale)
    return(df)
  }
  
  # Select and rename climate columns
  prism_subset <- prism_data |>
    dplyr::select(hex_id, dplyr::all_of(climate_cols))
  
  # Remove scale suffix from column names
  names(prism_subset) <- gsub(paste0("_", hex_scale, "$"), "", names(prism_subset))
  
  # Join
  df_with_climate <- df |>
    dplyr::left_join(prism_subset, by = setNames("hex_id", hex_col))
  
  n_climate <- sum(names(prism_subset) != "hex_id")
  message("  Added ", n_climate, " climate variables")
  
  df_with_climate
}

# ============================================================================
# Main Entry Point
# ============================================================================

#' Build all modeling datasets
#' 
#' @param project_dir Project root
#' @param overwrite Overwrite existing outputs
#' @param aggregate_reps Aggregate FIA across jitter replicates
#' @param hex_scale Hex scale for climate lookup (NULL = auto)
build_ndvi_plot_dataset <- function(project_dir = ".",
                                     overwrite = FALSE,
                                     aggregate_reps = TRUE,
                                     hex_scale = NULL) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Build NDVI Plot-Level Modeling Dataset                  ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Paths
  ndvi_dir <- fs::path(project_dir, "data", "processed", "ndvi")
  
  fia_ndvi_path <- fs::path(ndvi_dir, "fia_ndvi_extracted.csv")
  nefin_ndvi_path <- fs::path(ndvi_dir, "nefin_ndvi_extracted.csv")
  
  fia_out <- fs::path(ndvi_dir, "modeling_dataset_fia.csv")
  nefin_out <- fs::path(ndvi_dir, "modeling_dataset_nefin.csv")
  combined_out <- fs::path(ndvi_dir, "modeling_dataset_combined.csv")
  summary_out <- fs::path(ndvi_dir, "dataset_summary.yml")
  
  # Check if already done
  if (!overwrite && fs::file_exists(fia_out) && fs::file_exists(summary_out)) {
    message("✓ Modeling datasets already exist")
    message("  FIA: ", fia_out)
    if (fs::file_exists(nefin_out)) message("  NEFIN: ", nefin_out)
    if (fs::file_exists(combined_out)) message("  Combined: ", combined_out)
    message("  Use --overwrite to regenerate")
    return(invisible(NULL))
  }
  
  # Check for NDVI extraction
  if (!fs::file_exists(fia_ndvi_path)) {
    stop("FIA NDVI extraction not found: ", fia_ndvi_path,
         "\n  Run: Rscript R/extract_ndvi_to_plots.R")
  }
  
  # Compute FIA AGLB
  message("→ Computing FIA plot-level AGLB...")
  fia_aglb <- compute_fia_plot_aglb(project_dir)
  message("  FIA plots with AGLB: ", nrow(fia_aglb))
  
  # Build FIA dataset
  fia_dataset <- build_fia_dataset(fia_ndvi_path, fia_aglb, aggregate_reps)
  fia_dataset$source <- "FIA"
  
  # Add climate
  fia_dataset <- add_climate_covariates(fia_dataset, hex_scale = hex_scale)
  
  # Write FIA dataset
  readr::write_csv(fia_dataset, fia_out)
  message("\n✓ Wrote FIA dataset: ", fia_out)
  
  # Build NEFIN dataset
  nefin_dataset <- NULL
  if (fs::file_exists(nefin_ndvi_path)) {
    nefin_dataset <- build_nefin_dataset(nefin_ndvi_path)
    
    if (!is.null(nefin_dataset)) {
      nefin_dataset <- add_climate_covariates(nefin_dataset, hex_scale = hex_scale)
      readr::write_csv(nefin_dataset, nefin_out)
      message("✓ Wrote NEFIN dataset: ", nefin_out)
    }
  }
  
  # Build combined dataset
  if (!is.null(nefin_dataset)) {
    message("\n→ Building combined FIA+NEFIN dataset...")
    
    # Ensure column consistency
    common_cols <- intersect(names(fia_dataset), names(nefin_dataset))
    
    combined_dataset <- dplyr::bind_rows(
      fia_dataset |> dplyr::select(dplyr::all_of(common_cols)),
      nefin_dataset |> dplyr::select(dplyr::all_of(common_cols))
    )
    
    readr::write_csv(combined_dataset, combined_out)
    message("✓ Wrote combined dataset: ", combined_out)
  }
  
  # Summary statistics
  summary_stats <- list(
    created = as.character(Sys.time()),
    aggregate_reps = aggregate_reps,
    fia = list(
      path = as.character(fia_out),
      n_plots = nrow(fia_dataset),
      n_with_aglb = sum(!is.na(fia_dataset$aglb_Mg_per_ha)),
      n_with_modis = sum(!is.na(fia_dataset$ndvi_modis)),
      n_with_s2 = sum(!is.na(fia_dataset$ndvi_s2)),
      aglb_mean = mean(fia_dataset$aglb_Mg_per_ha, na.rm = TRUE),
      aglb_sd = sd(fia_dataset$aglb_Mg_per_ha, na.rm = TRUE),
      ndvi_modis_mean = mean(fia_dataset$ndvi_modis, na.rm = TRUE),
      ndvi_s2_mean = mean(fia_dataset$ndvi_s2, na.rm = TRUE)
    )
  )
  
  if (!is.null(nefin_dataset)) {
    summary_stats$nefin <- list(
      path = as.character(nefin_out),
      n_plots = nrow(nefin_dataset),
      n_with_aglb = sum(!is.na(nefin_dataset$aglb_Mg_per_ha)),
      n_with_modis = sum(!is.na(nefin_dataset$ndvi_modis)),
      n_with_s2 = sum(!is.na(nefin_dataset$ndvi_s2)),
      aglb_mean = mean(nefin_dataset$aglb_Mg_per_ha, na.rm = TRUE),
      aglb_sd = sd(nefin_dataset$aglb_Mg_per_ha, na.rm = TRUE),
      ndvi_modis_mean = mean(nefin_dataset$ndvi_modis, na.rm = TRUE),
      ndvi_s2_mean = mean(nefin_dataset$ndvi_s2, na.rm = TRUE)
    )
    
    summary_stats$combined <- list(
      path = as.character(combined_out),
      n_plots = nrow(combined_dataset),
      n_fia = sum(combined_dataset$source == "FIA"),
      n_nefin = sum(combined_dataset$source == "NEFIN")
    )
  }
  
  yaml::write_yaml(summary_stats, summary_out)
  message("\n✓ Wrote summary: ", summary_out)
  
  # Print summary
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("DATASET SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("\n")
  
  cat("FIA (fuzzed coordinates):\n")
  cat("  Plots: ", format(summary_stats$fia$n_plots, big.mark = ","), "\n")
  cat("  With AGLB: ", format(summary_stats$fia$n_with_aglb, big.mark = ","), "\n")
  cat("  With MODIS NDVI: ", format(summary_stats$fia$n_with_modis, big.mark = ","), "\n")
  cat("  With S2 NDVI: ", format(summary_stats$fia$n_with_s2, big.mark = ","), "\n")
  cat("  Mean AGLB: ", round(summary_stats$fia$aglb_mean, 2), " Mg/ha\n")
  
  if (!is.null(nefin_dataset)) {
    cat("\nNEFIN (true coordinates):\n")
    cat("  Plots: ", format(summary_stats$nefin$n_plots, big.mark = ","), "\n")
    cat("  With AGLB: ", format(summary_stats$nefin$n_with_aglb, big.mark = ","), "\n")
    cat("  With MODIS NDVI: ", format(summary_stats$nefin$n_with_modis, big.mark = ","), "\n")
    cat("  With S2 NDVI: ", format(summary_stats$nefin$n_with_s2, big.mark = ","), "\n")
    cat("  Mean AGLB: ", round(summary_stats$nefin$aglb_mean, 2), " Mg/ha\n")
    
    cat("\nCombined:\n")
    cat("  Total plots: ", format(summary_stats$combined$n_plots, big.mark = ","), "\n")
    cat("  FIA: ", format(summary_stats$combined$n_fia, big.mark = ","), "\n")
    cat("  NEFIN: ", format(summary_stats$combined$n_nefin, big.mark = ","), "\n")
  }
  
  cat("\n")
  cat("Next step:\n")
  cat("  Rscript R/ndvi_model_comparison.R\n")
  cat("\n")
  
  invisible(list(
    fia = fia_out,
    nefin = nefin_out,
    combined = combined_out
  ))
}

# ============================================================================
# CLI Entry Point
# ============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  overwrite <- "--overwrite" %in% args
  aggregate_reps <- !("--all-reps" %in% args)  # Default to aggregating
  
  build_ndvi_plot_dataset(
    project_dir = ".",
    overwrite = overwrite,
    aggregate_reps = aggregate_reps
  )
}
