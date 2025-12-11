#!/usr/bin/env Rscript
# ==============================================================================
# build_modeling_datasets.R
# ==============================================================================
# Purpose:
#   Assemble unified plot-level modeling datasets by combining:
#   - FIA AGLB (computed from tree data)
#   - NEFIN AGLB (from nefin_processed.csv)
#   - Extracted NDVI (MODIS, S2)
#   - Extracted PRISM climate (tmean, ppt)
#   - Hex assignments at multiple scales
#
# Creates four modeling datasets:
#   1. FIA fuzzed (single replicate or aggregated)
#   2. FIA mean across replicates (uncertainty quantified)
#   3. NEFIN only (true coordinates)
#   4. FIA + NEFIN combined
#
# Inputs:
#   - data/processed/climate_at_plots/fia_climate.csv
#   - data/processed/climate_at_plots/nefin_climate.csv
#   - data/processed/ndvi_at_plots/fia_ndvi.csv
#   - data/processed/ndvi_at_plots/nefin_ndvi.csv
#   - data/interim/fia_region/tree.csv, plot.csv
#   - data/processed/nefin_processed.csv
#   - data/processed/plot_hex_assignments.csv
#   - data/processed/nefin_hex_assignments.csv
#
# Outputs:
#   - data/processed/modeling/dataset_fia_fuzzed.csv
#   - data/processed/modeling/dataset_fia_mean.csv
#   - data/processed/modeling/dataset_nefin.csv
#   - data/processed/modeling/dataset_combined.csv
#   - data/processed/modeling/dataset_summary.yml
#
# Usage:
#   Rscript R/build_modeling_datasets.R [--overwrite]
#
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(fs)
  library(yaml)
  library(purrr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==============================================================================
# Compute FIA Plot-Level AGLB
# ==============================================================================

#' Compute plot-level AGLB from FIA tree data
#' 
#' Uses DRYBIO_AG (above-ground dry biomass in lbs) and TPA_UNADJ
#' Converts to Mg/ha
#' 
#' @param tree_path Path to tree.csv
#' @param plot_path Path to plot.csv
#' @return data.frame with CN, STATECD, MEASYEAR, aglb_Mg_per_ha
compute_fia_aglb <- function(tree_path, plot_path) {
  
  message("Computing FIA plot-level AGLB...")
  
  if (!fs::file_exists(tree_path)) {
    stop("Tree data not found: ", tree_path)
  }
  if (!fs::file_exists(plot_path)) {
    stop("Plot data not found: ", plot_path)
  }
  
  tree <- read_csv(tree_path, show_col_types = FALSE)
  plot <- read_csv(plot_path, show_col_types = FALSE)
  
  message("  Trees: ", format(nrow(tree), big.mark = ","))
  message("  Plots: ", format(nrow(plot), big.mark = ","))
  
  # Compute AGLB per plot
  # DRYBIO_AG is lbs, TPA_UNADJ is trees per acre
  # Convert: lbs * 0.000453592 (to Mg) * TPA * 2.471054 (acre to ha)
  
  plot_aglb <- tree %>%
    filter(!is.na(DRYBIO_AG), !is.na(TPA_UNADJ)) %>%
    filter(STATUSCD == 1) %>%  # Live trees only
    mutate(
      # Biomass in Mg/ha for this tree
      aglb_contrib = DRYBIO_AG * 0.000453592 * TPA_UNADJ * 2.471054
    ) %>%
    group_by(PLT_CN) %>%
    summarise(
      aglb_Mg_per_ha = sum(aglb_contrib, na.rm = TRUE),
      n_trees = n(),
      .groups = "drop"
    )
  
  # Join with plot metadata
  result <- plot %>%
    select(CN, STATECD, MEASYEAR) %>%
    left_join(plot_aglb, by = c("CN" = "PLT_CN")) %>%
    filter(!is.na(aglb_Mg_per_ha))
  
  message("  Plots with AGLB: ", nrow(result))
  
  result
}

# ==============================================================================
# Build FIA Dataset (with replicate aggregation)
# ==============================================================================

#' Build FIA modeling dataset
#' 
#' Combines AGLB with extracted covariates, optionally aggregating across
#' jitter replicates to get mean and SD of covariate values.
#' 
#' @param fia_aglb FIA plot-level AGLB
#' @param fia_ndvi_path Path to FIA NDVI extractions
#' @param fia_climate_path Path to FIA climate extractions
#' @param hex_path Path to plot hex assignments
#' @param aggregate If TRUE, compute mean/SD across replicates
#' @return data.frame for modeling
build_fia_dataset <- function(fia_aglb, 
                               fia_ndvi_path = "data/processed/ndvi_at_plots/fia_ndvi.csv",
                               fia_climate_path = "data/processed/climate_at_plots/fia_climate.csv",
                               hex_path = "data/processed/plot_hex_assignments.csv",
                               aggregate = TRUE) {
  
  message("\nBuilding FIA dataset (aggregate = ", aggregate, ")...")
  
  # Start with AGLB
  fia <- fia_aglb %>%
    select(CN, STATECD, MEASYEAR, aglb_Mg_per_ha)
  
  # Load and join NDVI
  if (fs::file_exists(fia_ndvi_path)) {
    ndvi <- read_csv(fia_ndvi_path, show_col_types = FALSE)
    message("  NDVI records: ", format(nrow(ndvi), big.mark = ","))
    
    if (aggregate) {
      # Aggregate NDVI across replicates
      ndvi_agg <- ndvi %>%
        group_by(CN) %>%
        summarise(
          ndvi_modis_mean = mean(ndvi_modis, na.rm = TRUE),
          ndvi_modis_sd = sd(ndvi_modis, na.rm = TRUE),
          ndvi_s2_mean = mean(ndvi_s2, na.rm = TRUE),
          ndvi_s2_sd = sd(ndvi_s2, na.rm = TRUE),
          n_reps = n(),
          .groups = "drop"
        ) %>%
        # Rename for simpler modeling
        rename(ndvi_modis = ndvi_modis_mean, ndvi_s2 = ndvi_s2_mean)
      
      fia <- fia %>%
        left_join(ndvi_agg, by = "CN")
    } else {
      # Use first replicate only
      ndvi_first <- ndvi %>%
        filter(replicate_id == 1) %>%
        select(CN, ndvi_modis, ndvi_s2)
      
      fia <- fia %>%
        left_join(ndvi_first, by = "CN")
    }
  } else {
    message("  NDVI data not found - skipping")
    fia$ndvi_modis <- NA_real_
    fia$ndvi_s2 <- NA_real_
  }
  
  # Load and join climate
  if (fs::file_exists(fia_climate_path)) {
    climate <- read_csv(fia_climate_path, show_col_types = FALSE)
    message("  Climate records: ", format(nrow(climate), big.mark = ","))
    
    if (aggregate) {
      climate_agg <- climate %>%
        group_by(CN) %>%
        summarise(
          tmean_mean = mean(tmean, na.rm = TRUE),
          tmean_sd = sd(tmean, na.rm = TRUE),
          ppt_mean = mean(ppt, na.rm = TRUE),
          ppt_sd = sd(ppt, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        rename(tmean = tmean_mean, ppt = ppt_mean)
      
      fia <- fia %>%
        left_join(climate_agg, by = "CN")
    } else {
      climate_first <- climate %>%
        filter(replicate_id == 1) %>%
        select(CN, tmean, ppt)
      
      fia <- fia %>%
        left_join(climate_first, by = "CN")
    }
  } else {
    message("  Climate data not found - skipping")
    fia$tmean <- NA_real_
    fia$ppt <- NA_real_
  }
  
  # Load and join hex assignments
  if (fs::file_exists(hex_path)) {
    hex <- read_csv(hex_path, show_col_types = FALSE) %>%
      select(CN, starts_with("hex_id_"))
    
    fia <- fia %>%
      left_join(hex, by = "CN")
  }
  
  # Add source identifier
  fia <- fia %>%
    mutate(source = "FIA")
  
  message("  Final FIA dataset: ", nrow(fia), " plots")
  
  fia
}

# ==============================================================================
# Build NEFIN Dataset
# ==============================================================================

#' Build NEFIN modeling dataset
#' 
#' @param nefin_path Path to NEFIN processed data
#' @param nefin_ndvi_path Path to NEFIN NDVI extractions
#' @param nefin_climate_path Path to NEFIN climate extractions
#' @param nefin_hex_path Path to NEFIN hex assignments
#' @return data.frame for modeling
build_nefin_dataset <- function(nefin_path = "data/processed/nefin_processed.csv",
                                 nefin_ndvi_path = "data/processed/ndvi_at_plots/nefin_ndvi.csv",
                                 nefin_climate_path = "data/processed/climate_at_plots/nefin_climate.csv",
                                 nefin_hex_path = "data/processed/nefin_hex_assignments.csv") {
  
  message("\nBuilding NEFIN dataset...")
  
  if (!fs::file_exists(nefin_path)) {
    message("  NEFIN processed data not found")
    return(NULL)
  }
  
  # Load NEFIN base data
  nefin <- read_csv(nefin_path, show_col_types = FALSE) %>%
    select(CN, STATECD, MEASYEAR, aglb_Mg_per_ha)
  
  message("  Base NEFIN: ", nrow(nefin), " plots")
  
  # Load and join NDVI
  if (fs::file_exists(nefin_ndvi_path)) {
    ndvi <- read_csv(nefin_ndvi_path, show_col_types = FALSE) %>%
      select(CN, ndvi_modis, ndvi_s2)
    
    nefin <- nefin %>%
      left_join(ndvi, by = "CN")
  } else {
    message("  NDVI data not found - skipping")
    nefin$ndvi_modis <- NA_real_
    nefin$ndvi_s2 <- NA_real_
  }
  
  # Load and join climate
  if (fs::file_exists(nefin_climate_path)) {
    climate <- read_csv(nefin_climate_path, show_col_types = FALSE) %>%
      select(CN, tmean, ppt)
    
    nefin <- nefin %>%
      left_join(climate, by = "CN")
  } else {
    message("  Climate data not found - skipping")
    nefin$tmean <- NA_real_
    nefin$ppt <- NA_real_
  }
  
  # Load and join hex assignments
  if (fs::file_exists(nefin_hex_path)) {
    hex <- read_csv(nefin_hex_path, show_col_types = FALSE) %>%
      select(CN, starts_with("hex_id_"))
    
    nefin <- nefin %>%
      left_join(hex, by = "CN")
  }
  
  # Add source
  nefin <- nefin %>%
    mutate(source = "NEFIN")
  
  message("  Final NEFIN dataset: ", nrow(nefin), " plots")
  
  nefin
}

# ==============================================================================
# Main Entry Point
# ==============================================================================

main <- function() {
  
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  BUILD MODELING DATASETS\n")
  cat(strrep("=", 70), "\n\n")
  
  # Parse arguments
  args <- commandArgs(trailingOnly = TRUE)
  overwrite <- "--overwrite" %in% args
  
  # Output directory
  out_dir <- fs::path("data", "processed", "modeling")
  fs::dir_create(out_dir, recurse = TRUE)
  
  # Check for existing outputs
  fia_fuzzed_out <- fs::path(out_dir, "dataset_fia_fuzzed.csv")
  fia_mean_out <- fs::path(out_dir, "dataset_fia_mean.csv")
  nefin_out <- fs::path(out_dir, "dataset_nefin.csv")
  combined_out <- fs::path(out_dir, "dataset_combined.csv")
  
  if (!overwrite && fs::file_exists(combined_out)) {
    message("Datasets already exist. Use --overwrite to regenerate.")
    message("  ", combined_out)
    return(invisible(NULL))
  }
  
  # Find FIA data
  fia_root <- if (fs::dir_exists("data/interim/fia_region")) {
    "data/interim/fia_region"
  } else {
    "data/interim/fia"
  }
  
  tree_path <- fs::path(fia_root, "tree.csv")
  plot_path <- fs::path(fia_root, "plot.csv")
  
  # Compute FIA AGLB
  fia_aglb <- compute_fia_aglb(tree_path, plot_path)
  
  # Build FIA datasets
  # 1. Single replicate (fuzzed)
  message("\n--- Dataset 1: FIA Fuzzed (single replicate) ---")
  fia_fuzzed <- build_fia_dataset(fia_aglb, aggregate = FALSE)
  write_csv(fia_fuzzed, fia_fuzzed_out)
  message("Saved: ", fia_fuzzed_out)
  
  # 2. Mean across replicates
  message("\n--- Dataset 2: FIA Mean (aggregated replicates) ---")
  fia_mean <- build_fia_dataset(fia_aglb, aggregate = TRUE)
  write_csv(fia_mean, fia_mean_out)
  message("Saved: ", fia_mean_out)
  
  # Build NEFIN dataset
  message("\n--- Dataset 3: NEFIN (true coordinates) ---")
  nefin <- build_nefin_dataset()
  
  if (!is.null(nefin)) {
    write_csv(nefin, nefin_out)
    message("Saved: ", nefin_out)
  }
  
  # Build combined dataset
  message("\n--- Dataset 4: Combined (FIA + NEFIN) ---")
  
  # Use FIA mean for combination
  if (!is.null(nefin)) {
    # Ensure consistent columns
    common_cols <- intersect(names(fia_mean), names(nefin))
    
    combined <- bind_rows(
      fia_mean %>% select(all_of(common_cols)),
      nefin %>% select(all_of(common_cols))
    )
    
    message("  Combined: ", nrow(combined), " plots")
    message("    FIA: ", sum(combined$source == "FIA"))
    message("    NEFIN: ", sum(combined$source == "NEFIN"))
    
    write_csv(combined, combined_out)
    message("Saved: ", combined_out)
  } else {
    combined <- fia_mean
    write_csv(combined, combined_out)
    message("Combined is FIA-only (no NEFIN data)")
  }
  
  # Write summary
  summary <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    datasets = list(
      fia_fuzzed = list(
        path = as.character(fia_fuzzed_out),
        n_plots = nrow(fia_fuzzed),
        description = "FIA with single jitter replicate"
      ),
      fia_mean = list(
        path = as.character(fia_mean_out),
        n_plots = nrow(fia_mean),
        description = "FIA with mean across jitter replicates"
      ),
      nefin = list(
        path = as.character(nefin_out),
        n_plots = if (!is.null(nefin)) nrow(nefin) else 0,
        description = "NEFIN with true coordinates"
      ),
      combined = list(
        path = as.character(combined_out),
        n_plots = nrow(combined),
        n_fia = sum(combined$source == "FIA"),
        n_nefin = sum(combined$source == "NEFIN"),
        description = "FIA (mean) + NEFIN combined"
      )
    ),
    columns = list(
      response = "aglb_Mg_per_ha",
      ndvi = c("ndvi_modis", "ndvi_s2"),
      climate = c("tmean", "ppt"),
      identifiers = c("CN", "STATECD", "MEASYEAR", "source"),
      hex_scales = names(fia_mean)[grepl("^hex_id_", names(fia_mean))]
    )
  )
  
  summary_path <- fs::path(out_dir, "dataset_summary.yml")
  write_yaml(summary, summary_path)
  message("\nSaved summary: ", summary_path)
  
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  DATASET BUILD COMPLETE\n")
  cat(strrep("=", 70), "\n\n")
}

# Run if called directly
if (!interactive()) {
  main()
}
