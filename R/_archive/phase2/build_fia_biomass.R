#!/usr/bin/env Rscript
# R/build_fia_biomass.R
# Compute FIA plot-level biomass using Chojnacky et al. (2014) allometric equations
# Then create complete FIA dataset with biomass + covariates
# Author: Soren Donisvitch
# Date: December 2024

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(fs)
})

# =============================================================================
# LOAD ALLOMETRIC MODEL
# =============================================================================

# Source the biomass allometric model function
# From: Chojnacky, D.C., Heath, L.S., Jenkins, J.C., 2014
source("R/Biomass_Allometric_Models.R")

# Vectorized version for efficiency
AG_Biomass_vec <- Vectorize(AG_Biomass)

# =============================================================================
# MAIN FUNCTION
# =============================================================================

build_fia_biomass <- function(
    tree_dir = "data/interim/fia/states",
    fia_cov_path = "data/processed/fia_covariates.csv",
    plot_hex_path = "data/processed/plot_hex_assignments.csv",
    output_path = "data/processed/fia_complete.csv"
) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  Build FIA Biomass Dataset                                                ║\n")
  cat("║  Using Chojnacky et al. (2014) Allometric Equations                       ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # ===========================================================================
  # 1. LOAD TREE DATA
  # ===========================================================================
  
  cat("═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 1: Loading FIA Tree Data\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  tree_files <- list.files(tree_dir, pattern = "tree\\.csv$", 
                           recursive = TRUE, full.names = TRUE)
  cat("  Tree files found:", length(tree_files), "\n")
  
  if (length(tree_files) == 0) {
    stop("No tree files found in ", tree_dir)
  }
  
  # Load all tree data
  cat("  Loading tree records...\n")
  trees <- lapply(tree_files, function(f) {
    state <- basename(dirname(f))
    df <- read_csv(f, show_col_types = FALSE)
    df$state_file <- state
    return(df)
  }) %>% bind_rows()
  
  cat("  Total tree records:", format(nrow(trees), big.mark = ","), "\n")
  
  # Check required columns
  required_cols <- c("PLT_CN", "SPCD", "DIA")
  missing_cols <- setdiff(required_cols, names(trees))
  if (length(missing_cols) > 0) {
    cat("  ⚠ Missing columns:", paste(missing_cols, collapse = ", "), "\n")
    cat("  Available columns:", paste(head(names(trees), 20), collapse = ", "), "...\n")
  }
  
  # ===========================================================================
  # 2. COMPUTE TREE-LEVEL BIOMASS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 2: Computing Tree-Level Biomass (Chojnacky et al. 2014)\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Filter to live trees with valid measurements
  # STATUSCD = 1 is live tree
  trees_live <- trees %>%
    filter(
      !is.na(DIA), DIA > 0,
      !is.na(SPCD)
    )
  
  if ("STATUSCD" %in% names(trees_live)) {
    trees_live <- trees_live %>% filter(STATUSCD == 1)
  }
  
  cat("  Live trees with valid DBH:", format(nrow(trees_live), big.mark = ","), "\n")
  
  # Convert DIA from inches to cm (FIA stores in inches)
  trees_live <- trees_live %>%
    mutate(DIA_cm = DIA * 2.54)
  
  # Compute biomass using allometric equations
  cat("  Applying allometric equations by species...\n")
  
  trees_live <- trees_live %>%
    mutate(
      biomass_kg = AG_Biomass_vec(SPCD, DIA_cm)
    )
  
  # Check results
  valid_biomass <- sum(!is.na(trees_live$biomass_kg))
  cat("  Trees with valid biomass:", format(valid_biomass, big.mark = ","), 
      "(", round(100 * valid_biomass / nrow(trees_live), 1), "%)\n")
  
  cat("  Biomass range:", round(min(trees_live$biomass_kg, na.rm = TRUE), 1), "-",
      round(max(trees_live$biomass_kg, na.rm = TRUE), 1), "kg\n")
  
  # ===========================================================================
  # 3. AGGREGATE TO PLOT LEVEL
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 3: Aggregating to Plot Level\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Get TPA adjustment factor
  # TPA_UNADJ is trees per acre represented by this tree
  if ("TPA_UNADJ" %in% names(trees_live)) {
    trees_live <- trees_live %>%
      mutate(TPA = TPA_UNADJ)
  } else {
    # Default expansion factor if not available
    cat("  ⚠ TPA_UNADJ not found, using default expansion\n")
    trees_live <- trees_live %>%
      mutate(TPA = 6.018046)  # Standard FIA subplot expansion
  }
  
  # Compute plot-level biomass
  # biomass_kg * TPA = kg/acre
  # kg/acre * 0.001 * 2.471 = Mg/ha
  plot_biomass <- trees_live %>%
    filter(!is.na(biomass_kg), !is.na(TPA)) %>%
    group_by(PLT_CN) %>%
    summarise(
      aglb_Mg_per_ha = sum(biomass_kg * TPA * 0.001 * 2.471, na.rm = TRUE),
      n_trees = n(),
      n_species = n_distinct(SPCD),
      mean_dbh_cm = mean(DIA_cm, na.rm = TRUE),
      ba_m2_ha = sum((pi * (DIA_cm/200)^2) * TPA * 2.471, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(CN = PLT_CN)
  
  cat("  Plots with computed biomass:", format(nrow(plot_biomass), big.mark = ","), "\n")
  cat("  Biomass range:", round(min(plot_biomass$aglb_Mg_per_ha), 1), "-",
      round(max(plot_biomass$aglb_Mg_per_ha), 1), "Mg/ha\n")
  cat("  Biomass mean:", round(mean(plot_biomass$aglb_Mg_per_ha), 1), "Mg/ha\n")
  cat("  Biomass median:", round(median(plot_biomass$aglb_Mg_per_ha), 1), "Mg/ha\n")
  
  # ===========================================================================
  # 4. JOIN WITH COVARIATES AND COORDINATES
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 4: Joining with Covariates and Coordinates\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Load FIA covariates
  fia_cov <- read_csv(fia_cov_path, show_col_types = FALSE)
  cat("  FIA covariates:", nrow(fia_cov), "rows\n")
  
  # Remove existing biomass column if present
  if ("aglb_Mg_per_ha" %in% names(fia_cov)) {
    fia_cov <- fia_cov %>% select(-aglb_Mg_per_ha)
  }
  
  # Load plot coordinates
  if (fs::file_exists(plot_hex_path)) {
    plot_hex <- read_csv(plot_hex_path, show_col_types = FALSE)
    cat("  Plot hex assignments:", nrow(plot_hex), "rows\n")
    
    # Get coordinate columns
    coord_cols <- c("CN", "lon_original", "lat_original")
    if (all(coord_cols %in% names(plot_hex))) {
      coords <- plot_hex %>% 
        select(all_of(coord_cols)) %>%
        distinct(CN, .keep_all = TRUE)
    } else {
      coords <- NULL
    }
  } else {
    coords <- NULL
  }
  
  # Join biomass to covariates
  fia_complete <- fia_cov %>%
    left_join(plot_biomass, by = "CN")
  
  # Add coordinates if not present
  if (!is.null(coords) && !"lon_original" %in% names(fia_complete)) {
    fia_complete <- fia_complete %>%
      left_join(coords, by = "CN")
  }
  
  cat("\n  Join results:\n")
  cat("    Total plots:", nrow(fia_complete), "\n")
  cat("    With biomass:", sum(!is.na(fia_complete$aglb_Mg_per_ha)), "\n")
  cat("    With NDVI (modis):", sum(!is.na(fia_complete$ndvi_modis)), "\n")
  cat("    With NDVI (S2):", sum(!is.na(fia_complete$ndvi_s2)), "\n")
  cat("    With climate:", sum(!is.na(fia_complete$tmean)), "\n")
  
  if ("lon_original" %in% names(fia_complete)) {
    cat("    With coordinates:", sum(!is.na(fia_complete$lon_original)), "\n")
  }
  
  # ===========================================================================
  # 5. FILTER AND SUMMARIZE
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 5: Summary by Time Window\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Complete cases
  fia_modeling <- fia_complete %>%
    filter(
      !is.na(aglb_Mg_per_ha),
      !is.na(ndvi_modis),
      !is.na(tmean),
      !is.na(ppt)
    )
  
  cat("  Plots ready for modeling (biomass + all covariates):", nrow(fia_modeling), "\n\n")
  
  # By year
  cat("  By MEASYEAR:\n")
  year_summary <- fia_modeling %>%
    count(MEASYEAR) %>%
    arrange(MEASYEAR)
  print(year_summary, n = 20)
  
  # 2020-2024 window
  fia_2020 <- fia_modeling %>% filter(MEASYEAR >= 2020, MEASYEAR <= 2024)
  cat("\n  Plots in 2020-2024 window:", nrow(fia_2020), "\n")
  
  # ===========================================================================
  # 6. SAVE
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 6: Saving Complete FIA Dataset\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  write_csv(fia_complete, output_path)
  cat("  ✓ Saved:", output_path, "\n")
  
  # Also save the plot biomass separately for reference
  biomass_path <- gsub("\\.csv$", "_biomass_only.csv", output_path)
  write_csv(plot_biomass, biomass_path)
  cat("  ✓ Saved:", biomass_path, "\n")
  
  # ===========================================================================
  # SUMMARY
  # ===========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  FIA BIOMASS DATASET COMPLETE                                             ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  cat("Summary:\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  cat("  Total FIA plots:", nrow(fia_complete), "\n")
  cat("  Plots with biomass:", sum(!is.na(fia_complete$aglb_Mg_per_ha)), "\n")
  cat("  Plots ready for modeling:", nrow(fia_modeling), "\n")
  cat("  Plots in 2020-2024:", nrow(fia_2020), "\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  cat("  Mean biomass:", round(mean(fia_complete$aglb_Mg_per_ha, na.rm = TRUE), 1), "Mg/ha\n")
  cat("  Median biomass:", round(median(fia_complete$aglb_Mg_per_ha, na.rm = TRUE), 1), "Mg/ha\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  
  cat("\nOutput files:\n")
  cat("  •", output_path, "\n")
  cat("  •", biomass_path, "\n")
  
  invisible(fia_complete)
}

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  build_fia_biomass()
}
