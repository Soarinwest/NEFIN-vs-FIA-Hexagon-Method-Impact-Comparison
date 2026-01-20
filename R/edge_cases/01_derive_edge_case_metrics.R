#!/usr/bin/env Rscript
# =============================================================================
# 01_derive_edge_case_metrics.R
# 
# Derive plot-level metrics for large-tree structure and mortality
# from FIA and NEFIN tree-level data
#
# GOAL: Characterize distributional differences, NOT spatial accuracy
# 
# Metrics derived:
#   - Large-tree structure: max_dbh, max_tree_biomass, p95_dbh, etc.
#   - Mortality: dead_tree_count, dead_biomass_fraction, mortality_ratio
#
# Author: Soren Donisvitch
# Date: December 2024
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(fs)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

config <- list(
  # FIA tree-level data
  fia_tree_path = "data/interim/fia_region/tree.csv",
  fia_cond_path = "data/interim/fia_region/cond.csv",
  fia_plot_path = "data/processed/fia_complete.csv",
  

  # NEFIN tree-level data
  nefin_tree_path = "data/raw/nefin/TREE_PLOT_DATA.csv",
  nefin_plot_path = "data/processed/nefin_processed.csv",
  
  # Alternative paths to check
  fia_tree_alt = "data/interim/fia/states",  # May need to combine state files
  
  # DBH threshold for "large tree" (cm)
  large_dbh_threshold = 40,  # ~16 inches
  
  # Output
  output_dir = "data/processed/edge_case_metrics",
  output_file = "plot_edge_case_metrics.csv"
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Compute plot-level metrics from tree data
compute_plot_metrics <- function(tree_df, 
                                  dbh_col = "DIA",
                                  biomass_col = "DRYBIO_AG",
                                  status_col = "STATUSCD",
                                  plot_id_col = "PLT_CN",
                                  tpa_col = "TPA_UNADJ",
                                  large_dbh = 40) {
  
  # Standardize column names
  tree_df <- tree_df %>%
    rename(
      dbh = all_of(dbh_col),
      plot_id = all_of(plot_id_col)
    )
  
  # Add biomass column if exists
  if (biomass_col %in% names(tree_df)) {
    tree_df <- tree_df %>% rename(biomass = all_of(biomass_col))
  } else {
    # Estimate biomass from DBH if not available (rough approximation)
    # Using general hardwood allometry: biomass = exp(-2.5 + 2.4 * ln(dbh))
    tree_df <- tree_df %>%
      mutate(biomass = exp(-2.5 + 2.4 * log(dbh)))
    warning("Biomass column not found, using DBH-based approximation")
  }
  
  # Add status column (1 = live, 2 = dead in FIA)
  if (status_col %in% names(tree_df)) {
    tree_df <- tree_df %>% 
      rename(status = all_of(status_col)) %>%
      mutate(is_dead = status == 2)
  } else {
    tree_df <- tree_df %>% mutate(is_dead = FALSE)
    warning("Status column not found, assuming all trees live")
  }
  
  # Add TPA if exists (for proper expansion)
  if (tpa_col %in% names(tree_df)) {
    tree_df <- tree_df %>% rename(tpa = all_of(tpa_col))
  } else {
    tree_df <- tree_df %>% mutate(tpa = 1)
  }
  
  # Filter to valid trees
  tree_df <- tree_df %>%
    filter(!is.na(dbh), dbh > 0)
  
  # Compute plot-level metrics
  plot_metrics <- tree_df %>%
    group_by(plot_id) %>%
    summarize(
      # Sample size
      n_trees = n(),
      n_trees_weighted = sum(tpa, na.rm = TRUE),
      
      # ===== LARGE-TREE STRUCTURE =====
      # Maximum values (most extreme tree)
      max_dbh = max(dbh, na.rm = TRUE),
      max_tree_biomass = max(biomass, na.rm = TRUE),
      
      # 95th percentile (robust extreme)
      p95_dbh = quantile(dbh, 0.95, na.rm = TRUE),
      p95_tree_biomass = quantile(biomass, 0.95, na.rm = TRUE),
      
      # Mean and median for context
      mean_dbh = mean(dbh, na.rm = TRUE),
      median_dbh = median(dbh, na.rm = TRUE),
      
      # Large tree counts
      n_large_trees = sum(dbh >= large_dbh, na.rm = TRUE),
      pct_large_trees = 100 * sum(dbh >= large_dbh, na.rm = TRUE) / n(),
      
      # Basal area of large trees
      ba_total = sum(pi * (dbh/200)^2, na.rm = TRUE),  # m² 
      ba_large = sum(pi * (dbh/200)^2 * (dbh >= large_dbh), na.rm = TRUE),
      pct_ba_large = 100 * ba_large / ba_total,
      
      # Biomass concentration (inequality)
      total_biomass = sum(biomass, na.rm = TRUE),
      top1_biomass_pct = 100 * max(biomass, na.rm = TRUE) / sum(biomass, na.rm = TRUE),
      top3_biomass_pct = 100 * sum(sort(biomass, decreasing = TRUE)[1:min(3, n())]) / sum(biomass, na.rm = TRUE),
      
      # ===== MORTALITY / DISTURBANCE =====
      # Dead tree counts
      n_dead = sum(is_dead, na.rm = TRUE),
      n_live = sum(!is_dead, na.rm = TRUE),
      mortality_ratio = n_dead / n_trees,
      
      # Dead biomass
      dead_biomass = sum(biomass * is_dead, na.rm = TRUE),
      live_biomass = sum(biomass * (!is_dead), na.rm = TRUE),
      dead_biomass_fraction = dead_biomass / (dead_biomass + live_biomass + 1e-10),
      
      # Dead large trees (combination metric)
      n_dead_large = sum(is_dead & dbh >= large_dbh, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  plot_metrics
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

derive_edge_case_metrics <- function(cfg = config) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
  cat("║                                                                              ║\n")
  cat("║   Edge Case Metrics: Large-Tree Structure & Mortality                        ║\n")
  cat("║   FIA vs NEFIN Distributional Comparison                                     ║\n")
  cat("║                                                                              ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  fs::dir_create(cfg$output_dir, recurse = TRUE)
  
  results <- list()
  
  # ===========================================================================
  # 1. PROCESS FIA TREE DATA
  # ===========================================================================
  
  cat("═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 1: Processing FIA Tree-Level Data\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  if (fs::file_exists(cfg$fia_tree_path)) {
    cat("Loading FIA tree data:", cfg$fia_tree_path, "\n")
    fia_tree <- read_csv(cfg$fia_tree_path, show_col_types = FALSE)
    cat("  Rows:", nrow(fia_tree), "\n")
    cat("  Columns:", paste(names(fia_tree), collapse = ", "), "\n")
    
    # Identify column names (FIA uses DIA for DBH, STATUSCD for live/dead)
    dbh_col <- intersect(names(fia_tree), c("DIA", "DBH", "DIAM", "dia"))[1]
    biomass_col <- intersect(names(fia_tree), c("DRYBIO_AG", "CARBON_AG", "BIOMASS", "drybio_ag"))[1]
    status_col <- intersect(names(fia_tree), c("STATUSCD", "TREE_STATUS", "STATUS", "statuscd"))[1]
    plot_col <- intersect(names(fia_tree), c("PLT_CN", "CN", "PLOT_CN", "plt_cn"))[1]
    tpa_col <- intersect(names(fia_tree), c("TPA_UNADJ", "TPA", "tpa_unadj"))[1]
    
    cat("\n  Identified columns:\n")
    cat("    DBH:", ifelse(is.na(dbh_col), "NOT FOUND", dbh_col), "\n")
    cat("    Biomass:", ifelse(is.na(biomass_col), "NOT FOUND", biomass_col), "\n")
    cat("    Status:", ifelse(is.na(status_col), "NOT FOUND", status_col), "\n")
    cat("    Plot ID:", ifelse(is.na(plot_col), "NOT FOUND", plot_col), "\n")
    cat("    TPA:", ifelse(is.na(tpa_col), "NOT FOUND", tpa_col), "\n")
    
    if (!is.na(dbh_col) && !is.na(plot_col)) {
      cat("\n  Computing plot-level metrics...\n")
      
      fia_metrics <- compute_plot_metrics(
        fia_tree,
        dbh_col = dbh_col,
        biomass_col = ifelse(is.na(biomass_col), "NONE", biomass_col),
        status_col = ifelse(is.na(status_col), "NONE", status_col),
        plot_id_col = plot_col,
        tpa_col = ifelse(is.na(tpa_col), "NONE", tpa_col),
        large_dbh = cfg$large_dbh_threshold
      )
      
      fia_metrics <- fia_metrics %>%
        mutate(source = "FIA")
      
      cat("  ✓ FIA plots with metrics:", nrow(fia_metrics), "\n")
      results$fia <- fia_metrics
    } else {
      cat("  ✗ Missing required columns\n")
    }
    
  } else {
    cat("FIA tree file not found:", cfg$fia_tree_path, "\n")
    cat("Checking alternative paths...\n")
    
    # Try to combine state files
    state_dirs <- list.dirs(cfg$fia_tree_alt, recursive = FALSE)
    if (length(state_dirs) > 0) {
      cat("  Found state directories:", length(state_dirs), "\n")
      # Would need to combine - complex
    }
  }
  
  # ===========================================================================
  # 2. PROCESS NEFIN TREE DATA
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 2: Processing NEFIN Tree-Level Data\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  if (fs::file_exists(cfg$nefin_tree_path)) {
    cat("Loading NEFIN tree data:", cfg$nefin_tree_path, "\n")
    nefin_tree <- read_csv(cfg$nefin_tree_path, show_col_types = FALSE)
    cat("  Rows:", nrow(nefin_tree), "\n")
    cat("  Columns:", paste(head(names(nefin_tree), 20), collapse = ", "), "...\n")
    
    # Identify column names (NEFIN may use different conventions)
    dbh_col <- intersect(names(nefin_tree), c("DBH", "DIA", "dbh", "DIAM", "dia_cm"))[1]
    biomass_col <- intersect(names(nefin_tree), c("BIOMASS", "DRYBIO_AG", "AGB", "biomass", "tree_biomass"))[1]
    status_col <- intersect(names(nefin_tree), c("STATUS", "STATUSCD", "TREE_STATUS", "status", "live_dead"))[1]
    plot_col <- intersect(names(nefin_tree), c("PLOT_ID", "CN", "PLT_CN", "plot_id", "PLOT"))[1]
    
    cat("\n  Identified columns:\n")
    cat("    DBH:", ifelse(is.na(dbh_col), "NOT FOUND", dbh_col), "\n")
    cat("    Biomass:", ifelse(is.na(biomass_col), "NOT FOUND", biomass_col), "\n")
    cat("    Status:", ifelse(is.na(status_col), "NOT FOUND", status_col), "\n")
    cat("    Plot ID:", ifelse(is.na(plot_col), "NOT FOUND", plot_col), "\n")
    
    if (!is.na(dbh_col) && !is.na(plot_col)) {
      cat("\n  Computing plot-level metrics...\n")
      
      nefin_metrics <- compute_plot_metrics(
        nefin_tree,
        dbh_col = dbh_col,
        biomass_col = ifelse(is.na(biomass_col), "NONE", biomass_col),
        status_col = ifelse(is.na(status_col), "NONE", status_col),
        plot_id_col = plot_col,
        tpa_col = "NONE",
        large_dbh = cfg$large_dbh_threshold
      )
      
      nefin_metrics <- nefin_metrics %>%
        mutate(source = "NEFIN")
      
      cat("  ✓ NEFIN plots with metrics:", nrow(nefin_metrics), "\n")
      results$nefin <- nefin_metrics
    } else {
      cat("  ✗ Missing required columns\n")
    }
    
  } else {
    cat("NEFIN tree file not found:", cfg$nefin_tree_path, "\n")
    
    # Try to derive from plot-level if available
    if (fs::file_exists(cfg$nefin_plot_path)) {
      cat("Checking plot-level data for existing metrics...\n")
      nefin_plot <- read_csv(cfg$nefin_plot_path, show_col_types = FALSE)
      cat("  NEFIN plot columns:", paste(names(nefin_plot), collapse = ", "), "\n")
      
      # Check if QMD or other tree-size proxies exist
      if ("QMD" %in% names(nefin_plot)) {
        cat("  Found QMD (quadratic mean diameter) - can use as large-tree proxy\n")
      }
    }
  }
  
  # ===========================================================================
  # 3. COMBINE AND SAVE
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 3: Combining Results\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  if (length(results) > 0) {
    combined <- bind_rows(results)
    
    cat("Combined dataset:\n")
    cat("  Total plots:", nrow(combined), "\n")
    cat("  FIA plots:", sum(combined$source == "FIA"), "\n")
    cat("  NEFIN plots:", sum(combined$source == "NEFIN"), "\n")
    
    # Save
    output_path <- fs::path(cfg$output_dir, cfg$output_file)
    write_csv(combined, output_path)
    cat("\n  ✓ Saved to:", output_path, "\n")
    
    # Quick summary
    cat("\n")
    cat("═══════════════════════════════════════════════════════════════════════════════\n")
    cat("QUICK SUMMARY (Large-Tree Metrics)\n")
    cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
    
    summary_stats <- combined %>%
      group_by(source) %>%
      summarize(
        n_plots = n(),
        
        # Max DBH
        max_dbh_median = median(max_dbh, na.rm = TRUE),
        max_dbh_p95 = quantile(max_dbh, 0.95, na.rm = TRUE),
        max_dbh_max = max(max_dbh, na.rm = TRUE),
        
        # Mortality
        mortality_ratio_median = median(mortality_ratio, na.rm = TRUE),
        dead_biomass_frac_median = median(dead_biomass_fraction, na.rm = TRUE),
        
        # Large trees
        pct_with_large_trees = 100 * mean(n_large_trees > 0, na.rm = TRUE),
        
        .groups = "drop"
      )
    
    print(summary_stats)
    
    return(combined)
    
  } else {
    cat("No data processed. Check file paths.\n")
    return(NULL)
  }
}

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  derive_edge_case_metrics()
}
