#!/usr/bin/env Rscript
# =============================================================================
# run_pipeline.R - Consolidated Master Pipeline Controller
# =============================================================================
# Author: Soren Donisvitch
# Date: December 2024
# 
# This is the reorganized pipeline with scripts grouped by function:
#   00_utils/        - Shared utility functions
#   01_data_prep/    - Data acquisition and preparation
#   02_uncertainty/  - Monte Carlo positional uncertainty
#   03_comparison/   - FIA vs FIA+NEFIN comparison
#   04_analysis/     - Statistical analysis
#   05_visualization/- Publication figures
#
# Usage:
#   Rscript run_pipeline.R --all              # Full pipeline
#   Rscript run_pipeline.R --data             # Data prep only
#   Rscript run_pipeline.R --uncertainty      # MC uncertainty only
#   Rscript run_pipeline.R --compare          # Comparison only
#   Rscript run_pipeline.R --analyze          # Analysis only
#   Rscript run_pipeline.R --visualize        # Figures only
#   Rscript run_pipeline.R --help             # Show this help
# =============================================================================

suppressPackageStartupMessages({
  library(fs)
  library(yaml)
})

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
SCRIPT_DIR <- "R"
CONFIG_FILE <- "process.yml"

# Load config if exists
if (file.exists(CONFIG_FILE)) {
  config <- yaml::read_yaml(CONFIG_FILE)
} else {
  config <- list()
  warning("No process.yml found - using defaults")
}

# -----------------------------------------------------------------------------
# Argument Parsing
# -----------------------------------------------------------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  list(
    help = "--help" %in% args || "-h" %in% args,
    all = "--all" %in% args,
    data = "--data" %in% args,
    uncertainty = "--uncertainty" %in% args || "--mc" %in% args,
    compare = "--compare" %in% args,
    analyze = "--analyze" %in% args,
    visualize = "--visualize" %in% args || "--viz" %in% args,
    phase2 = "--phase2" %in% args,
    validate = "--validate" %in% args,
    force = "--force" %in% args,
    verbose = "--verbose" %in% args || "-v" %in% args
  )
}

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------
run_script <- function(path, description = NULL) {
  if (!file.exists(path)) {
    cat("  ⚠ Script not found:", path, "\n")
    return(FALSE)
  }
  
  desc <- if (!is.null(description)) description else basename(path)
  cat("  → Running:", desc, "\n")
  
  tryCatch({
    source(path, local = new.env())
    cat("  ✓ Complete:", desc, "\n")
    TRUE
  }, error = function(e) {
    cat("  ✗ Failed:", desc, "\n")
    cat("    Error:", conditionMessage(e), "\n")
    FALSE
  })
}

print_header <- function(stage_name) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n")
  cat(" ", stage_name, "\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n")
}

print_help <- function() {
  cat("
FIA-NEFIN Pipeline Controller
=============================

Usage: Rscript run_pipeline.R [OPTIONS]

OPTIONS:
  --all           Run complete Phase 1 pipeline (stages 1-5)
  --data          Stage 1: Data preparation only
  --uncertainty   Stage 2: Monte Carlo uncertainty analysis
  --compare       Stage 3: FIA vs FIA+NEFIN comparison
  --analyze       Stage 4: Statistical analysis
  --visualize     Stage 5: Generate publication figures
  --phase2        Stage 6: Phase 2 spatial modeling (run separately)
  --validate      Run data validation checks
  --force         Force re-run even if outputs exist
  --verbose, -v   Show detailed output
  --help, -h      Show this help message

EXAMPLES:
  # Full Phase 1 pipeline from scratch
  Rscript run_pipeline.R --all

  # Just rerun comparison and analysis (after fixing a bug)
  Rscript run_pipeline.R --compare --analyze --visualize

  # Run Phase 2 spatial modeling
  Rscript run_pipeline.R --phase2

  # Validate data before running
  Rscript run_pipeline.R --validate

PIPELINE STAGES:
  1. DATA PREP      (01_data_prep/)
     - Initialize project structure
     - Create/filter hex grids
     - Pull FIA data
     - Process NEFIN data
     - Assign plots to hexes
     - Extract covariates

  2. UNCERTAINTY    (02_uncertainty/)
     - Build Monte Carlo jitter library
     - Compute metrics with positional uncertainty

  3. COMPARISON     (03_comparison/)
     - Compare FIA-only vs FIA+NEFIN estimates
     - Process all scales
     - Consolidate results

  4. ANALYSIS       (04_analysis/)
     - Error decomposition
     - Advanced statistics
     - Phase 1 hypothesis tests
     - NEFIN dominance/bias analysis

  5. VISUALIZATION  (05_visualization/)
     - NEFIN impact figures
     - Consolidated results dashboard
     - Spatial visualizations

  6. PHASE 2        (_archive/phase2/) - Run separately with --phase2
     - Spatial model comparison (XGBoost)
     - Fuzzing effect analysis
     - Sensor resolution comparison
     - Landscape heterogeneity analysis
     - Prediction maps
")
}

# -----------------------------------------------------------------------------
# Stage Functions
# -----------------------------------------------------------------------------

run_data_prep <- function(verbose = FALSE) {
  print_header("STAGE 1: DATA PREPARATION")
  
  scripts <- c(
    "R/01_data_prep/01_init_project.R",
    "R/01_data_prep/02_create_hex_grid.R",
    "R/01_data_prep/03_filter_hex_grid.R",
    "R/01_data_prep/04_fia_pull.R",
    "R/01_data_prep/05_process_nefin.R",
    "R/01_data_prep/06_assign_plots.R",
    "R/01_data_prep/07_assign_nefin.R",
    "R/01_data_prep/08_extract_covariates.R"
  )
  
  success <- TRUE
  for (script in scripts) {
    if (!run_script(script)) success <- FALSE
  }
  
  return(success)
}

run_uncertainty <- function(verbose = FALSE) {
  print_header("STAGE 2: MONTE CARLO UNCERTAINTY")
  
  scripts <- c(
    "R/02_uncertainty/01_build_jitter_library.R",
    "R/02_uncertainty/02_compute_metrics.R"
  )
  
  success <- TRUE
  for (script in scripts) {
    if (!run_script(script)) success <- FALSE
  }
  
  return(success)
}

run_comparison <- function(verbose = FALSE) {
  print_header("STAGE 3: FIA vs FIA+NEFIN COMPARISON")
  
  scripts <- c(
    "R/03_comparison/03_process_all_scales.R",
    "R/03_comparison/01_compare_fia_nefin.R",
    "R/03_comparison/02_consolidate_results.R"
  )
  
  success <- TRUE
  for (script in scripts) {
    if (!run_script(script)) success <- FALSE
  }
  
  return(success)
}

run_analysis <- function(verbose = FALSE) {
  print_header("STAGE 4: STATISTICAL ANALYSIS")
  
  scripts <- c(
    "R/04_analysis/01_error_analysis.R",
    "R/04_analysis/02_advanced_analysis.R",
    "R/04_analysis/03_phase1_hypothesis_tests.R",
    "R/04_analysis/04_nefin_dominance_analysis.R"
  )
  
  success <- TRUE
  for (script in scripts) {
    if (!run_script(script)) success <- FALSE
  }
  
  return(success)
}

run_visualization <- function(verbose = FALSE) {
  print_header("STAGE 5: VISUALIZATION")
  
  scripts <- c(
    "R/05_visualization/01_nefin_impact_figures.R",
    "R/05_visualization/02_visualize_results.R",
    "R/05_visualization/03_spatial_visualizations.R"
  )
  
  success <- TRUE
  for (script in scripts) {
    if (!run_script(script)) success <- FALSE
  }
  
  return(success)
}

run_validation <- function(verbose = FALSE) {
  print_header("DATA VALIDATION")
  run_script("R/validate_dashboard_data.R")
}

run_phase2 <- function(verbose = FALSE) {
  print_header("STAGE 6: PHASE 2 - SPATIAL MODELING")
  
  cat("  Phase 2 tests whether true coordinates improve spatial predictions.\n\n")
  
  # Find Phase 2 master script
  phase2_paths <- c(
    "run_phase2_complete_v2.R",
    "R/run_phase2_complete.R",
    "R/_archive/phase2/run_phase2_complete.R"
  )
  
  for (p in phase2_paths) {
    if (file.exists(p)) {
      cat("  Using:", p, "\n\n")
      return(run_script(p, "Phase 2 Complete Analysis"))
    }
  }
  
  # If no master script, run individual scripts
  cat("  Running individual Phase 2 scripts...\n\n")
  
  scripts <- c(
    "R/_archive/phase2/spatial_model_comparison_v2.R",
    "R/_archive/phase2/fuzzing_effect_analysis.R",
    "R/_archive/phase2/sensor_resolution_comparison.R",
    "R/_archive/phase2/landscape_heterogeneity_analysis.R"
  )
  
  success <- TRUE
  for (script in scripts) {
    if (file.exists(script)) {
      if (!run_script(script)) success <- FALSE
    } else {
      cat("  ⚠ Not found:", script, "\n")
    }
  }
  
  return(success)
}

# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------

main <- function() {
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  FIA-NEFIN Hexagon Analysis Pipeline (Consolidated)                      ║\n")
  cat("║  Version: 4.1                                                            ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
  
  # Source utilities
  cat("\nLoading utilities...\n")
  source("R/00_utils/utils_spatial.R")
  source("R/00_utils/utils_metrics.R")
  source("R/00_utils/utils_scale_names.R")
  cat("✓ Utilities loaded\n")
  
  args <- parse_args()
  
  if (args$help) {
    print_help()
    return(invisible(NULL))
  }
  
  if (args$validate) {
    run_validation(args$verbose)
    return(invisible(NULL))
  }
  
  # Track timing
  start_time <- Sys.time()
  
  # Determine which stages to run
  run_all <- args$all || (!args$data && !args$uncertainty && !args$compare && 
                          !args$analyze && !args$visualize && !args$phase2)
  
  # Execute stages
  if (run_all || args$data) {
    run_data_prep(args$verbose)
  }
  
  if (run_all || args$uncertainty) {
    run_uncertainty(args$verbose)
  }
  
  if (run_all || args$compare) {
    run_comparison(args$verbose)
  }
  
  if (run_all || args$analyze) {
    run_analysis(args$verbose)
  }
  
  if (run_all || args$visualize) {
    run_visualization(args$verbose)
  }
  
  # Phase 2 is NOT included in --all (run separately with --phase2)
  if (args$phase2) {
    run_phase2(args$verbose)
  }
  
  # Summary
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n")
  cat("PIPELINE COMPLETE\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n")
  cat(sprintf("Total time: %.1f minutes\n", as.numeric(elapsed)))
  cat("Outputs: runs/consolidated_*/\n")
  if (args$phase2) {
    cat("Phase 2:  runs/phase2_complete/\n")
  }
  cat("\n")
}

# Run if executed directly
if (!interactive()) {
  main()
}
