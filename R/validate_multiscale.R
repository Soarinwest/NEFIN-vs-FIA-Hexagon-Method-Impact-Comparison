# R/validate_multiscale.R
# Validate multi-scale jitter library setup

suppressPackageStartupMessages({
  library(fs); library(yaml); library(readr); library(dplyr)
})

validate_multiscale_setup <- function(project_dir = ".") {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║       Multi-Scale Setup Validation                       ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  all_ok <- TRUE
  
  # ========================================================================
  # 1. Check config file
  # ========================================================================
  cat("1. Checking configs/process.yml...\n")
  
  cfg_path <- "configs/process.yml"
  if (!file.exists(cfg_path)) {
    cat("   ✗ Missing: ", cfg_path, "\n")
    all_ok <- FALSE
  } else {
    cfg <- yaml::read_yaml(cfg_path)
    
    if (is.null(cfg$hex_grids)) {
      cat("   ✗ No hex_grids defined in config\n")
      cat("      Add hex_grids section with multiple grids\n")
      all_ok <- FALSE
    } else {
      grid_names <- sapply(cfg$hex_grids, function(x) x$name)
      cat("   ✓ Found ", length(cfg$hex_grids), " hex grids: ", 
          paste(grid_names, collapse = ", "), "\n")
      
      # Check each grid file exists
      for (grid in cfg$hex_grids) {
        if (!file.exists(grid$path)) {
          cat("   ✗ Grid file missing: ", grid$path, "\n")
          cat("      Run: Rscript R/create_hex_grid.R\n")
          all_ok <- FALSE
        } else {
          cat("   ✓ ", grid$name, ": ", grid$path, "\n")
        }
      }
    }
  }
  
  cat("\n")
  
  # ========================================================================
  # 2. Check jitter library
  # ========================================================================
  cat("2. Checking jitter library...\n")
  
  jitter_dir <- fs::path(project_dir, "data", "processed", "mc_jitter_library")
  manifest_file <- fs::path(jitter_dir, "manifest.yml")
  
  if (!fs::dir_exists(jitter_dir)) {
    cat("   ✗ Jitter library not found: ", jitter_dir, "\n")
    cat("      Run: Rscript run_pipeline.R --jitter\n")
    all_ok <- FALSE
  } else if (!fs::file_exists(manifest_file)) {
    cat("   ✗ Manifest missing: ", manifest_file, "\n")
    cat("      Run: Rscript run_pipeline.R --jitter --overwrite\n")
    all_ok <- FALSE
  } else {
    manifest <- yaml::read_yaml(manifest_file)
    
    cat("   ✓ Library created: ", manifest$created, "\n")
    cat("   ✓ Expected replicates: ", manifest$n_replicates, "\n")
    
    # Check actual replicates
    reps_dir <- fs::path(jitter_dir, "replicates")
    rep_files <- list.files(reps_dir, pattern = "^rep_\\d{4}\\.csv$")
    cat("   ✓ Available replicates: ", length(rep_files), "\n")
    
    if (length(rep_files) < manifest$n_replicates) {
      cat("   ⚠ Partial library (", length(rep_files), " of ", 
          manifest$n_replicates, ")\n")
      cat("      Resume with: Rscript run_pipeline.R --jitter\n")
    }
    
    # Check if multi-scale
    if (is.null(manifest$hex_grids)) {
      cat("   ✗ Library is SINGLE-scale (old format)\n")
      cat("      Rebuild with: Rscript run_pipeline.R --jitter --overwrite\n")
      all_ok <- FALSE
    } else {
      lib_grids <- sapply(manifest$hex_grids, function(x) x$name)
      cat("   ✓ Multi-scale library with grids: ", paste(lib_grids, collapse = ", "), "\n")
      
      # Check if config grids match library grids
      if (!is.null(cfg$hex_grids)) {
        cfg_grids <- sapply(cfg$hex_grids, function(x) x$name)
        
        if (!setequal(cfg_grids, lib_grids)) {
          cat("   ⚠ Config grids don't match library grids\n")
          cat("      Config: ", paste(cfg_grids, collapse = ", "), "\n")
          cat("      Library: ", paste(lib_grids, collapse = ", "), "\n")
          cat("      Rebuild library to match config\n")
          all_ok <- FALSE
        } else {
          cat("   ✓ Config and library grids match\n")
        }
      }
      
      # Check first replicate for column structure
      if (length(rep_files) > 0) {
        first_rep <- fs::path(reps_dir, rep_files[1])
        rep_data <- suppressMessages(readr::read_csv(first_rep, n_max = 5, show_col_types = FALSE))
        
        cat("\n   Checking replicate file structure...\n")
        
        expected_cols <- paste0("hex_id_", lib_grids)
        found_cols <- intersect(expected_cols, names(rep_data))
        
        if (length(found_cols) != length(expected_cols)) {
          cat("   ✗ Missing hex_id columns in replicates\n")
          cat("      Expected: ", paste(expected_cols, collapse = ", "), "\n")
          cat("      Found: ", paste(found_cols, collapse = ", "), "\n")
          all_ok <- FALSE
        } else {
          cat("   ✓ All hex_id columns present: ", paste(found_cols, collapse = ", "), "\n")
          
          # Show example assignments
          cat("\n   Example assignments (first 3 plots):\n")
          example <- rep_data[1:min(3, nrow(rep_data)), c("CN", found_cols)]
          print(as.data.frame(example), row.names = FALSE)
        }
        
        # Check for coordinate columns
        if (all(c("lon_jittered", "lat_jittered") %in% names(rep_data))) {
          cat("   ✓ Jittered coordinates present\n")
        } else {
          cat("   ✗ Missing jittered coordinates\n")
          all_ok <- FALSE
        }
      }
    }
  }
  
  cat("\n")
  
  # ========================================================================
  # 3. Check plot assignments
  # ========================================================================
  cat("3. Checking plot assignments...\n")
  
  assignments_file <- fs::path(project_dir, "data", "processed", "plot_hex_assignments.csv")
  
  if (!fs::file_exists(assignments_file)) {
    cat("   ✗ Missing: ", assignments_file, "\n")
    cat("      Run: Rscript run_pipeline.R --assign\n")
    all_ok <- FALSE
  } else {
    assignments <- suppressMessages(readr::read_csv(assignments_file, show_col_types = FALSE))
    cat("   ✓ Plot assignments exist\n")
    cat("   ✓ Total plots: ", nrow(assignments), "\n")
    cat("   ✓ Plots with hex: ", sum(!is.na(assignments$hex_id)), "\n")
    
    if (!all(c("lat_original", "lon_original") %in% names(assignments))) {
      cat("   ✗ Missing original coordinates\n")
      cat("      Re-run: Rscript run_pipeline.R --assign --overwrite\n")
      all_ok <- FALSE
    } else {
      cat("   ✓ Original coordinates preserved\n")
    }
  }
  
  cat("\n")
  
  # ========================================================================
  # 4. Summary
  # ========================================================================
  cat("═══════════════════════════════════════════════════════════\n")
  
  if (all_ok) {
    cat("✓ VALIDATION PASSED\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("\nReady to process all scales:\n")
    cat("  Rscript R/08_process_all_scales.R\n\n")
    cat("Or process individual grids:\n")
    cat("  Rscript run_pipeline.R --compute --analyze --grid=1.5k\n")
    cat("  Rscript run_pipeline.R --compute --analyze --grid=3k\n")
    cat("  Rscript run_pipeline.R --compute --analyze --grid=6k\n\n")
  } else {
    cat("✗ VALIDATION FAILED\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("\nFix issues above, then re-run validation:\n")
    cat("  Rscript R/validate_multiscale.R\n\n")
  }
  
  invisible(all_ok)
}

# Run validation
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  result <- validate_multiscale_setup(".")
  if (!result) quit(status = 1)
}