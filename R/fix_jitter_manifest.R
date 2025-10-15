# R/fix_jitter_manifest.R
# Recreate manifest.yml for existing jitter library

suppressPackageStartupMessages({
  library(fs); library(yaml); library(readr); library(dplyr)
})

fix_jitter_manifest <- function(project_dir = ".") {
  
  message("\n═══════════════════════════════════════════════════════════")
  message("FIXING JITTER LIBRARY MANIFEST")
  message("═══════════════════════════════════════════════════════════\n")
  
  jitter_dir <- fs::path(project_dir, "data", "processed", "mc_jitter_library")
  replicates_dir <- fs::path(jitter_dir, "replicates")
  manifest_file <- fs::path(jitter_dir, "manifest.yml")
  
  # Check if replicates exist
  if (!fs::dir_exists(replicates_dir)) {
    stop("No replicates directory found at: ", replicates_dir)
  }
  
  # Find replicate files
  rep_files <- list.files(replicates_dir, pattern = "^rep_\\d{4}\\.csv$", full.names = TRUE)
  
  if (length(rep_files) == 0) {
    stop("No replicate CSV files found in: ", replicates_dir)
  }
  
  message("→ Found ", length(rep_files), " replicate files")
  
  # Read first replicate to get structure
  message("→ Analyzing replicate structure...")
  first_rep <- readr::read_csv(rep_files[1], show_col_types = FALSE)
  
  message("  Columns: ", paste(names(first_rep), collapse = ", "))
  message("  Rows: ", nrow(first_rep))
  
  # Detect hex grid columns
  hex_cols <- names(first_rep)[grepl("^hex_id_", names(first_rep))]
  
  if (length(hex_cols) == 0) {
    message("  ⚠ No hex_id_* columns found - this is SINGLE-SCALE jitter library")
    message("    (Old format with just 'hex_id_jittered' column)")
    
    # Single-scale manifest
    manifest <- list(
      created = as.character(Sys.time()),
      n_replicates = length(rep_files),
      n_plots = nrow(first_rep),
      radius_m = 1609.34,  # Default
      format = "csv",
      constrained = TRUE,
      used_hex_union = TRUE,
      used_state_constraint = TRUE,
      hex_grids = NULL,  # Single-scale (no multi-grid info)
      replicates_dir = "replicates",
      columns = names(first_rep)
    )
    
  } else {
    message("  ✓ Multi-scale jitter library detected")
    message("    Hex grids: ", paste(gsub("^hex_id_", "", hex_cols), collapse = ", "))
    
    # Load config to get hex grid paths
    cfg_path <- "configs/process.yml"
    hex_grids <- if (file.exists(cfg_path)) {
      cfg <- yaml::read_yaml(cfg_path)
      cfg$hex_grids
    } else {
      NULL
    }
    
    if (is.null(hex_grids)) {
      message("  ⚠ configs/process.yml not found - creating manifest without grid paths")
      
      # Create minimal grid info from column names
      grid_names <- gsub("^hex_id_", "", hex_cols)
      hex_grids <- lapply(grid_names, function(name) {
        list(
          name = name,
          path = paste0("data/hex/hex_grid_", name, ".geojson")
        )
      })
    }
    
    # Multi-scale manifest
    manifest <- list(
      created = as.character(Sys.time()),
      n_replicates = length(rep_files),
      n_plots = nrow(first_rep),
      radius_m = 1609.34,  # Default
      format = "csv",
      constrained = TRUE,
      used_hex_union = TRUE,
      used_state_constraint = TRUE,
      hex_grids = hex_grids,
      replicates_dir = "replicates",
      columns = names(first_rep)
    )
  }
  
  # Check for coordinate columns
  has_coords <- all(c("lon_jittered", "lat_jittered") %in% names(first_rep))
  if (!has_coords) {
    message("  ⚠ Missing jittered coordinates (lon_jittered, lat_jittered)")
  } else {
    message("  ✓ Jittered coordinates present")
  }
  
  # Write manifest
  message("\n→ Writing manifest...")
  yaml::write_yaml(manifest, manifest_file)
  message("  ✓ Wrote: ", manifest_file)
  
  # Summary
  message("\n═══════════════════════════════════════════════════════════")
  message("MANIFEST CREATED")
  message("═══════════════════════════════════════════════════════════")
  message("  Replicates: ", manifest$n_replicates)
  message("  Plots: ", manifest$n_plots)
  message("  Format: ", manifest$format)
  
  if (!is.null(manifest$hex_grids)) {
    message("  Hex grids: ", length(manifest$hex_grids))
    for (grid in manifest$hex_grids) {
      message("    - ", grid$name)
    }
  } else {
    message("  Type: Single-scale (legacy format)")
  }
  
  message("\n✓ Jitter library is now ready to use!")
  message("  Next step: Rscript R/06_compute_metrics.R")
  message("═══════════════════════════════════════════════════════════\n")
  
  invisible(manifest_file)
}

# Run if called directly
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  fix_jitter_manifest(".")
}