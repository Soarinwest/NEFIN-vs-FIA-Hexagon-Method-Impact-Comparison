# R/fix_jitter_na_assignments.R
# Fix NA hex assignments by snapping to nearest hex while preserving jittered coordinates

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(fs); library(yaml)
  library(parallel); library(pbapply)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

fix_jitter_na_assignments <- function(project_dir = ".",
                                      fix_original = TRUE,
                                      fix_jitter_library = TRUE,
                                      overwrite = FALSE) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Fix NA Hex Assignments (Snap to Nearest)                ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Load config
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  hex_grids <- cfg$hex_grids
  if (is.null(hex_grids)) {
    stop("No hex_grids defined in configs/process.yml")
  }
  
  # Turn off spherical geometry
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  crs_ref <- sf::st_crs(5070)
  
  # Load all hex grids once
  cat("→ Loading hex grids...\n")
  hex_grids_sf <- list()
  
  for (grid in hex_grids) {
    cat("  Loading ", grid$name, "...\n")
    hx <- if (is.null(grid$layer)) {
      sf::st_read(grid$path, quiet = TRUE)
    } else {
      sf::st_read(grid$path, layer = grid$layer, quiet = TRUE)
    }
    
    if (!("hex_id" %in% names(hx))) {
      if ("ID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = ID)
      else if ("OBJECTID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = OBJECTID)
      else hx$hex_id <- seq_len(nrow(hx))
    }
    
    hx$hex_id <- as.character(hx$hex_id)
    hx <- sf::st_make_valid(sf::st_zm(hx, drop = TRUE, what = "ZM"))
    hx_5070 <- sf::st_transform(hx, crs_ref)
    sf::st_crs(hx_5070) <- crs_ref
    
    hex_grids_sf[[grid$name]] <- hx_5070
  }
  
  # Helper function to snap points to nearest hex
  snap_to_nearest_hex <- function(points_sf, hex_sf, hex_col_name) {
    # Find rows with NA assignments
    na_idx <- which(is.na(points_sf[[hex_col_name]]))
    
    if (length(na_idx) == 0) {
      return(points_sf)  # No NAs to fix
    }
    
    cat("    Found ", length(na_idx), " NAs to fix...\n")
    
    # Get just the points that need fixing
    na_points <- points_sf[na_idx, ]
    
    # Find nearest hex for each NA point
    nearest_idx <- sf::st_nearest_feature(na_points, hex_sf)
    
    # Assign the hex_id from nearest hex
    points_sf[[hex_col_name]][na_idx] <- hex_sf$hex_id[nearest_idx]
    
    cat("    ✓ Fixed ", length(na_idx), " NAs\n")
    
    return(points_sf)
  }
  
  # Fix original assignments
  if (fix_original) {
    cat("\n→ Fixing original plot assignments...\n")
    
    assignments_file <- fs::path(project_dir, "data", "processed", "plot_hex_assignments.csv")
    
    if (!fs::file_exists(assignments_file)) {
      cat("  ⚠ Original assignments not found: ", assignments_file, "\n")
    } else {
      # Backup original
      if (!overwrite && !fs::file_exists(paste0(assignments_file, ".backup"))) {
        fs::file_copy(assignments_file, paste0(assignments_file, ".backup"))
        cat("  Created backup: ", paste0(assignments_file, ".backup"), "\n")
      }
      
      assignments <- readr::read_csv(assignments_file, show_col_types = FALSE)
      
      # Convert to spatial
      valid_coords <- !is.na(assignments$lat_original) & !is.na(assignments$lon_original) &
        is.finite(assignments$lat_original) & is.finite(assignments$lon_original)
      
      if (any(!valid_coords)) {
        cat("  ⚠ Skipping ", sum(!valid_coords), " rows with invalid coordinates\n")
      }
      
      assignments_valid <- assignments[valid_coords, ]
      
      pts <- sf::st_as_sf(assignments_valid, 
                          coords = c("lon_original", "lat_original"), 
                          crs = 4326, remove = FALSE)
      pts_5070 <- sf::st_transform(pts, crs_ref)
      sf::st_crs(pts_5070) <- crs_ref
      
      # Fix NAs for each grid
      for (grid_name in names(hex_grids_sf)) {
        hex_col <- paste0("hex_id_", grid_name)
        
        if (hex_col %in% names(pts_5070)) {
          cat("  Fixing ", grid_name, " assignments...\n")
          n_na_before <- sum(is.na(pts_5070[[hex_col]]))
          
          if (n_na_before > 0) {
            pts_5070 <- snap_to_nearest_hex(pts_5070, hex_grids_sf[[grid_name]], hex_col)
          } else {
            cat("    No NAs found\n")
          }
        }
      }
      
      # Convert back to dataframe
      assignments_fixed <- sf::st_drop_geometry(pts_5070)
      
      # Add back any rows with invalid coordinates (keep them as-is)
      if (any(!valid_coords)) {
        assignments_fixed <- dplyr::bind_rows(
          assignments_fixed,
          assignments[!valid_coords, ]
        )
      }
      
      # Write fixed file
      readr::write_csv(assignments_fixed, assignments_file)
      cat("  ✓ Updated: ", assignments_file, "\n")
    }
  }
  
  # Fix jitter library
  if (fix_jitter_library) {
    cat("\n→ Fixing jitter library assignments...\n")
    
    jitter_dir <- fs::path(project_dir, "data", "processed", "mc_jitter_library")
    replicates_dir <- fs::path(jitter_dir, "replicates")
    
    if (!fs::dir_exists(replicates_dir)) {
      cat("  ⚠ Jitter library not found: ", replicates_dir, "\n")
    } else {
      rep_files <- list.files(replicates_dir, pattern = "^rep_\\d{4}\\.csv$", full.names = TRUE)
      
      cat("  Found ", length(rep_files), " replicate files to fix\n")
      
      # Process each replicate file
      pb <- pbapply::pblapply(rep_files, function(rep_file) {
        # Read replicate
        rep_data <- readr::read_csv(rep_file, show_col_types = FALSE)
        
        # Convert to spatial using JITTERED coordinates
        valid_coords <- !is.na(rep_data$lat_jittered) & !is.na(rep_data$lon_jittered) &
          is.finite(rep_data$lat_jittered) & is.finite(rep_data$lon_jittered)
        
        if (any(!valid_coords)) {
          warning("Replicate ", basename(rep_file), " has ", sum(!valid_coords), 
                  " rows with invalid jittered coordinates")
        }
        
        rep_valid <- rep_data[valid_coords, ]
        
        pts <- sf::st_as_sf(rep_valid, 
                            coords = c("lon_jittered", "lat_jittered"), 
                            crs = 4326, remove = FALSE)
        pts_5070 <- sf::st_transform(pts, crs_ref)
        sf::st_crs(pts_5070) <- crs_ref
        
        # Fix NAs for each grid
        for (grid_name in names(hex_grids_sf)) {
          hex_col <- paste0("hex_id_", grid_name)
          
          if (hex_col %in% names(pts_5070)) {
            n_na_before <- sum(is.na(pts_5070[[hex_col]]))
            
            if (n_na_before > 0) {
              pts_5070 <- snap_to_nearest_hex(pts_5070, hex_grids_sf[[grid_name]], hex_col)
            }
          }
        }
        
        # Convert back to dataframe
        rep_fixed <- sf::st_drop_geometry(pts_5070)
        
        # Add back any rows with invalid coordinates
        if (any(!valid_coords)) {
          rep_fixed <- dplyr::bind_rows(
            rep_fixed,
            rep_data[!valid_coords, ]
          )
        }
        
        # Overwrite the file
        readr::write_csv(rep_fixed, rep_file)
        
        return(basename(rep_file))
      })
      
      cat("\n  ✓ Fixed all ", length(rep_files), " replicate files\n")
    }
  }
  
  # Summary report
  cat("\n═══════════════════════════════════════════════════════════")
  cat("\nSUMMARY")
  cat("\n═══════════════════════════════════════════════════════════\n")
  
  if (fix_original) {
    cat("✓ Fixed original plot assignments\n")
  }
  
  if (fix_jitter_library) {
    cat("✓ Fixed jitter library assignments\n")
  }
  
  cat("\nAll NA hex assignments have been fixed by snapping to nearest hex\n")
  cat("while preserving the jittered coordinates.\n")
  
  invisible(TRUE)
}

# Quick check function to report NA statistics
check_na_stats <- function(project_dir = ".") {
  cat("\n")
  cat("══════════════════════════════════════════════════════════\n")
  cat("NA Assignment Statistics\n")
  cat("══════════════════════════════════════════════════════════\n")
  
  # Check original assignments
  assignments_file <- fs::path(project_dir, "data", "processed", "plot_hex_assignments.csv")
  
  if (fs::file_exists(assignments_file)) {
    cat("\n→ Original assignments:\n")
    assignments <- readr::read_csv(assignments_file, show_col_types = FALSE)
    
    hex_cols <- names(assignments)[grepl("^hex_id_", names(assignments))]
    
    for (col in hex_cols) {
      n_na <- sum(is.na(assignments[[col]]))
      pct_na <- 100 * n_na / nrow(assignments)
      cat("  ", col, ": ", n_na, " NAs (", round(pct_na, 2), "%)\n", sep = "")
    }
  }
  
  # Check one replicate as example
  jitter_dir <- fs::path(project_dir, "data", "processed", "mc_jitter_library", "replicates")
  
  if (fs::dir_exists(jitter_dir)) {
    rep_files <- list.files(jitter_dir, pattern = "^rep_0001\\.csv$", full.names = TRUE)
    
    if (length(rep_files) > 0) {
      cat("\n→ Example replicate (rep_0001.csv):\n")
      rep_data <- readr::read_csv(rep_files[1], show_col_types = FALSE)
      
      hex_cols <- names(rep_data)[grepl("^hex_id_", names(rep_data))]
      
      for (col in hex_cols) {
        n_na <- sum(is.na(rep_data[[col]]))
        pct_na <- 100 * n_na / nrow(rep_data)
        cat("  ", col, ": ", n_na, " NAs (", round(pct_na, 2), "%)\n", sep = "")
      }
    }
  }
  
  invisible(TRUE)
}

# CLI
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if ("--check" %in% args) {
    check_na_stats(".")
  } else {
    fix_jitter_na_assignments(
      project_dir = ".",
      fix_original = !("--skip-original" %in% args),
      fix_jitter_library = !("--skip-jitter" %in% args),
      overwrite = "--overwrite" %in% args
    )
  }
}