# R/05_build_jitter_library.R (MULTI-SCALE)
# Pre-generate N jittered coordinate sets with MULTIPLE hex grid assignments

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(fs); library(yaml); library(glue)
})

source("R/utils_spatial.R")

`%||%` <- function(a,b) if (!is.null(a)) a else b

stage4_build_jitter_library <- function(project_dir = ".",
                                        hex_grids = NULL,  # List of hex grid specs
                                        n_replicates = 100,
                                        radius_m = 1609.34,
                                        use_constraints = TRUE,
                                        overwrite = FALSE) {
  
  assignments_file <- fs::path(project_dir, "data", "processed", "plot_hex_assignments.csv")
  if (!fs::file_exists(assignments_file)) {
    stop("Run R/04_assign_plots.R first. Missing: ", assignments_file)
  }
  
  out_dir <- fs::path(project_dir, "data", "processed", "mc_jitter_library")
  replicates_dir <- fs::path(out_dir, "replicates")
  manifest_file <- fs::path(out_dir, "manifest.yml")
  
  # Load config to get hex grids if not provided
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  if (is.null(hex_grids)) {
    hex_grids <- cfg$hex_grids
    if (is.null(hex_grids)) {
      stop("No hex_grids specified. Set in function call or configs/process.yml")
    }
  }
  
  # Validate hex grids
  for (i in seq_along(hex_grids)) {
    grid <- hex_grids[[i]]
    if (is.null(grid$name)) stop("hex_grids[[", i, "]] missing 'name'")
    if (is.null(grid$path)) stop("hex_grids[[", i, "]] missing 'path'")
    if (!fs::file_exists(grid$path)) {
      stop("Hex grid not found: ", grid$path, "\n  Run --create-hex first!")
    }
  }
  
  grid_names <- sapply(hex_grids, function(x) x$name)
  message("→ Multi-scale jitter library with ", length(hex_grids), " hex grids:")
  for (grid in hex_grids) {
    message("    ", grid$name, ": ", grid$path)
  }
  
  # Check for existing library
  if (fs::dir_exists(out_dir) && fs::file_exists(manifest_file) && !overwrite) {
    manifest <- yaml::read_yaml(manifest_file)
    
    # Check if grids match
    existing_grids <- sapply(manifest$hex_grids, function(x) x$name)
    if (!setequal(existing_grids, grid_names)) {
      message("⚠ Existing library has different hex grids:")
      message("    Existing: ", paste(existing_grids, collapse = ", "))
      message("    Requested: ", paste(grid_names, collapse = ", "))
      message("  Use overwrite=TRUE to regenerate")
      stop("Hex grid mismatch")
    }
    
    message("✓ Jitter library already exists:")
    message("    Created: ", manifest$created)
    message("    Replicates: ", manifest$n_replicates)
    message("    Plots: ", manifest$n_plots)
    message("    Hex grids: ", paste(existing_grids, collapse = ", "))
    message("  Use overwrite=TRUE to regenerate")
    return(invisible(out_dir))
  }
  
  # Setup directories
  fs::dir_create(out_dir, recurse = TRUE)
  fs::dir_create(replicates_dir, recurse = TRUE)
  
  message("→ Reading plot assignments...")
  assignments <- readr::read_csv(assignments_file, show_col_types = FALSE)
  
  # FIXED: Filter for plots that have at least one hex assignment
  # Check all hex_id_* columns
  hex_id_cols <- names(assignments)[grepl("^hex_id_", names(assignments))]
  if (length(hex_id_cols) == 0) {
    stop("No hex_id_* columns found in assignments file!")
  }
  
  # Create a filter for plots with at least one valid hex assignment
  has_any_hex <- rowSums(!is.na(assignments[hex_id_cols])) > 0
  
  valid_plots <- assignments |>
    dplyr::filter(has_any_hex,
                  is.finite(lat_original), is.finite(lon_original))
  
  message("  Total plots with hex assignment: ", nrow(valid_plots))
  
  if (!nrow(valid_plots)) stop("No valid plots to jitter!")
  
  use_hex_union <- use_constraints && isTRUE(cfg$mask$use_hex_union)
  use_state <- use_constraints && isTRUE(cfg$mask$use_state_constraint)
  
  # Setup spatial constraints (do once, reuse for all replicates)
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  crs_ref <- sf::st_crs(5070)
  
  # Use FIRST (finest) grid for constraint boundary
  first_grid <- hex_grids[[1]]
  hex_union_5070 <- if (use_hex_union) {
    message("  Building hex union for jitter constraints (using ", first_grid$name, ")...")
    hu <- build_hex_union_5070(first_grid$path, first_grid$layer)
    sf::st_crs(hu) <- crs_ref
    hu
  } else NULL
  
  state_polys_5070 <- if (use_state) {
    message("  Loading state boundaries for jitter constraints...")
    sp <- build_state_polys_5070(cfg$mask$state_geo_path %||% "", 
                                 cfg$mask$state_field %||% "STATEFP")
    sf::st_crs(sp) <- crs_ref
    sp
  } else NULL
  
  # Load ALL hex grids for assignment
  message("→ Loading all hex grids for multi-scale assignment...")
  hex_grids_sf <- lapply(hex_grids, function(grid) {
    message("    Loading ", grid$name, "...")
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
    
    list(name = grid$name, sf = hx_5070)
  })
  
  message("→ Converting plots to spatial points (once)...")
  pts <- sf::st_as_sf(valid_plots, 
                      coords = c("lon_original", "lat_original"), 
                      crs = 4326, remove = FALSE)
  pts_5070 <- sf::st_transform(pts, crs_ref)
  sf::st_crs(pts_5070) <- crs_ref
  
  # Check which replicates already exist
  existing_files <- list.files(replicates_dir, pattern = "^rep_\\d{4}\\.csv$")
  existing_reps <- if (length(existing_files)) {
    as.integer(gsub("^rep_(\\d{4})\\.csv$", "\\1", existing_files))
  } else {
    integer(0)
  }
  
  completed <- length(existing_reps)
  remaining <- setdiff(1:n_replicates, existing_reps)
  
  if (completed > 0 && !overwrite) {
    message(sprintf("→ Found %d existing replicates, resuming from replicate %d", 
                    completed, min(remaining)))
  } else if (overwrite && completed > 0) {
    message("→ Overwrite=TRUE, deleting existing replicates and starting fresh")
    fs::dir_delete(replicates_dir)
    fs::dir_create(replicates_dir)
    remaining <- 1:n_replicates
    completed <- 0
  }
  
  if (length(remaining) == 0) {
    message("✓ All replicates already complete!")
  } else {
    message(sprintf("→ Generating %d jittered coordinate sets with %d hex grids...", 
                    length(remaining), length(hex_grids)))
    message("  Each replicate: jitter → assign to ALL grids → save CSV")
    
    start_time <- Sys.time()
    report_every <- max(1, floor(length(remaining) / 10))
    
    for (idx in seq_along(remaining)) {
      r <- remaining[idx]
      
      if (idx %% report_every == 0 || idx == 1) {
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        if (idx > 1) {
          rate <- idx / elapsed
          eta <- (length(remaining) - idx) / rate
          message(sprintf("  Progress: %d/%d (%.1f%%) - ETA: %.1f min", 
                          completed + idx, n_replicates, 
                          100 * (completed + idx) / n_replicates, 
                          eta / 60))
        } else {
          message(sprintf("  Starting replicate %d/%d", r, n_replicates))
        }
      }
      
      # Generate jitter for this replicate
      if (!is.null(hex_union_5070) || !is.null(state_polys_5070)) {
        pts_j <- constrained_jitter_once(pts_5070, radius_m, hex_union_5070, 
                                         state_polys_5070, 
                                         max_reroll = as.integer(cfg$mask$max_reroll %||% 20))
      } else {
        u <- runif(nrow(pts_5070))
        rr <- sqrt(u) * radius_m
        theta <- runif(nrow(pts_5070), 0, 2 * pi)
        dx <- rr * cos(theta)
        dy <- rr * sin(theta)
        pts_j <- sf::st_set_geometry(pts_5070, sf::st_geometry(pts_5070) + cbind(dx, dy))
        sf::st_crs(pts_j) <- crs_ref
      }
      
      # CRITICAL: Assign jittered points to ALL hex grids
      rep_data <- sf::st_drop_geometry(pts_j) |>
        dplyr::select(CN, STATECD, MEASYEAR, UNITCD, COUNTYCD, PLOT) |>
        dplyr::mutate(replicate_id = r)
      
      for (grid_info in hex_grids_sf) {
        grid_name <- grid_info$name
        grid_sf <- grid_info$sf
        
        # Spatial join for this grid
        joined <- sf::st_join(pts_j, grid_sf["hex_id"], left = TRUE, join = sf::st_intersects)
        
        # Extract hex_id for this grid
        hex_col_name <- paste0("hex_id_", grid_name)
        
        if ("hex_id.y" %in% names(joined)) {
          rep_data[[hex_col_name]] <- joined$hex_id.y
        } else if ("hex_id.x" %in% names(joined)) {
          rep_data[[hex_col_name]] <- joined$hex_id.x
        } else {
          rep_data[[hex_col_name]] <- joined$hex_id
        }
      }
      
      # Extract jittered coordinates back to lat/lon
      pts_j_4326 <- sf::st_transform(pts_j, 4326)
      coords_j <- sf::st_coordinates(pts_j_4326)
      
      rep_data$lon_jittered <- coords_j[, 1]
      rep_data$lat_jittered <- coords_j[, 2]
      
      # Save this replicate as CSV
      rep_file <- fs::path(replicates_dir, sprintf("rep_%04d.csv", r))
      readr::write_csv(rep_data, rep_file)
    }
    
    message(sprintf("✓ Jittering complete in %.1f minutes", 
                    as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
  }
  
  # Write manifest
  manifest <- list(
    created = as.character(Sys.time()),
    n_replicates = n_replicates,
    n_plots = nrow(valid_plots),
    radius_m = radius_m,
    format = "csv",
    constrained = use_constraints,
    used_hex_union = !is.null(hex_union_5070),
    used_state_constraint = !is.null(state_polys_5070),
    hex_grids = lapply(hex_grids, function(g) list(name = g$name, path = g$path)),
    replicates_dir = "replicates",
    columns = c("CN", "STATECD", "MEASYEAR", "UNITCD", "COUNTYCD", "PLOT", 
                paste0("hex_id_", grid_names),
                "replicate_id", "lon_jittered", "lat_jittered")
  )
  
  yaml::write_yaml(manifest, manifest_file)
  message("✓ Wrote: ", manifest_file)
  
  message("\n=== Multi-Scale Jitter Library Summary ===")
  message("  Replicates:     ", n_replicates)
  message("  Plots per rep:  ", nrow(valid_plots))
  message("  Hex grids:      ", length(hex_grids))
  for (grid in hex_grids) {
    message("    - ", grid$name, " (", grid$path, ")")
  }
  message("  Format:         CSV (individual files)")
  message("  Location:       ", replicates_dir)
  message("  Constrained:    ", use_constraints)
  message("\nColumns in each CSV:")
  message("  ", paste(manifest$columns, collapse = ", "))
  
  invisible(list(
    library_dir = out_dir,
    replicates_dir = replicates_dir,
    manifest = manifest_file
  ))
}

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  overwrite <- "--overwrite" %in% args
  
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  stage4_build_jitter_library(
    project_dir = cfg$project_dir %||% ".",
    hex_grids = cfg$hex_grids,  # Load from config
    n_replicates = cfg$mc_reps %||% 100,
    radius_m = cfg$jitter_radius_m %||% 1609.34, 
    use_constraints = TRUE,
    overwrite = overwrite
  )
}