# R/05_build_jitter_library.R
# Pre-generate N jittered coordinate sets with hex assignment and incremental saving

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(fs); library(yaml); library(glue)
})

source("R/utils_spatial.R")

`%||%` <- function(a,b) if (!is.null(a)) a else b

stage4_build_jitter_library <- function(project_dir = ".",
                                        hex_path = "data/hex/hex_grid.geojson",
                                        hex_layer = NULL,
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
  
  # Check for existing library
  if (fs::dir_exists(out_dir) && fs::file_exists(manifest_file) && !overwrite) {
    manifest <- yaml::read_yaml(manifest_file)
    message("✓ Jitter library already exists:")
    message("    Created: ", manifest$created)
    message("    Replicates: ", manifest$n_replicates)
    message("    Plots: ", manifest$n_plots)
    message("  Use overwrite=TRUE to regenerate")
    return(invisible(out_dir))
  }
  
  # Setup directories
  fs::dir_create(out_dir, recurse = TRUE)
  fs::dir_create(replicates_dir, recurse = TRUE)
  
  message("→ Reading plot assignments...")
  assignments <- readr::read_csv(assignments_file, show_col_types = FALSE)
  
  valid_plots <- assignments |>
    dplyr::filter(!is.na(hex_id), 
                  is.finite(lat_original), is.finite(lon_original))
  
  message("  Total plots with hex assignment: ", nrow(valid_plots))
  
  if (!nrow(valid_plots)) stop("No valid plots to jitter!")
  
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  use_hex_union <- use_constraints && isTRUE(cfg$mask$use_hex_union)
  use_state <- use_constraints && isTRUE(cfg$mask$use_state_constraint)
  
  # Setup spatial constraints (do once, reuse for all replicates)
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  crs_ref <- sf::st_crs(5070)
  
  hex_union_5070 <- if (use_hex_union) {
    message("  Building hex union for jitter constraints...")
    hu <- build_hex_union_5070(hex_path, hex_layer)
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
  
  # Read hex grid for spatial joins
  message("→ Reading hex grid for assignment...")
  hx <- if (is.null(hex_layer)) {
    sf::st_read(hex_path, quiet = TRUE)
  } else {
    sf::st_read(hex_path, layer = hex_layer, quiet = TRUE)
  }
  
  if (!("hex_id" %in% names(hx))) {
    if ("ID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = ID)
    else if ("OBJECTID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = OBJECTID)
    else hx$hex_id <- seq_len(nrow(hx))
  }
  hx$hex_id <- as.character(hx$hex_id)
  hx <- sf::st_make_valid(sf::st_zm(hx, drop = TRUE, what = "ZM"))
  hx_5070 <- sf::st_transform(hx, crs_ref)
  
  message("→ Converting to spatial points (once)...")
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
    message(sprintf("→ Generating %d jittered coordinate sets...", length(remaining)))
    message("  Each replicate: jitter → assign hex → save CSV")
    
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
      
      # CRITICAL: Assign jittered points to hexes (may be different hex!)
      joined <- sf::st_join(pts_j, hx_5070["hex_id"], left = TRUE, join = sf::st_intersects)
      
      # Handle potential duplicate hex_id columns from join
      if ("hex_id.y" %in% names(joined)) {
        joined <- joined |>
          dplyr::mutate(hex_id_jittered = hex_id.y) |>
          dplyr::select(-hex_id.y, -hex_id.x)
      } else if ("hex_id.x" %in% names(joined)) {
        joined <- joined |>
          dplyr::mutate(hex_id_jittered = hex_id.x) |>
          dplyr::select(-hex_id.x)
      } else {
        joined <- joined |>
          dplyr::mutate(hex_id_jittered = hex_id)
      }
      
      # Extract jittered coordinates back to lat/lon
      pts_j_4326 <- sf::st_transform(joined, 4326)
      coords_j <- sf::st_coordinates(pts_j_4326)
      
      # Create replicate data (NO GEOMETRY - just tabular CSV)
      rep_data <- sf::st_drop_geometry(pts_j_4326) |>
        dplyr::select(CN, STATECD, MEASYEAR, UNITCD, COUNTYCD, PLOT, hex_id_jittered) |>
        dplyr::mutate(
          replicate_id = r,
          lon_jittered = coords_j[, 1],
          lat_jittered = coords_j[, 2]
        ) |>
        as.data.frame()
      
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
    replicates_dir = "replicates",
    columns = c("CN", "STATECD", "MEASYEAR", "UNITCD", "COUNTYCD", "PLOT", 
                "hex_id_jittered", "replicate_id", "lon_jittered", "lat_jittered")
  )
  
  yaml::write_yaml(manifest, manifest_file)
  message("✓ Wrote: ", manifest_file)
  
  message("\n=== Jitter Library Summary ===")
  message("  Replicates:     ", n_replicates)
  message("  Plots per rep:  ", nrow(valid_plots))
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
    hex_path = cfg$hex_path %||% "data/hex/hex_grid.geojson",
    hex_layer = cfg$hex_layer %||% NULL,
    n_replicates = cfg$mc_reps %||% 100,
    radius_m = cfg$jitter_radius_m %||% 1609.34,
    use_constraints = TRUE,
    overwrite = overwrite
  )
}