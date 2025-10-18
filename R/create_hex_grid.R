# R/create_hex_grid.R
# Generate hexagonal grid at specified size and clip to study region

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(fs); library(units)
  library(parallel); library(doParallel); library(foreach)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

# Estimate processing time based on area
estimate_hex_time <- function(area_acres, n_features = 1000) {
  # Rough estimates based on hex size (larger hexes = fewer hexes = faster)
  # These are approximations and will vary by system
  if (area_acres < 1000) {
    base_time <- 180  # 3 minutes for small hexes
  } else if (area_acres < 10000) {
    base_time <- 120  # 2 minutes for medium hexes
  } else if (area_acres < 100000) {
    base_time <- 60   # 1 minute for large hexes
  } else if (area_acres < 1000000) {
    base_time <- 30   # 30 seconds for very large hexes
  } else {
    base_time <- 15   # 15 seconds for huge hexes
  }
  
  # Adjust for number of features in clip region
  complexity_factor <- sqrt(n_features / 1000)
  estimated_seconds <- base_time * complexity_factor
  
  return(estimated_seconds)
}

format_time <- function(seconds) {
  if (seconds < 60) {
    return(paste0(round(seconds), " seconds"))
  } else if (seconds < 3600) {
    return(paste0(round(seconds / 60, 1), " minutes"))
  } else {
    return(paste0(round(seconds / 3600, 1), " hours"))
  }
}

#' Create hexagonal grid
#' 
#' @param area_acres Target area per hexagon in acres
#' @param clip_to Path to polygon to clip grid to (GeoJSON, shapefile, etc.)
#' @param clip_layer Layer name if clip_to is multi-layer (e.g., geopackage)
#' @param exclude_water Path to waterbodies polygon to exclude (optional)
#' @param buffer_miles Buffer around clipping region (to catch edge hexes) - set to NULL for no buffer
#' @param crs Target CRS (default: 5070 - Albers Equal Area)
#' @param out_path Output path for hex grid GeoJSON
#' @param overwrite Overwrite existing output
#' 
#' @return Path to output file
create_hex_grid <- function(area_acres = 6000,
                            clip_to = "data/hex/hex_grid.geojson",
                            clip_layer = NULL,
                            exclude_water = NULL,
                            buffer_miles = NULL,  # Default to NULL (no buffer)
                            crs = 5070,
                            out_path = NULL,
                            overwrite = FALSE) {
  
  start_time <- Sys.time()
  
  message("\n═══════════════════════════════════════════════════════════")
  message("HEXAGONAL GRID GENERATION")
  message("═══════════════════════════════════════════════════════════\n")
  
  # Auto-generate output path if not provided
  if (is.null(out_path)) {
    size_label <- if (area_acres >= 1000) {
      paste0(round(area_acres / 1000), "k")
    } else {
      area_acres
    }
    out_path <- fs::path("data", "hex", paste0("hex_grid_", size_label, "ac.geojson"))
  }
  
  if (fs::file_exists(out_path) && !overwrite) {
    message("✓ Output already exists: ", out_path)
    message("  Use overwrite=TRUE to regenerate")
    return(invisible(out_path))
  }
  
  # Turn off spherical geometry for better performance
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  # Read clipping region
  message("→ Reading clipping region: ", clip_to)
  clip_poly <- if (is.null(clip_layer)) {
    sf::st_read(clip_to, quiet = TRUE)
  } else {
    sf::st_read(clip_to, layer = clip_layer, quiet = TRUE)
  }
  
  clip_poly <- sf::st_make_valid(clip_poly)
  message("  Loaded ", nrow(clip_poly), " features")
  
  # Estimate processing time
  est_time <- estimate_hex_time(area_acres, nrow(clip_poly))
  message("\n⏱ Estimated processing time: ", format_time(est_time))
  
  # Transform to target CRS
  message("\n→ Transforming to target CRS...")
  clip_poly_proj <- sf::st_transform(clip_poly, crs)
  
  # Union all features into single boundary
  message("→ Creating boundary union...")
  boundary <- sf::st_union(clip_poly_proj) |> sf::st_make_valid()
  
  # Handle waterbodies exclusion
  water_mask <- NULL
  if (!is.null(exclude_water) && nzchar(exclude_water)) {
    if (fs::file_exists(exclude_water)) {
      message("→ Loading waterbodies to exclude: ", exclude_water)
      water <- sf::st_read(exclude_water, quiet = TRUE)
      water <- sf::st_make_valid(water)
      message("  Loaded ", nrow(water), " waterbody features")
      
      # Transform to target CRS
      water_proj <- sf::st_transform(water, crs)
      
      # Union all water polygons
      message("→ Creating waterbodies mask...")
      water_mask <- sf::st_union(water_proj) |> sf::st_make_valid()
      
      # Get water area stats
      water_area_m2 <- as.numeric(sf::st_area(water_mask))
      water_area_acres <- water_area_m2 / 4046.86
      message("  Total water area: ", format(round(water_area_acres), big.mark = ","), " acres")
      
      # Subtract water from land boundary
      message("→ Subtracting waterbodies from land area...")
      boundary <- sf::st_difference(boundary, water_mask) |> sf::st_make_valid()
      
      land_area_m2 <- as.numeric(sf::st_area(boundary))
      land_area_acres <- land_area_m2 / 4046.86
      message("  Remaining land area: ", format(round(land_area_acres), big.mark = ","), " acres")
      message("  Water coverage: ", round(100 * water_area_acres / (water_area_acres + land_area_acres), 1), "%")
    } else {
      message("⚠ Waterbodies file not found: ", exclude_water)
      message("  Continuing without water exclusion...")
    }
  }
  
  # Add buffer only if buffer_miles is not NULL and greater than 0
  if (!is.null(buffer_miles) && buffer_miles > 0) {
    buffer_m <- buffer_miles * 1609.34
    message("→ Buffering by ", buffer_miles, " miles (", round(buffer_m), " meters)")
    boundary_buffered <- sf::st_buffer(boundary, buffer_m)
  } else {
    message("→ No buffer applied (buffer_miles = ", buffer_miles, ")")
    boundary_buffered <- boundary
  }
  
  # Calculate hex cell size
  area_m2 <- area_acres * 4046.86  # acres to m²
  side_m <- sqrt(area_m2 / (3 * sqrt(3) / 2))
  cellsize_m <- side_m * sqrt(3)  # flat-topped hex width
  
  # Convert to miles for reporting
  side_miles <- side_m / 1609.34
  cellsize_miles <- cellsize_m / 1609.34
  
  # Estimate number of hexagons
  bbox_area <- as.numeric(sf::st_area(sf::st_as_sfc(sf::st_bbox(boundary_buffered))))
  est_n_hexes <- round(bbox_area / area_m2 * 1.2)  # 1.2 factor for packing efficiency
  
  message("\n→ Hexagon parameters:")
  message("  Target area: ", format(area_acres, big.mark = ","), " acres")
  message("  Hex side: ", round(side_m, 1), " m (", round(side_miles, 3), " miles)")
  message("  Cell size (width): ", round(cellsize_m, 1), " m (", round(cellsize_miles, 3), " miles)")
  message("  Estimated hexagons: ~", format(est_n_hexes, big.mark = ","))
  
  # Create hex grid
  message("\n→ Generating hexagonal grid...")
  hex_grid <- sf::st_make_grid(
    boundary_buffered,
    cellsize = cellsize_m,
    square = FALSE,
    flat_topped = TRUE
  )
  
  message("  Generated ", format(length(hex_grid), big.mark = ","), " hexagons (before clipping)")
  
  # Convert to sf object
  hex_sf <- sf::st_sf(
    hex_id = sprintf("H%06d", seq_along(hex_grid)),
    geometry = hex_grid,
    crs = crs
  )
  
  # Clip to original boundary (no buffer) - this already has water removed
  message("→ Clipping to study region (land only)...")
  hex_clipped <- sf::st_intersection(hex_sf, boundary)
  
  # Keep only hexes with >50% overlap (avoid slivers)
  message("→ Filtering sliver hexagons...")
  hex_clipped$area <- as.numeric(sf::st_area(hex_clipped))
  target_area_m2 <- area_m2
  
  initial_count <- nrow(hex_clipped)
  hex_clipped <- hex_clipped |>
    dplyr::filter(area > target_area_m2 * 0.5) |>
    dplyr::select(-area)
  
  removed_slivers <- initial_count - nrow(hex_clipped)
  if (removed_slivers > 0) {
    message("  Removed ", format(removed_slivers, big.mark = ","), " sliver hexes (<50% coverage)")
  }
  message("  Kept ", format(nrow(hex_clipped), big.mark = ","), " hexagons after clipping")
  
  # Re-assign sequential IDs
  hex_clipped$hex_id <- sprintf("H%06d", seq_len(nrow(hex_clipped)))
  
  # Transform to WGS84 for output
  message("→ Transforming to WGS84 for GeoJSON output...")
  hex_out <- sf::st_transform(hex_clipped, 4326)
  
  # Write output
  fs::dir_create(fs::path_dir(out_path), recurse = TRUE)
  message("→ Writing to: ", out_path)
  sf::st_write(hex_out, out_path, delete_dsn = TRUE, quiet = TRUE)
  
  # Summary statistics
  actual_areas_m2 <- as.numeric(sf::st_area(sf::st_transform(hex_out, crs)))
  actual_areas_acres <- actual_areas_m2 / 4046.86
  
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  message("\n═══════════════════════════════════════════════════════════")
  message("SUMMARY")
  message("═══════════════════════════════════════════════════════════")
  message("  Output: ", out_path)
  message("  Hexagons: ", format(nrow(hex_out), big.mark = ","))
  if (!is.null(water_mask)) {
    message("  Waterbodies excluded: YES")
  }
  message("  Buffer applied: ", if(!is.null(buffer_miles) && buffer_miles > 0) paste0(buffer_miles, " miles") else "NO")
  message("  Target area: ", format(area_acres, big.mark = ","), " acres")
  message("  Actual area (mean): ", round(mean(actual_areas_acres), 1), " acres")
  message("  Actual area (range): ", round(min(actual_areas_acres), 1), " - ", 
          round(max(actual_areas_acres), 1), " acres")
  message("  Hex spacing: ~", round(cellsize_miles, 2), " miles")
  message("  Processing time: ", format_time(elapsed_time))
  message("═══════════════════════════════════════════════════════════\n")
  
  invisible(out_path)
}

# Helper: Create multiple grid sizes in parallel
create_multiple_grids_parallel <- function(hex_sizes_specs,
                                           n_cores = NULL,
                                           overwrite = FALSE) {
  
  # Define the %||% function here so it's available in this scope
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  
  if (is.null(n_cores)) {
    n_cores <- min(parallel::detectCores() - 1, length(hex_sizes_specs), 4)
  }
  
  message("\n════════════════════════════════════════════════════════════")
  message("PARALLEL HEX GRID CREATION")
  message("════════════════════════════════════════════════════════════")
  message("  Grids to create: ", length(hex_sizes_specs))
  message("  Using ", n_cores, " cores")
  
  # Estimate total time
  total_est_time <- 0
  for (spec in hex_sizes_specs) {
    total_est_time <- total_est_time + estimate_hex_time(spec$area_acres)
  }
  parallel_est_time <- total_est_time / n_cores
  
  message("  Estimated time (sequential): ", format_time(total_est_time))
  message("  Estimated time (parallel): ", format_time(parallel_est_time))
  message("════════════════════════════════════════════════════════════\n")
  
  # Check which grids already exist
  to_create <- list()
  already_exist <- list()
  
  for (spec in hex_sizes_specs) {
    if (fs::file_exists(spec$out_path) && !overwrite) {
      already_exist[[length(already_exist) + 1]] <- spec
    } else {
      to_create[[length(to_create) + 1]] <- spec
    }
  }
  
  if (length(already_exist) > 0) {
    message("Skipping ", length(already_exist), " existing grids:")
    for (spec in already_exist) {
      message("  ✓ ", spec$out_path)
    }
    message("")
  }
  
  if (length(to_create) == 0) {
    message("All grids already exist! Use overwrite=TRUE to regenerate.")
    return(invisible(NULL))
  }
  
  message("Creating ", length(to_create), " new grids...\n")
  
  # Set up parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  
  # Export the %||% function to workers
  clusterEvalQ(cl, {
    `%||%` <- function(a,b) if (!is.null(a)) a else b
  })
  
  # Run in parallel
  results <- foreach(spec = to_create,
                     .packages = c("sf", "dplyr", "fs", "units"),
                     .export = c("create_hex_grid", "estimate_hex_time", 
                                 "format_time")) %dopar% {
                                   
                                   create_hex_grid(
                                     area_acres = spec$area_acres,
                                     clip_to = spec$clip_to,
                                     clip_layer = spec$clip_layer %||% NULL,
                                     exclude_water = spec$exclude_water %||% NULL,
                                     buffer_miles = spec$buffer_miles %||% NULL,  # Use NULL if not specified
                                     crs = spec$crs %||% 5070,
                                     out_path = spec$out_path,
                                     overwrite = overwrite
                                   )
                                 }
  
  message("\n════════════════════════════════════════════════════════════")
  message("PARALLEL PROCESSING COMPLETE")
  message("════════════════════════════════════════════════════════════")
  message("  Created: ", length(results), " grids")
  message("════════════════════════════════════════════════════════════\n")
  
  invisible(results)
}

# Helper: Create multiple grid sizes at once (sequential)
create_multiple_grids <- function(area_acres_list = c(1500, 3000, 6000, 12000),
                                  clip_to = "data/hex/hex_grid.geojson",
                                  clip_layer = NULL,
                                  exclude_water = NULL,
                                  buffer_miles = NULL,
                                  crs = 5070,
                                  overwrite = FALSE) {
  
  results <- list()
  
  for (area in area_acres_list) {
    message("\n\n")
    message("████████████████████████████████████████████████████████████")
    message("Creating grid: ", format(area, big.mark = ","), " acres")
    message("████████████████████████████████████████████████████████████")
    
    out_path <- create_hex_grid(
      area_acres = area,
      clip_to = clip_to,
      clip_layer = clip_layer,
      exclude_water = exclude_water,
      buffer_miles = buffer_miles,
      crs = crs,
      overwrite = overwrite
    )
    
    results[[as.character(area)]] <- out_path
  }
  
  message("\n\n")
  message("═══════════════════════════════════════════════════════════")
  message("ALL GRIDS CREATED")
  message("═══════════════════════════════════════════════════════════")
  for (area in names(results)) {
    message("  ", area, " acres: ", results[[area]])
  }
  message("═══════════════════════════════════════════════════════════\n")
  
  invisible(results)
}

# CLI interface
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  get_arg <- function(flag, default = NULL) {
    hit <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
  }
  
  # Check for parallel mode
  if ("--parallel" %in% args) {
    # Load config
    if (file.exists("configs/process.yml")) {
      suppressPackageStartupMessages(library(yaml))
      cfg <- yaml::read_yaml("configs/process.yml")
      if (!is.null(cfg$hex_sizes)) {
        n_cores <- as.numeric(get_arg("--cores", NULL))
        create_multiple_grids_parallel(
          hex_sizes_specs = cfg$hex_sizes,
          n_cores = n_cores,
          overwrite = "--overwrite" %in% args
        )
      } else {
        stop("No hex_sizes defined in configs/process.yml")
      }
    } else {
      stop("configs/process.yml not found")
    }
  } else if ("--multiple" %in% args) {
    # Sequential multi-grid mode
    areas_str <- get_arg("--areas", "1500,3000,6000,12000")
    areas <- as.numeric(strsplit(areas_str, ",")[[1]])
    
    buffer_arg <- get_arg("--buffer", NULL)
    buffer_miles <- if (!is.null(buffer_arg)) as.numeric(buffer_arg) else NULL
    
    create_multiple_grids(
      area_acres_list = areas,
      clip_to = get_arg("--clip", "data/hex/hex_grid.geojson"),
      clip_layer = get_arg("--layer", NULL),
      exclude_water = get_arg("--exclude-water", NULL),
      buffer_miles = buffer_miles,
      overwrite = "--overwrite" %in% args
    )
  } else {
    # Single grid mode
    area <- as.numeric(get_arg("--area", "6000"))
    
    buffer_arg <- get_arg("--buffer", NULL)
    buffer_miles <- if (!is.null(buffer_arg)) as.numeric(buffer_arg) else NULL
    
    create_hex_grid(
      area_acres = area,
      clip_to = get_arg("--clip", "data/hex/hex_grid.geojson"),
      clip_layer = get_arg("--layer", NULL),
      exclude_water = get_arg("--exclude-water", NULL),
      buffer_miles = buffer_miles,
      out_path = get_arg("--out", NULL),
      overwrite = "--overwrite" %in% args
    )
  }
}