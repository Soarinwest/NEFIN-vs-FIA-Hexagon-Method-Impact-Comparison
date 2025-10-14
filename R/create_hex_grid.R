# R/create_hex_grid.R
# Generate hexagonal grid at specified size and clip to study region

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(fs); library(units)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

#' Create hexagonal grid
#' 
#' @param area_acres Target area per hexagon in acres
#' @param clip_to Path to polygon to clip grid to (GeoJSON, shapefile, etc.)
#' @param clip_layer Layer name if clip_to is multi-layer (e.g., geopackage)
#' @param exclude_water Path to waterbodies polygon to exclude (optional)
#' @param buffer_miles Buffer around clipping region (to catch edge hexes)
#' @param crs Target CRS (default: 5070 - Albers Equal Area)
#' @param out_path Output path for hex grid GeoJSON
#' @param overwrite Overwrite existing output
#' 
#' @return Path to output file
create_hex_grid <- function(area_acres = 6000,
                            clip_to = "data/hex/hex_grid.geojson",
                            clip_layer = NULL,
                            exclude_water = NULL,
                            buffer_miles = 5,
                            crs = 5070,
                            out_path = NULL,
                            overwrite = FALSE) {
  
  message("\n═══════════════════════════════════════════════════════════")
  message("HEXAGONAL GRID GENERATION")
  message("═══════════════════════════════════════════════════════════\n")
  
  # Auto-generate output path if not provided
  if (is.null(out_path)) {
    size_label <- if (area_acres >= 1000) {
      paste0(round(area_acres / 1000, 1), "k")
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
  
  # Transform to target CRS
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
  
  # Add buffer
  buffer_m <- buffer_miles * 1609.34
  message("→ Buffering by ", buffer_miles, " miles (", round(buffer_m), " meters)")
  boundary_buffered <- sf::st_buffer(boundary, buffer_m)
  
  # Calculate hex cell size
  # For regular hexagon: area = (3 * sqrt(3) / 2) * side^2
  # We want flat-topped hexagons, cellsize is the distance across (2 * side * cos(30°))
  
  area_m2 <- area_acres * 4046.86  # acres to m²
  side_m <- sqrt(area_m2 / (3 * sqrt(3) / 2))
  cellsize_m <- side_m * sqrt(3)  # flat-topped hex width
  
  # Convert to miles for reporting
  side_miles <- side_m / 1609.34
  cellsize_miles <- cellsize_m / 1609.34
  
  message("\n→ Hexagon parameters:")
  message("  Target area: ", format(area_acres, big.mark = ","), " acres")
  message("  Hex side: ", round(side_m, 1), " m (", round(side_miles, 3), " miles)")
  message("  Cell size (width): ", round(cellsize_m, 1), " m (", round(cellsize_miles, 3), " miles)")
  
  # Create hex grid
  message("\n→ Generating hexagonal grid...")
  hex_grid <- sf::st_make_grid(
    boundary_buffered,
    cellsize = cellsize_m,
    square = FALSE,
    flat_topped = TRUE
  )
  
  message("  Generated ", length(hex_grid), " hexagons (before clipping)")
  
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
  hex_clipped$area <- as.numeric(sf::st_area(hex_clipped))
  target_area_m2 <- area_m2
  
  initial_count <- nrow(hex_clipped)
  hex_clipped <- hex_clipped |>
    dplyr::filter(area > target_area_m2 * 0.5) |>
    dplyr::select(-area)
  
  removed_slivers <- initial_count - nrow(hex_clipped)
  if (removed_slivers > 0) {
    message("  Removed ", removed_slivers, " sliver hexes (<50% coverage)")
  }
  message("  Kept ", nrow(hex_clipped), " hexagons after clipping")
  
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
  
  message("\n═══════════════════════════════════════════════════════════")
  message("SUMMARY")
  message("═══════════════════════════════════════════════════════════")
  message("  Output: ", out_path)
  message("  Hexagons: ", nrow(hex_out))
  if (!is.null(water_mask)) {
    message("  Waterbodies excluded: YES")
  }
  message("  Target area: ", format(area_acres, big.mark = ","), " acres")
  message("  Actual area (mean): ", round(mean(actual_areas_acres), 1), " acres")
  message("  Actual area (range): ", round(min(actual_areas_acres), 1), " - ", 
          round(max(actual_areas_acres), 1), " acres")
  message("  Hex spacing: ~", round(cellsize_miles, 2), " miles")
  message("═══════════════════════════════════════════════════════════\n")
  
  invisible(out_path)
}

# Helper: Create multiple grid sizes at once
create_multiple_grids <- function(area_acres_list = c(1500, 3000, 6000, 12000),
                                  clip_to = "data/hex/hex_grid.geojson",
                                  clip_layer = NULL,
                                  exclude_water = NULL,
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
  
  # Check for multi-grid mode
  if ("--multiple" %in% args) {
    # Parse area list
    areas_str <- get_arg("--areas", "1500,3000,6000,12000")
    areas <- as.numeric(strsplit(areas_str, ",")[[1]])
    
    create_multiple_grids(
      area_acres_list = areas,
      clip_to = get_arg("--clip", "data/hex/hex_grid.geojson"),
      clip_layer = get_arg("--layer", NULL),
      exclude_water = get_arg("--exclude-water", NULL),
      overwrite = "--overwrite" %in% args
    )
  } else {
    # Single grid
    area <- as.numeric(get_arg("--area", "6000"))
    
    create_hex_grid(
      area_acres = area,
      clip_to = get_arg("--clip", "data/hex/hex_grid.geojson"),
      clip_layer = get_arg("--layer", NULL),
      exclude_water = get_arg("--exclude-water", NULL),
      buffer_miles = as.numeric(get_arg("--buffer", "5")),
      out_path = get_arg("--out", NULL),
      overwrite = "--overwrite" %in% args
    )
  }
}