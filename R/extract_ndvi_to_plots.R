# R/extract_ndvi_to_plots.R
# ============================================================================
# Purpose: Extract NDVI values from MODIS and Sentinel-2 rasters to plot locations
#          Handles both FIA (jittered coordinates) and NEFIN (true coordinates)
#
# Inputs:
#   - NDVI GeoTIFFs (MODIS and/or S2) in data/processed/ndvi/
#   - FIA jitter library: data/processed/mc_jitter_library/
#   - NEFIN processed: data/processed/nefin_processed.csv
#   - Plot hex assignments: data/processed/plot_hex_assignments.csv
#   - Configuration: configs/process.yml
#
# Outputs:
#   - data/processed/ndvi/fia_ndvi_extracted.csv (FIA with NDVI by replicate)
#   - data/processed/ndvi/nefin_ndvi_extracted.csv (NEFIN with NDVI)
#   - data/processed/ndvi/extraction_manifest.yml (metadata)
#
# Usage:
#   Rscript R/extract_ndvi_to_plots.R [--overwrite] [--modis-only] [--s2-only]
#
# Pipeline integration:
#   Run after: NDVI rasters exported from GEE, jitter library built
#   Run before: R/build_ndvi_plot_dataset.R
# ============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(readr)
  library(fs)
  library(yaml)
  library(glue)
  library(tidyr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ============================================================================
# Configuration Loading
# ============================================================================

load_ndvi_config <- function(cfg_path = "configs/process.yml") {
  if (!fs::file_exists(cfg_path)) {
    stop("Configuration file not found: ", cfg_path)
  }
  
  cfg <- yaml::read_yaml(cfg_path)
  
  # Set defaults for NDVI section if missing
  ndvi_cfg <- cfg$ndvi %||% list()
  
  ndvi_cfg$use_modis <- ndvi_cfg$use_modis %||% TRUE
  ndvi_cfg$use_s2 <- ndvi_cfg$use_s2 %||% TRUE
  ndvi_cfg$modis_dir <- ndvi_cfg$modis_dir %||% "data/processed/ndvi/modis"
  ndvi_cfg$s2_dir <- ndvi_cfg$s2_dir %||% "data/processed/ndvi/s2"
  ndvi_cfg$modis_pattern <- ndvi_cfg$modis_pattern %||% "MODIS_NDVI_5yr_blocked_{start}_{end}.tif"
  ndvi_cfg$s2_pattern <- ndvi_cfg$s2_pattern %||% "S2_NDVI_10m_5yr_blocked_{start}_{end}.tif"
  ndvi_cfg$extraction_buffer_m <- ndvi_cfg$extraction_buffer_m %||% 0
  ndvi_cfg$modis_extraction_method <- ndvi_cfg$modis_extraction_method %||% "point"
  ndvi_cfg$s2_extraction_method <- ndvi_cfg$s2_extraction_method %||% "point"
  
  list(
    ndvi = ndvi_cfg,
    years = cfg$years,
    level_window = cfg$level_window %||% 5,
    hex_grids = cfg$hex_grids,
    project_dir = cfg$project_dir %||% "."
  )
}

# ============================================================================
# NDVI Raster Discovery
# ============================================================================

#' Find available NDVI rasters matching configuration
#' 
#' @param ndvi_cfg NDVI configuration list
#' @return data.frame with columns: product, start_year, end_year, path, exists
find_ndvi_rasters <- function(ndvi_cfg) {
  
  rasters <- list()
  
  # Check MODIS rasters
  if (ndvi_cfg$use_modis && !is.null(ndvi_cfg$years_to_windows)) {
    for (window_name in names(ndvi_cfg$years_to_windows)) {
      window <- ndvi_cfg$years_to_windows[[window_name]]
      
      if (isTRUE(window$has_modis)) {
        filename <- glue::glue(ndvi_cfg$modis_pattern,
                               start = window$start,
                               end = window$end)
        filepath <- fs::path(ndvi_cfg$modis_dir, filename)
        
        rasters[[length(rasters) + 1]] <- data.frame(
          product = "modis",
          start_year = window$start,
          end_year = window$end,
          window_name = window_name,
          path = as.character(filepath),
          exists = fs::file_exists(filepath),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Check Sentinel-2 rasters
  if (ndvi_cfg$use_s2 && !is.null(ndvi_cfg$years_to_windows)) {
    for (window_name in names(ndvi_cfg$years_to_windows)) {
      window <- ndvi_cfg$years_to_windows[[window_name]]
      
      if (isTRUE(window$has_s2)) {
        filename <- glue::glue(ndvi_cfg$s2_pattern,
                               start = window$start,
                               end = window$end)
        filepath <- fs::path(ndvi_cfg$s2_dir, filename)
        
        rasters[[length(rasters) + 1]] <- data.frame(
          product = "s2",
          start_year = window$start,
          end_year = window$end,
          window_name = window_name,
          path = as.character(filepath),
          exists = fs::file_exists(filepath),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(rasters) == 0) {
    return(data.frame(
      product = character(),
      start_year = integer(),
      end_year = integer(),
      window_name = character(),
      path = character(),
      exists = logical()
    ))
  }
  
  dplyr::bind_rows(rasters)
}

#' Map measurement year to appropriate NDVI window
#' 
#' @param measyear Measurement year (integer)
#' @param level_window Window size in years
#' @param years_to_windows Configuration mapping
#' @return window_name string or NA
map_year_to_window <- function(measyear, level_window, years_to_windows) {
  if (is.null(years_to_windows)) return(NA_character_)
  
  # Find window that contains this measurement year
  for (window_name in names(years_to_windows)) {
    window <- years_to_windows[[window_name]]
    if (measyear >= window$start && measyear <= window$end) {
      return(window_name)
    }
  }
  
  NA_character_
}

# ============================================================================
# NDVI Extraction Functions
# ============================================================================

#' Extract NDVI values from a raster to points
#' 
#' @param points_sf sf object with point geometries
#' @param raster_path Path to GeoTIFF
#' @param method "point" or "buffer"
#' @param buffer_m Buffer radius in meters (if method="buffer")
#' @return numeric vector of NDVI values
extract_ndvi_values <- function(points_sf, raster_path, method = "point", buffer_m = 0) {
  
  if (!fs::file_exists(raster_path)) {
    warning("Raster not found: ", raster_path)
    return(rep(NA_real_, nrow(points_sf)))
  }
  
  # Load raster
  r <- terra::rast(raster_path)
  
  # Convert sf to terra vect
  pts_vect <- terra::vect(points_sf)
  
  # Ensure same CRS
  if (!terra::same.crs(pts_vect, r)) {
    pts_vect <- terra::project(pts_vect, terra::crs(r))
  }
  
  # Extract values
  if (method == "buffer" && buffer_m > 0) {
    # Create buffers and extract mean
    vals <- terra::extract(r, pts_vect, buffer = buffer_m, fun = mean, na.rm = TRUE)
  } else {
    # Point extraction
    vals <- terra::extract(r, pts_vect)
  }
  
  # Return values (second column is the value, first is ID)
  if (ncol(vals) >= 2) {
    return(as.numeric(vals[[2]]))
  } else {
    return(rep(NA_real_, nrow(points_sf)))
  }
}

# ============================================================================
# FIA Extraction (with jitter replicates)
# ============================================================================

#' Extract NDVI for FIA plots using jitter library
#' 
#' @param cfg Configuration list
#' @param available_rasters data.frame from find_ndvi_rasters()
#' @param max_reps Maximum number of replicates to process (NULL = all)
#' @return data.frame with columns: CN, STATECD, MEASYEAR, replicate_id, 
#'         ndvi_modis, ndvi_s2, window_name, hex_id_*
extract_fia_ndvi <- function(cfg, available_rasters, max_reps = NULL) {
  
  message("\n═══════════════════════════════════════════════════════════")
  message("EXTRACTING NDVI FOR FIA PLOTS (JITTERED COORDINATES)")
  message("═══════════════════════════════════════════════════════════\n")
  
  ndvi_cfg <- cfg$ndvi
  jitter_dir <- fs::path(cfg$project_dir, "data", "processed", "mc_jitter_library")
  replicates_dir <- fs::path(jitter_dir, "replicates")
  
  if (!fs::dir_exists(replicates_dir)) {
    stop("Jitter library not found: ", replicates_dir,
         "\n  Run: Rscript run_pipeline.R --jitter")
  }
  
  # Get replicate files
  rep_files <- list.files(replicates_dir, pattern = "^rep_\\d{4}\\.csv$", full.names = TRUE)
  
  if (!is.null(max_reps) && max_reps < length(rep_files)) {
    rep_files <- rep_files[1:max_reps]
  }
  
  message("→ Processing ", length(rep_files), " jitter replicates")
  
  # Filter to existing rasters only
  available_rasters <- available_rasters |> dplyr::filter(exists)
  
  if (nrow(available_rasters) == 0) {
    warning("No NDVI rasters found. Check paths in configs/process.yml")
    return(NULL)
  }
  
  message("→ Available NDVI products:")
  for (i in seq_len(nrow(available_rasters))) {
    r <- available_rasters[i, ]
    message("    ", r$product, " ", r$start_year, "-", r$end_year, ": ", r$path)
  }
  
  # Load plot-level FIA data for MEASYEAR
  fia_root <- if (fs::dir_exists(fs::path(cfg$project_dir, "data", "interim", "fia_region"))) {
    fs::path(cfg$project_dir, "data", "interim", "fia_region")
  } else {
    fs::path(cfg$project_dir, "data", "interim", "fia")
  }
  
  plot_csv <- fs::path(fia_root, "plot.csv")
  if (!fs::file_exists(plot_csv)) {
    stop("FIA plot data not found: ", plot_csv)
  }
  
  plot_data <- readr::read_csv(plot_csv, show_col_types = FALSE) |>
    dplyr::select(CN, STATECD, MEASYEAR) |>
    dplyr::distinct()
  
  message("→ FIA plots with measurement years: ", nrow(plot_data))
  
  # Turn off s2 for sf operations

old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  # Process each replicate
  all_results <- list()
  
  start_time <- Sys.time()
  
  for (idx in seq_along(rep_files)) {
    rep_file <- rep_files[idx]
    
    if (idx %% 10 == 1 || idx == length(rep_files)) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      if (idx > 1) {
        rate <- idx / elapsed
        eta <- (length(rep_files) - idx) / rate
        message(sprintf("  Replicate %d/%d (%.1f%%) - ETA: %.1f min",
                        idx, length(rep_files),
                        100 * idx / length(rep_files),
                        eta / 60))
      } else {
        message(sprintf("  Replicate %d/%d", idx, length(rep_files)))
      }
    }
    
    # Load replicate
    rep_data <- readr::read_csv(rep_file, show_col_types = FALSE)
    
    # Join with MEASYEAR
    rep_data <- rep_data |>
      dplyr::left_join(plot_data |> dplyr::select(CN, MEASYEAR),
                       by = "CN")
    
    # Map each plot to its NDVI window
    rep_data <- rep_data |>
      dplyr::rowwise() |>
      dplyr::mutate(
        window_name = map_year_to_window(MEASYEAR, cfg$level_window, ndvi_cfg$years_to_windows)
      ) |>
      dplyr::ungroup()
    
    # Convert to spatial using jittered coordinates
    valid_coords <- !is.na(rep_data$lat_jittered) & !is.na(rep_data$lon_jittered) &
      is.finite(rep_data$lat_jittered) & is.finite(rep_data$lon_jittered)
    
    if (sum(valid_coords) == 0) {
      warning("No valid jittered coordinates in ", basename(rep_file))
      next
    }
    
    rep_spatial <- rep_data[valid_coords, ]
    
    pts_sf <- sf::st_as_sf(rep_spatial,
                           coords = c("lon_jittered", "lat_jittered"),
                           crs = 4326, remove = FALSE)
    
    # Initialize NDVI columns
    rep_spatial$ndvi_modis <- NA_real_
    rep_spatial$ndvi_s2 <- NA_real_
    
    # Extract NDVI for each window/product combination
    for (window_name in unique(rep_spatial$window_name)) {
      if (is.na(window_name)) next
      
      window_idx <- rep_spatial$window_name == window_name & !is.na(rep_spatial$window_name)
      window_pts <- pts_sf[window_idx, ]
      
      # MODIS extraction
      modis_raster <- available_rasters |>
        dplyr::filter(product == "modis", window_name == !!window_name)
      
      if (nrow(modis_raster) == 1) {
        modis_vals <- extract_ndvi_values(
          window_pts,
          modis_raster$path[1],
          method = ndvi_cfg$modis_extraction_method,
          buffer_m = ndvi_cfg$modis_buffer_m %||% 0
        )
        rep_spatial$ndvi_modis[window_idx] <- modis_vals
      }
      
      # Sentinel-2 extraction
      s2_raster <- available_rasters |>
        dplyr::filter(product == "s2", window_name == !!window_name)
      
      if (nrow(s2_raster) == 1) {
        s2_vals <- extract_ndvi_values(
          window_pts,
          s2_raster$path[1],
          method = ndvi_cfg$s2_extraction_method,
          buffer_m = ndvi_cfg$s2_buffer_m %||% 0
        )
        rep_spatial$ndvi_s2[window_idx] <- s2_vals
      }
    }
    
    # Keep relevant columns
    hex_cols <- names(rep_spatial)[grepl("^hex_id_", names(rep_spatial))]
    
    result <- rep_spatial |>
      dplyr::select(CN, STATECD, MEASYEAR, replicate_id,
                    lon_jittered, lat_jittered,
                    window_name, ndvi_modis, ndvi_s2,
                    dplyr::all_of(hex_cols))
    
    all_results[[idx]] <- result
  }
  
  message("\n→ Combining results...")
  combined <- dplyr::bind_rows(all_results)
  
  message("  Total rows: ", format(nrow(combined), big.mark = ","))
  message("  MODIS values extracted: ", sum(!is.na(combined$ndvi_modis)))
  message("  S2 values extracted: ", sum(!is.na(combined$ndvi_s2)))
  
  combined
}

# ============================================================================
# NEFIN Extraction (true coordinates)
# ============================================================================

#' Extract NDVI for NEFIN plots using true coordinates
#' 
#' @param cfg Configuration list
#' @param available_rasters data.frame from find_ndvi_rasters()
#' @return data.frame with NDVI values
extract_nefin_ndvi <- function(cfg, available_rasters) {
  
  message("\n═══════════════════════════════════════════════════════════")
  message("EXTRACTING NDVI FOR NEFIN PLOTS (TRUE COORDINATES)")
  message("═══════════════════════════════════════════════════════════\n")
  
  ndvi_cfg <- cfg$ndvi
  nefin_file <- fs::path(cfg$project_dir, "data", "processed", "nefin_processed.csv")
  
  if (!fs::file_exists(nefin_file)) {
    message("⚠ NEFIN data not found: ", nefin_file)
    message("  Run: Rscript R/process_nefin_data.R")
    return(NULL)
  }
  
  # Load NEFIN data
  nefin <- readr::read_csv(nefin_file, show_col_types = FALSE)
  message("→ NEFIN plots: ", nrow(nefin))
  
  # Filter to existing rasters
  available_rasters <- available_rasters |> dplyr::filter(exists)
  
  if (nrow(available_rasters) == 0) {
    warning("No NDVI rasters found")
    return(NULL)
  }
  
  # Turn off s2
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  # Map each plot to its NDVI window
  nefin <- nefin |>
    dplyr::rowwise() |>
    dplyr::mutate(
      window_name = map_year_to_window(MEASYEAR, cfg$level_window, ndvi_cfg$years_to_windows)
    ) |>
    dplyr::ungroup()
  
  # Filter valid coordinates
  valid_coords <- !is.na(nefin$lat_public) & !is.na(nefin$lon_public) &
    is.finite(nefin$lat_public) & is.finite(nefin$lon_public)
  
  nefin <- nefin[valid_coords, ]
  
  # Convert to spatial
  pts_sf <- sf::st_as_sf(nefin,
                         coords = c("lon_public", "lat_public"),
                         crs = 4326, remove = FALSE)
  
  # Initialize NDVI columns
  nefin$ndvi_modis <- NA_real_
  nefin$ndvi_s2 <- NA_real_
  
  # Extract NDVI for each window/product
  for (window_name in unique(nefin$window_name)) {
    if (is.na(window_name)) next
    
    message("  Processing window: ", window_name)
    
    window_idx <- nefin$window_name == window_name & !is.na(nefin$window_name)
    window_pts <- pts_sf[window_idx, ]
    
    # MODIS extraction
    modis_raster <- available_rasters |>
      dplyr::filter(product == "modis", window_name == !!window_name)
    
    if (nrow(modis_raster) == 1) {
      modis_vals <- extract_ndvi_values(
        window_pts,
        modis_raster$path[1],
        method = ndvi_cfg$modis_extraction_method,
        buffer_m = ndvi_cfg$modis_buffer_m %||% 0
      )
      nefin$ndvi_modis[window_idx] <- modis_vals
    }
    
    # Sentinel-2 extraction
    s2_raster <- available_rasters |>
      dplyr::filter(product == "s2", window_name == !!window_name)
    
    if (nrow(s2_raster) == 1) {
      s2_vals <- extract_ndvi_values(
        window_pts,
        s2_raster$path[1],
        method = ndvi_cfg$s2_extraction_method,
        buffer_m = ndvi_cfg$s2_buffer_m %||% 0
      )
      nefin$ndvi_s2[window_idx] <- s2_vals
    }
  }
  
  message("\n  MODIS values extracted: ", sum(!is.na(nefin$ndvi_modis)))
  message("  S2 values extracted: ", sum(!is.na(nefin$ndvi_s2)))
  
  nefin
}

# ============================================================================
# Main Entry Point
# ============================================================================

#' Extract NDVI to all plot locations
#' 
#' @param project_dir Project root directory
#' @param overwrite Overwrite existing outputs
#' @param modis_only Only process MODIS
#' @param s2_only Only process Sentinel-2
#' @param max_reps Maximum FIA jitter replicates to process (NULL = all)
extract_ndvi_to_plots <- function(project_dir = ".",
                                  overwrite = FALSE,
                                  modis_only = FALSE,
                                  s2_only = FALSE,
                                  max_reps = NULL) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  NDVI Extraction to Plot Locations                       ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Load configuration
  cfg <- load_ndvi_config("configs/process.yml")
  cfg$project_dir <- project_dir
  
  # Adjust config based on flags
  if (modis_only) cfg$ndvi$use_s2 <- FALSE
  if (s2_only) cfg$ndvi$use_modis <- FALSE
  
  # Create output directory
  out_dir <- fs::path(project_dir, "data", "processed", "ndvi")
  fs::dir_create(out_dir, recurse = TRUE)
  
  fia_out <- fs::path(out_dir, "fia_ndvi_extracted.csv")
  nefin_out <- fs::path(out_dir, "nefin_ndvi_extracted.csv")
  manifest_out <- fs::path(out_dir, "extraction_manifest.yml")
  
  # Check if already processed
  if (!overwrite && fs::file_exists(fia_out) && fs::file_exists(manifest_out)) {
    message("✓ NDVI extraction already complete")
    message("  FIA: ", fia_out)
    if (fs::file_exists(nefin_out)) {
      message("  NEFIN: ", nefin_out)
    }
    message("  Use --overwrite to regenerate")
    return(invisible(list(fia = fia_out, nefin = nefin_out)))
  }
  
  # Find available NDVI rasters
  message("→ Discovering NDVI rasters...")
  available_rasters <- find_ndvi_rasters(cfg$ndvi)
  
  n_exist <- sum(available_rasters$exists)
  n_total <- nrow(available_rasters)
  
  message("  Found: ", n_exist, "/", n_total, " expected rasters")
  
  if (n_exist == 0) {
    message("\n⚠ No NDVI rasters found!")
    message("  Expected locations:")
    message("    MODIS: ", cfg$ndvi$modis_dir)
    message("    S2: ", cfg$ndvi$s2_dir)
    message("\n  Export rasters from Google Earth Engine first.")
    return(invisible(NULL))
  }
  
  # Show what's available
  message("\n→ Available rasters:")
  for (i in seq_len(nrow(available_rasters))) {
    r <- available_rasters[i, ]
    status <- if (r$exists) "✓" else "✗"
    message("  ", status, " ", r$product, " ", r$start_year, "-", r$end_year)
  }
  
  # Extract FIA NDVI
  fia_ndvi <- extract_fia_ndvi(cfg, available_rasters, max_reps = max_reps)
  
  if (!is.null(fia_ndvi) && nrow(fia_ndvi) > 0) {
    readr::write_csv(fia_ndvi, fia_out)
    message("\n✓ Wrote FIA NDVI: ", fia_out)
  }
  
  # Extract NEFIN NDVI
  nefin_ndvi <- extract_nefin_ndvi(cfg, available_rasters)
  
  if (!is.null(nefin_ndvi) && nrow(nefin_ndvi) > 0) {
    readr::write_csv(nefin_ndvi, nefin_out)
    message("\n✓ Wrote NEFIN NDVI: ", nefin_out)
  }
  
  # Write manifest
  manifest <- list(
    created = as.character(Sys.time()),
    config = cfg$ndvi,
    rasters_found = n_exist,
    rasters_expected = n_total,
    available_rasters = as.list(available_rasters),
    fia = list(
      output = as.character(fia_out),
      rows = if (!is.null(fia_ndvi)) nrow(fia_ndvi) else 0,
      modis_extracted = if (!is.null(fia_ndvi)) sum(!is.na(fia_ndvi$ndvi_modis)) else 0,
      s2_extracted = if (!is.null(fia_ndvi)) sum(!is.na(fia_ndvi$ndvi_s2)) else 0
    ),
    nefin = list(
      output = as.character(nefin_out),
      rows = if (!is.null(nefin_ndvi)) nrow(nefin_ndvi) else 0,
      modis_extracted = if (!is.null(nefin_ndvi)) sum(!is.na(nefin_ndvi$ndvi_modis)) else 0,
      s2_extracted = if (!is.null(nefin_ndvi)) sum(!is.na(nefin_ndvi$ndvi_s2)) else 0
    )
  )
  
  yaml::write_yaml(manifest, manifest_out)
  message("\n✓ Wrote manifest: ", manifest_out)
  
  # Summary
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("EXTRACTION SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("\n")
  
  if (!is.null(fia_ndvi)) {
    cat("FIA (jittered coordinates):\n")
    cat("  Total plot-replicates: ", format(nrow(fia_ndvi), big.mark = ","), "\n")
    cat("  MODIS values: ", format(sum(!is.na(fia_ndvi$ndvi_modis)), big.mark = ","), "\n")
    cat("  S2 values: ", format(sum(!is.na(fia_ndvi$ndvi_s2)), big.mark = ","), "\n")
  }
  
  cat("\n")
  
  if (!is.null(nefin_ndvi)) {
    cat("NEFIN (true coordinates):\n")
    cat("  Total plots: ", format(nrow(nefin_ndvi), big.mark = ","), "\n")
    cat("  MODIS values: ", format(sum(!is.na(nefin_ndvi$ndvi_modis)), big.mark = ","), "\n")
    cat("  S2 values: ", format(sum(!is.na(nefin_ndvi$ndvi_s2)), big.mark = ","), "\n")
  }
  
  cat("\n")
  cat("Output files:\n")
  cat("  ", fia_out, "\n")
  if (fs::file_exists(nefin_out)) {
    cat("  ", nefin_out, "\n")
  }
  cat("  ", manifest_out, "\n")
  
  cat("\n")
  cat("Next step:\n")
  cat("  Rscript R/build_ndvi_plot_dataset.R\n")
  cat("\n")
  
  invisible(list(
    fia = fia_out,
    nefin = nefin_out,
    manifest = manifest_out
  ))
}

# ============================================================================
# CLI Entry Point
# ============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  overwrite <- "--overwrite" %in% args
  modis_only <- "--modis-only" %in% args
  s2_only <- "--s2-only" %in% args
  
  # Parse max_reps argument
  max_reps <- NULL
  max_reps_arg <- grep("^--max-reps=", args, value = TRUE)
  if (length(max_reps_arg)) {
    max_reps <- as.integer(sub("^--max-reps=", "", max_reps_arg[1]))
  }
  
  extract_ndvi_to_plots(
    project_dir = ".",
    overwrite = overwrite,
    modis_only = modis_only,
    s2_only = s2_only,
    max_reps = max_reps
  )
}
