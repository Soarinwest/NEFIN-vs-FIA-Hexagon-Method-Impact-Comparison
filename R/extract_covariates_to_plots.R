#!/usr/bin/env Rscript
# R/extract_covariates_to_plots.R
# Extract climate (PRISM) and NDVI covariates to plot locations
# Works with existing manually downloaded TIFFs

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(readr)
  library(fs)
  library(tidyr)
})

extract_covariates_to_plots <- function(
    plot_assignments_path = "data/processed/plot_hex_assignments.csv",
    prism_dir = "data/processed/prism",
    ndvi_dir = "data/processed/ndvi",
    output_dir = "data/processed",
    use_jittered = FALSE,
    jitter_rep = 1,
    overwrite = FALSE
) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Extract Covariates to Plot Locations                    ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # =========================================================================
  # Load plot data
  # =========================================================================
  
  if (use_jittered) {
    # Use jittered coordinates from MC library
    jitter_file <- sprintf("data/processed/mc_jitter_library/replicates/rep_%04d.csv", jitter_rep)
    if (!fs::file_exists(jitter_file)) {
      stop("Jitter replicate file not found: ", jitter_file)
    }
    cat("→ Loading jittered coordinates (replicate", jitter_rep, ")...\n")
    plots <- readr::read_csv(jitter_file, show_col_types = FALSE)
    lon_col <- "lon_jittered"
    lat_col <- "lat_jittered"
  } else {
    # Use original (fuzzed) FIA coordinates
    cat("→ Loading plot assignments...\n")
    plots <- readr::read_csv(plot_assignments_path, show_col_types = FALSE)
    lon_col <- "LON"
    lat_col <- "LAT"
  }
  
  cat("  Loaded", nrow(plots), "plots\n")
  
  # Check for coordinate columns
  if (!lon_col %in% names(plots) || !lat_col %in% names(plots)) {
    # Try alternative names
    if ("lon" %in% names(plots)) lon_col <- "lon"
    if ("lat" %in% names(plots)) lat_col <- "lat"
    if ("longitude" %in% names(plots)) lon_col <- "longitude"
    if ("latitude" %in% names(plots)) lat_col <- "latitude"
  }
  
  if (!lon_col %in% names(plots) || !lat_col %in% names(plots)) {
    stop("Cannot find longitude/latitude columns. Available: ", paste(names(plots), collapse = ", "))
  }
  
  cat("  Using coordinates:", lon_col, "/", lat_col, "\n")
  
  # Filter out rows with missing coordinates
  plots_valid <- plots %>%
    filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]]))
  
  cat("  Valid coordinates:", nrow(plots_valid), "of", nrow(plots), "\n")
  
  # Convert to spatial points (WGS84)
  cat("\n→ Converting to spatial points...\n")
  plots_sf <- st_as_sf(plots_valid, 
                       coords = c(lon_col, lat_col), 
                       crs = 4326)
  
  # Get bounding box for diagnostics
  bbox <- st_bbox(plots_sf)
  cat("  Bounding box (WGS84):\n")
  cat("    Lon:", round(bbox["xmin"], 2), "to", round(bbox["xmax"], 2), "\n")
  cat("    Lat:", round(bbox["ymin"], 2), "to", round(bbox["ymax"], 2), "\n")
  
  # =========================================================================
  # Extract PRISM climate data
  # =========================================================================
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("Extracting PRISM Climate Data\n")
  cat("═══════════════════════════════════════════════════════════\n")
  
  prism_files <- list.files(prism_dir, pattern = "\\.tif$", full.names = TRUE)
  
  if (length(prism_files) == 0) {
    cat("  ⚠ No PRISM TIF files found in", prism_dir, "\n")
    plots_valid$tmean <- NA_real_
    plots_valid$ppt <- NA_real_
  } else {
    cat("  Found", length(prism_files), "PRISM files:\n")
    for (f in prism_files) cat("    -", basename(f), "\n")
    
    # Determine which time period to use based on MEASYEAR
    # Map measurement years to PRISM windows
    plots_valid <- plots_valid %>%
      mutate(
        prism_window = case_when(
          MEASYEAR >= 2020 ~ "2020_2024",
          MEASYEAR >= 2015 ~ "2015_2019",
          MEASYEAR >= 2010 ~ "2010_2014",
          MEASYEAR >= 2005 ~ "2005_2009",
          TRUE ~ "2000_2004"
        )
      )
    
    # Initialize columns
    plots_valid$tmean <- NA_real_
    plots_valid$ppt <- NA_real_
    
    # Process each time window
    for (window in unique(plots_valid$prism_window)) {
      cat("\n  Processing window:", window, "\n")
      
      idx <- which(plots_valid$prism_window == window)
      if (length(idx) == 0) next
      
      cat("    Plots in window:", length(idx), "\n")
      
      # Find matching files
      tmean_file <- prism_files[grepl(paste0("tmean.*", window), prism_files)]
      ppt_file <- prism_files[grepl(paste0("ppt.*", window), prism_files)]
      
      # Extract tmean
      if (length(tmean_file) == 1) {
        cat("    Loading:", basename(tmean_file), "\n")
        r_tmean <- terra::rast(tmean_file)
        
        # Check CRS and reproject points if needed
        pts_window <- plots_sf[idx, ]
        if (!identical(st_crs(pts_window)$wkt, crs(r_tmean))) {
          cat("    Reprojecting points to raster CRS...\n")
          pts_window <- st_transform(pts_window, crs(r_tmean))
        }
        
        # Extract values
        vals <- terra::extract(r_tmean, terra::vect(pts_window))
        
        # Get the value column (might be named differently)
        val_col <- names(vals)[2]  # First col is ID
        if (!is.null(val_col)) {
          plots_valid$tmean[idx] <- vals[[val_col]]
          n_valid <- sum(!is.na(vals[[val_col]]))
          cat("    Extracted tmean:", n_valid, "valid values\n")
        }
      } else {
        cat("    ⚠ tmean file not found for window", window, "\n")
      }
      
      # Extract ppt
      if (length(ppt_file) == 1) {
        cat("    Loading:", basename(ppt_file), "\n")
        r_ppt <- terra::rast(ppt_file)
        
        pts_window <- plots_sf[idx, ]
        if (!identical(st_crs(pts_window)$wkt, crs(r_ppt))) {
          pts_window <- st_transform(pts_window, crs(r_ppt))
        }
        
        vals <- terra::extract(r_ppt, terra::vect(pts_window))
        val_col <- names(vals)[2]
        if (!is.null(val_col)) {
          plots_valid$ppt[idx] <- vals[[val_col]]
          n_valid <- sum(!is.na(vals[[val_col]]))
          cat("    Extracted ppt:", n_valid, "valid values\n")
        }
      } else {
        cat("    ⚠ ppt file not found for window", window, "\n")
      }
    }
    
    cat("\n  Climate extraction summary:\n")
    cat("    tmean: ", sum(!is.na(plots_valid$tmean)), "/", nrow(plots_valid), "valid\n")
    cat("    ppt:   ", sum(!is.na(plots_valid$ppt)), "/", nrow(plots_valid), "valid\n")
  }
  
  # =========================================================================
  # Extract NDVI data
  # =========================================================================
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("Extracting NDVI Data\n")
  cat("═══════════════════════════════════════════════════════════\n")
  
  # MODIS NDVI
  modis_dir <- fs::path(ndvi_dir, "modis")
  modis_files <- if (fs::dir_exists(modis_dir)) {
    list.files(modis_dir, pattern = "\\.tif$", full.names = TRUE)
  } else {
    list.files(ndvi_dir, pattern = "MODIS.*\\.tif$", full.names = TRUE)
  }
  
  plots_valid$ndvi_modis <- NA_real_
  
  if (length(modis_files) > 0) {
    cat("  Found", length(modis_files), "MODIS NDVI files:\n")
    for (f in modis_files) cat("    -", basename(f), "\n")
    
    # Map years to MODIS windows
    plots_valid <- plots_valid %>%
      mutate(
        modis_window = case_when(
          MEASYEAR >= 2020 ~ "2020_2024",
          MEASYEAR >= 2015 ~ "2015_2019",
          MEASYEAR >= 2010 ~ "2010_2014",
          MEASYEAR >= 2005 ~ "2005_2009",
          TRUE ~ "2000_2004"
        )
      )
    
    for (window in unique(plots_valid$modis_window)) {
      cat("\n  Processing MODIS window:", window, "\n")
      
      idx <- which(plots_valid$modis_window == window)
      if (length(idx) == 0) next
      
      # Find matching file
      modis_file <- modis_files[grepl(window, modis_files)]
      
      if (length(modis_file) == 1) {
        cat("    Loading:", basename(modis_file), "\n")
        r_ndvi <- terra::rast(modis_file)
        
        pts_window <- plots_sf[idx, ]
        if (!identical(st_crs(pts_window)$wkt, crs(r_ndvi))) {
          cat("    Reprojecting points to raster CRS...\n")
          pts_window <- st_transform(pts_window, crs(r_ndvi))
        }
        
        vals <- terra::extract(r_ndvi, terra::vect(pts_window))
        val_col <- names(vals)[2]
        if (!is.null(val_col)) {
          plots_valid$ndvi_modis[idx] <- vals[[val_col]]
          n_valid <- sum(!is.na(vals[[val_col]]))
          cat("    Extracted:", n_valid, "valid values\n")
        }
      } else {
        cat("    ⚠ MODIS file not found for window", window, "\n")
      }
    }
    
    cat("\n  MODIS NDVI: ", sum(!is.na(plots_valid$ndvi_modis)), "/", nrow(plots_valid), "valid\n")
  } else {
    cat("  ⚠ No MODIS NDVI files found\n")
  }
  
  # Sentinel-2 NDVI
  s2_dir <- fs::path(ndvi_dir, "s2")
  s2_files <- if (fs::dir_exists(s2_dir)) {
    list.files(s2_dir, pattern = "\\.tif$", full.names = TRUE)
  } else {
    list.files(ndvi_dir, pattern = "s2.*\\.tif$|sentinel.*\\.tif$", full.names = TRUE, ignore.case = TRUE)
  }
  
  plots_valid$ndvi_s2 <- NA_real_
  
  if (length(s2_files) > 0) {
    cat("\n  Found", length(s2_files), "Sentinel-2 NDVI files:\n")
    for (f in s2_files) cat("    -", basename(f), "\n")
    
    # For S2, just use the first/only file (typically covers 2016-2019)
    s2_file <- s2_files[1]
    cat("    Loading:", basename(s2_file), "\n")
    r_s2 <- terra::rast(s2_file)
    
    # Only extract for plots in relevant years (2016-2019)
    idx_s2 <- which(plots_valid$MEASYEAR >= 2016 & plots_valid$MEASYEAR <= 2019)
    
    if (length(idx_s2) > 0) {
      pts_s2 <- plots_sf[idx_s2, ]
      if (!identical(st_crs(pts_s2)$wkt, crs(r_s2))) {
        cat("    Reprojecting points to raster CRS...\n")
        pts_s2 <- st_transform(pts_s2, crs(r_s2))
      }
      
      vals <- terra::extract(r_s2, terra::vect(pts_s2))
      val_col <- names(vals)[2]
      if (!is.null(val_col)) {
        plots_valid$ndvi_s2[idx_s2] <- vals[[val_col]]
        n_valid <- sum(!is.na(vals[[val_col]]))
        cat("    Extracted:", n_valid, "valid values for 2016-2019 plots\n")
      }
    }
    
    cat("\n  Sentinel-2 NDVI: ", sum(!is.na(plots_valid$ndvi_s2)), "/", nrow(plots_valid), "valid\n")
  } else {
    cat("  ⚠ No Sentinel-2 NDVI files found\n")
  }
  
  # =========================================================================
  # Save results
  # =========================================================================
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("Saving Results\n")
  cat("═══════════════════════════════════════════════════════════\n")
  
  # Create output directories
  climate_dir <- fs::path(output_dir, "climate_at_plots")
  ndvi_out_dir <- fs::path(output_dir, "ndvi_at_plots")
  fs::dir_create(climate_dir, recurse = TRUE)
  fs::dir_create(ndvi_out_dir, recurse = TRUE)
  
  # Select output columns
  id_cols <- c("CN", "STATECD", "MEASYEAR")
  if (use_jittered) {
    id_cols <- c(id_cols, "replicate_id", "lon_jittered", "lat_jittered")
  }
  
  # Hex columns if present
  hex_cols <- names(plots_valid)[grepl("^hex_id_", names(plots_valid))]
  
  # Climate output
  climate_out <- plots_valid %>%
    select(all_of(c(id_cols, hex_cols, "tmean", "ppt")))
  
  climate_file <- fs::path(climate_dir, "fia_climate.csv")
  readr::write_csv(climate_out, climate_file)
  cat("  ✓ Wrote:", climate_file, "\n")
  
  # NDVI output
  ndvi_out <- plots_valid %>%
    select(all_of(c(id_cols, hex_cols, "ndvi_modis", "ndvi_s2")))
  
  ndvi_file <- fs::path(ndvi_out_dir, "fia_ndvi.csv")
  readr::write_csv(ndvi_out, ndvi_file)
  cat("  ✓ Wrote:", ndvi_file, "\n")
  
  # Combined output
  combined_out <- plots_valid %>%
    select(all_of(c(id_cols, hex_cols, "tmean", "ppt", "ndvi_modis", "ndvi_s2")))
  
  combined_file <- fs::path(output_dir, "fia_covariates.csv")
  readr::write_csv(combined_out, combined_file)
  cat("  ✓ Wrote:", combined_file, "\n")
  
  # =========================================================================
  # Summary
  # =========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  EXTRACTION COMPLETE                                      ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  cat("Summary:\n")
  cat("  Total plots:", nrow(plots_valid), "\n")
  cat("  tmean:      ", sum(!is.na(plots_valid$tmean)), "valid (", 
      round(100 * sum(!is.na(plots_valid$tmean)) / nrow(plots_valid), 1), "%)\n")
  cat("  ppt:        ", sum(!is.na(plots_valid$ppt)), "valid (",
      round(100 * sum(!is.na(plots_valid$ppt)) / nrow(plots_valid), 1), "%)\n")
  cat("  ndvi_modis: ", sum(!is.na(plots_valid$ndvi_modis)), "valid (",
      round(100 * sum(!is.na(plots_valid$ndvi_modis)) / nrow(plots_valid), 1), "%)\n")
  cat("  ndvi_s2:    ", sum(!is.na(plots_valid$ndvi_s2)), "valid (",
      round(100 * sum(!is.na(plots_valid$ndvi_s2)) / nrow(plots_valid), 1), "%)\n")
  
  cat("\nOutput files:\n")
  cat("  ", climate_file, "\n")
  cat("  ", ndvi_file, "\n")
  cat("  ", combined_file, "\n")
  
  invisible(combined_out)
}

# =============================================================================
# Also extract for NEFIN plots
# =============================================================================

extract_covariates_to_nefin <- function(
    nefin_path = "data/processed/nefin_processed.csv",
    prism_dir = "data/processed/prism",
    ndvi_dir = "data/processed/ndvi",
    output_dir = "data/processed"
) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("Extracting Covariates for NEFIN Plots\n")
  cat("═══════════════════════════════════════════════════════════\n")
  
  if (!fs::file_exists(nefin_path)) {
    cat("  ⚠ NEFIN data not found:", nefin_path, "\n")
    return(invisible(NULL))
  }
  
  nefin <- readr::read_csv(nefin_path, show_col_types = FALSE)
  cat("  Loaded", nrow(nefin), "NEFIN plots\n")
  
  # Find coordinate columns
  lon_col <- if ("LON" %in% names(nefin)) "LON" else if ("lon" %in% names(nefin)) "lon" else "longitude"
  lat_col <- if ("LAT" %in% names(nefin)) "LAT" else if ("lat" %in% names(nefin)) "lat" else "latitude"
  
  # Convert to sf
  nefin_sf <- st_as_sf(nefin, coords = c(lon_col, lat_col), crs = 4326)
  
  # Initialize covariate columns
  nefin$tmean <- NA_real_
  nefin$ppt <- NA_real_
  nefin$ndvi_modis <- NA_real_
  nefin$ndvi_s2 <- NA_real_
  
  # Extract PRISM - use 2015-2019 as default for NEFIN
  prism_files <- list.files(prism_dir, pattern = "\\.tif$", full.names = TRUE)
  
  tmean_file <- prism_files[grepl("tmean.*2015_2019", prism_files)]
  if (length(tmean_file) == 1) {
    r <- terra::rast(tmean_file)
    pts <- st_transform(nefin_sf, crs(r))
    vals <- terra::extract(r, terra::vect(pts))
    nefin$tmean <- vals[[2]]
    cat("  tmean extracted:", sum(!is.na(nefin$tmean)), "valid\n")
  }
  
  ppt_file <- prism_files[grepl("ppt.*2015_2019", prism_files)]
  if (length(ppt_file) == 1) {
    r <- terra::rast(ppt_file)
    pts <- st_transform(nefin_sf, crs(r))
    vals <- terra::extract(r, terra::vect(pts))
    nefin$ppt <- vals[[2]]
    cat("  ppt extracted:", sum(!is.na(nefin$ppt)), "valid\n")
  }
  
  # Extract MODIS NDVI
  modis_dir <- fs::path(ndvi_dir, "modis")
  modis_file <- list.files(modis_dir, pattern = "2015_2019.*\\.tif$", full.names = TRUE)
  if (length(modis_file) == 1) {
    r <- terra::rast(modis_file)
    pts <- st_transform(nefin_sf, crs(r))
    vals <- terra::extract(r, terra::vect(pts))
    nefin$ndvi_modis <- vals[[2]]
    cat("  ndvi_modis extracted:", sum(!is.na(nefin$ndvi_modis)), "valid\n")
  }
  
  # Save
  climate_file <- fs::path(output_dir, "climate_at_plots", "nefin_climate.csv")
  ndvi_file <- fs::path(output_dir, "ndvi_at_plots", "nefin_ndvi.csv")
  
  readr::write_csv(nefin %>% select(any_of(c("PLOT_ID", "LON", "LAT", "tmean", "ppt"))), climate_file)
  readr::write_csv(nefin %>% select(any_of(c("PLOT_ID", "LON", "LAT", "ndvi_modis", "ndvi_s2"))), ndvi_file)
  
  cat("  ✓ Wrote:", climate_file, "\n")
  cat("  ✓ Wrote:", ndvi_file, "\n")
  
  invisible(nefin)
}

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  use_jittered <- "--jittered" %in% args
  include_nefin <- "--nefin" %in% args || "--all" %in% args
  
  # Extract FIA
  extract_covariates_to_plots(use_jittered = use_jittered)
  
  # Extract NEFIN
  if (include_nefin) {
    extract_covariates_to_nefin()
  }
}
