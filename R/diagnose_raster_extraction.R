#!/usr/bin/env Rscript
# ==============================================================================
# diagnose_raster_extraction.R
# ==============================================================================
# Purpose: Debug why extraction returns all NAs
# Checks: CRS, extent, resolution, plot coverage, actual values
# ==============================================================================

library(terra)
library(dplyr)
library(readr)

cat("\n")
cat(strrep("=", 70), "\n")
cat("  RASTER EXTRACTION DIAGNOSTICS\n")
cat(strrep("=", 70), "\n\n")

# ------------------------------------------------------------------------------
# 1. Check Raster Files
# ------------------------------------------------------------------------------

cat("1. CHECKING RASTER FILES\n")
cat(strrep("-", 50), "\n\n")

raster_files <- c(
  prism_tmean = "data/processed/prism/prism_tmean_ne_2015_2019.tif",
  prism_ppt = "data/processed/prism/prism_ppt_ne_2015_2019.tif",
  modis = "data/processed/ndvi/modis/MODIS_NDVI_5yr_blocked_2015_2019.tif",
  s2 = "data/processed/ndvi/s2/ndvi_2016_2019_raw.tif"
)

for (name in names(raster_files)) {
  path <- raster_files[name]
  cat("Raster:", name, "\n")
  cat("  Path:", path, "\n")
  
  if (!file.exists(path)) {
    cat("  STATUS: FILE NOT FOUND\n\n")
    next
  }
  
  r <- tryCatch(rast(path), error = function(e) NULL)
  if (is.null(r)) {
    cat("  STATUS: FAILED TO LOAD\n\n")
    next
  }
  
  cat("  CRS:", crs(r, describe = TRUE)$name, "\n")
  cat("  EPSG:", crs(r, describe = TRUE)$code, "\n")
  cat("  Extent:\n")
  e <- ext(r)
  cat("    xmin:", e$xmin, " xmax:", e$xmax, "\n")
  cat("    ymin:", e$ymin, " ymax:", e$ymax, "\n")
  cat("  Resolution:", res(r)[1], "x", res(r)[2], "\n")
  cat("  Dimensions:", nrow(r), "rows x", ncol(r), "cols\n")
  
  # Check for actual values
  vals <- values(r, mat = FALSE)
  n_valid <- sum(!is.na(vals))
  n_total <- length(vals)
  cat("  Values: ", n_valid, "/", n_total, " non-NA (", 
      round(100 * n_valid / n_total, 1), "%)\n", sep = "")
  
  if (n_valid > 0) {
    cat("  Range:", min(vals, na.rm = TRUE), "to", max(vals, na.rm = TRUE), "\n")
  }
  cat("\n")
}

# ------------------------------------------------------------------------------
# 2. Check Plot Coordinates
# ------------------------------------------------------------------------------

cat("\n2. CHECKING PLOT COORDINATES\n")
cat(strrep("-", 50), "\n\n")

# FIA jitter replicate
jitter_file <- "data/processed/mc_jitter_library/replicates/rep_0001.csv"
if (file.exists(jitter_file)) {
  jitter <- read_csv(jitter_file, show_col_types = FALSE)
  cat("FIA Jitter Rep 1:\n")
  cat("  Rows:", nrow(jitter), "\n")
  cat("  Columns:", paste(names(jitter)[1:min(10, ncol(jitter))], collapse = ", "), "...\n")
  
  if ("lon_jittered" %in% names(jitter) && "lat_jittered" %in% names(jitter)) {
    lon <- jitter$lon_jittered
    lat <- jitter$lat_jittered
    cat("  lon_jittered range:", min(lon, na.rm = TRUE), "to", max(lon, na.rm = TRUE), "\n")
    cat("  lat_jittered range:", min(lat, na.rm = TRUE), "to", max(lat, na.rm = TRUE), "\n")
    cat("  Valid coords:", sum(!is.na(lon) & !is.na(lat)), "/", nrow(jitter), "\n")
  } else {
    cat("  WARNING: lon_jittered/lat_jittered columns not found!\n")
    cat("  Available columns:", paste(names(jitter), collapse = ", "), "\n")
  }
} else {
  cat("FIA jitter file not found:", jitter_file, "\n")
}

cat("\n")

# NEFIN
nefin_file <- "data/processed/nefin_processed.csv"
if (file.exists(nefin_file)) {
  nefin <- read_csv(nefin_file, show_col_types = FALSE)
  cat("NEFIN:\n")
  cat("  Rows:", nrow(nefin), "\n")
  cat("  Columns:", paste(names(nefin)[1:min(10, ncol(nefin))], collapse = ", "), "...\n")
  
  # Find coordinate columns
  lon_col <- if ("lon_public" %in% names(nefin)) "lon_public" else if ("lon" %in% names(nefin)) "lon" else NULL
  lat_col <- if ("lat_public" %in% names(nefin)) "lat_public" else if ("lat" %in% names(nefin)) "lat" else NULL
  
  if (!is.null(lon_col) && !is.null(lat_col)) {
    lon <- nefin[[lon_col]]
    lat <- nefin[[lat_col]]
    cat("  Using:", lon_col, "/", lat_col, "\n")
    cat("  Longitude range:", min(lon, na.rm = TRUE), "to", max(lon, na.rm = TRUE), "\n")
    cat("  Latitude range:", min(lat, na.rm = TRUE), "to", max(lat, na.rm = TRUE), "\n")
    cat("  Valid coords:", sum(!is.na(lon) & !is.na(lat)), "/", nrow(nefin), "\n")
  } else {
    cat("  WARNING: Cannot find coordinate columns!\n")
    cat("  Available columns:", paste(names(nefin), collapse = ", "), "\n")
  }
} else {
  cat("NEFIN file not found:", nefin_file, "\n")
}

# ------------------------------------------------------------------------------
# 3. Test Extraction
# ------------------------------------------------------------------------------

cat("\n\n3. TEST EXTRACTION\n")
cat(strrep("-", 50), "\n\n")

# Use the first available raster
test_raster <- NULL
for (path in raster_files) {
  if (file.exists(path)) {
    test_raster <- path
    break
  }
}

if (is.null(test_raster)) {
  cat("ERROR: No rasters found to test!\n")
} else {
  cat("Testing with:", basename(test_raster), "\n\n")
  
  r <- rast(test_raster)
  raster_crs <- crs(r)
  raster_ext <- ext(r)
  
  # Get sample points from jitter
  if (file.exists(jitter_file)) {
    jitter <- read_csv(jitter_file, show_col_types = FALSE)
    
    # Take first 10 valid points
    valid <- !is.na(jitter$lon_jittered) & !is.na(jitter$lat_jittered)
    sample_pts <- jitter[valid, ][1:min(10, sum(valid)), ]
    
    cat("Sample points (first 10):\n")
    for (i in 1:nrow(sample_pts)) {
      cat(sprintf("  Point %d: lon=%.4f, lat=%.4f\n", 
                  i, sample_pts$lon_jittered[i], sample_pts$lat_jittered[i]))
    }
    
    # Create points in WGS84
    pts <- vect(
      data.frame(x = sample_pts$lon_jittered, y = sample_pts$lat_jittered),
      geom = c("x", "y"),
      crs = "EPSG:4326"
    )
    
    cat("\nPoints CRS: EPSG:4326 (WGS84)\n")
    cat("Raster CRS:", crs(r, describe = TRUE)$name, "\n")
    
    # Check if raster is also in 4326
    raster_epsg <- crs(r, describe = TRUE)$code
    
    if (!is.na(raster_epsg) && raster_epsg == "4326") {
      cat("\nBoth in EPSG:4326, no projection needed.\n")
      pts_proj <- pts
    } else {
      cat("\nProjecting points to raster CRS...\n")
      pts_proj <- tryCatch({
        project(pts, raster_crs)
      }, error = function(e) {
        cat("  Projection error:", e$message, "\n")
        return(NULL)
      })
    }
    
    if (!is.null(pts_proj)) {
      # Check if points fall within raster extent
      pts_coords <- crds(pts_proj)
      cat("\nProjected point coordinates:\n")
      for (i in 1:nrow(pts_coords)) {
        in_extent <- pts_coords[i, 1] >= raster_ext$xmin & 
                     pts_coords[i, 1] <= raster_ext$xmax &
                     pts_coords[i, 2] >= raster_ext$ymin & 
                     pts_coords[i, 2] <= raster_ext$ymax
        status <- if (in_extent) "IN EXTENT" else "OUTSIDE"
        cat(sprintf("  Point %d: x=%.2f, y=%.2f [%s]\n", 
                    i, pts_coords[i, 1], pts_coords[i, 2], status))
      }
      
      # Extract
      cat("\nExtracting values...\n")
      vals <- extract(r, pts_proj)
      cat("Extracted values:\n")
      print(vals)
      
      n_valid <- sum(!is.na(vals[[2]]))
      cat("\nResult:", n_valid, "/", nrow(vals), "valid values extracted\n")
      
      if (n_valid == 0) {
        cat("\n*** ALL VALUES ARE NA ***\n")
        cat("\nPossible causes:\n")
        cat("  1. Points fall outside raster extent\n")
        cat("  2. Points fall in NoData areas of raster\n")
        cat("  3. CRS mismatch causing incorrect projection\n")
        cat("  4. Raster has no valid data in study area\n")
      }
    }
  }
}

# ------------------------------------------------------------------------------
# 4. Visual Check - Extent Overlap
# ------------------------------------------------------------------------------

cat("\n\n4. EXTENT COMPARISON\n")
cat(strrep("-", 50), "\n\n")

# Get plot extent
if (file.exists(jitter_file)) {
  jitter <- read_csv(jitter_file, show_col_types = FALSE)
  valid <- !is.na(jitter$lon_jittered) & !is.na(jitter$lat_jittered)
  
  plot_ext <- list(
    xmin = min(jitter$lon_jittered[valid]),
    xmax = max(jitter$lon_jittered[valid]),
    ymin = min(jitter$lat_jittered[valid]),
    ymax = max(jitter$lat_jittered[valid])
  )
  
  cat("Plot extent (WGS84):\n")
  cat(sprintf("  Longitude: %.4f to %.4f\n", plot_ext$xmin, plot_ext$xmax))
  cat(sprintf("  Latitude:  %.4f to %.4f\n", plot_ext$ymin, plot_ext$ymax))
  
  # Compare with each raster
  cat("\nRaster extents:\n")
  for (name in names(raster_files)) {
    path <- raster_files[name]
    if (!file.exists(path)) next
    
    r <- rast(path)
    e <- ext(r)
    epsg <- crs(r, describe = TRUE)$code
    
    cat(sprintf("\n  %s (EPSG:%s):\n", name, ifelse(is.na(epsg), "unknown", epsg)))
    cat(sprintf("    X: %.4f to %.4f\n", e$xmin, e$xmax))
    cat(sprintf("    Y: %.4f to %.4f\n", e$ymin, e$ymax))
    
    # If raster is in 4326, check overlap directly
    if (!is.na(epsg) && epsg == "4326") {
      overlap_x <- plot_ext$xmin <= e$xmax && plot_ext$xmax >= e$xmin
      overlap_y <- plot_ext$ymin <= e$ymax && plot_ext$ymax >= e$ymin
      
      if (overlap_x && overlap_y) {
        cat("    OVERLAP: YES - plots should intersect raster\n")
      } else {
        cat("    OVERLAP: NO - plots are OUTSIDE raster extent!\n")
        if (!overlap_x) cat("      - No X overlap\n")
        if (!overlap_y) cat("      - No Y overlap\n")
      }
    }
  }
}

cat("\n")
cat(strrep("=", 70), "\n")
cat("  DIAGNOSTICS COMPLETE\n")
cat(strrep("=", 70), "\n\n")
