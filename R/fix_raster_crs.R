# ==============================================================================
# fix_raster_crs.R
# ==============================================================================
# Purpose: Assign missing CRS to MODIS and S2 rasters exported from GEE
#
# The coordinate ranges suggest EPSG:5070 (NAD83 / Conus Albers):
#   MODIS extent: 1,319,500 to 2,264,000 (x), 2,149,000 to 3,013,250 (y)
#   These are typical Albers meter coordinates for the Northeast US
#
# ==============================================================================

library(terra)

cat("\n")
cat(strrep("=", 60), "\n")
cat("  FIX RASTER CRS\n")
cat(strrep("=", 60), "\n\n")

# Define the CRS to assign (NAD83 / Conus Albers)
# This is the standard projection for USGS/EPA continental US data
TARGET_CRS <- "EPSG:5070"

# Rasters that need CRS assignment
rasters_to_fix <- list(
  # MODIS NDVI
  list(
    path = "data/processed/ndvi/modis/MODIS_NDVI_5yr_blocked_2000_2004.tif",
    crs = TARGET_CRS
  ),
  list(
    path = "data/processed/ndvi/modis/MODIS_NDVI_5yr_blocked_2005_2009.tif",
    crs = TARGET_CRS
  ),
  list(
    path = "data/processed/ndvi/modis/MODIS_NDVI_5yr_blocked_2010_2014.tif",
    crs = TARGET_CRS
  ),
  list(
    path = "data/processed/ndvi/modis/MODIS_NDVI_5yr_blocked_2015_2019.tif",
    crs = TARGET_CRS
  ),
  list(
    path = "data/processed/ndvi/modis/MODIS_NDVI_5yr_blocked_2020_2024.tif",
    crs = TARGET_CRS
  ),
  # Sentinel-2 NDVI
  list(
    path = "data/processed/ndvi/s2/ndvi_2016_2019_raw.tif",
    crs = TARGET_CRS
  ),
  list(
    path = "data/processed/ndvi/s2/ndvi_2020_2024_raw.tif",
    crs = TARGET_CRS
  )
)

# Process each raster
for (item in rasters_to_fix) {
  path <- item$path
  
  if (!file.exists(path)) {
    cat("SKIP (not found):", path, "\n")
    next
  }
  
  cat("Processing:", basename(path), "\n")
  
  # Load raster
  r <- rast(path)
  
  # Check current CRS
  current_crs <- crs(r, describe = TRUE)$name
  cat("  Current CRS:", ifelse(is.na(current_crs), "NONE", current_crs), "\n")
  
  # Check extent to verify it's in projected coords
  e <- ext(r)
  if (e$xmin > 0 && e$xmin < 180 && e$ymin > 0 && e$ymin < 90) {
    cat("  Extent looks like geographic (lat/lon), skipping CRS assignment\n")
    next
  }
  
  # Assign CRS
  crs(r) <- item$crs
  cat("  Assigned CRS:", item$crs, "\n")
  
  # Create backup
  backup_path <- paste0(tools::file_path_sans_ext(path), "_backup.tif")
  if (!file.exists(backup_path)) {
    file.copy(path, backup_path)
    cat("  Backup created:", basename(backup_path), "\n")
  }
  
  # Save with CRS
  # Use writeRaster with overwrite
  writeRaster(r, path, overwrite = TRUE)
  cat("  Saved with CRS\n\n")
}

# ------------------------------------------------------------------------------
# Verify the fix
# ------------------------------------------------------------------------------

cat("\n")
cat(strrep("=", 60), "\n")
cat("  VERIFICATION\n")
cat(strrep("=", 60), "\n\n")

# Check one MODIS file
test_path <- "data/processed/ndvi/modis/MODIS_NDVI_5yr_blocked_2015_2019.tif"
if (file.exists(test_path)) {
  r <- rast(test_path)
  cat("MODIS 2015-2019:\n")
  cat("  CRS:", crs(r, describe = TRUE)$name, "\n")
  cat("  EPSG:", crs(r, describe = TRUE)$code, "\n")
  
  # Test extraction with a sample point
  # Burlington, VT in WGS84: -73.2121, 44.4759
  test_lon <- -73.2121
  test_lat <- 44.4759
  
  cat("\nTest extraction (Burlington, VT):\n")
  cat("  Input coords (WGS84):", test_lon, ",", test_lat, "\n")
  
  # Create point and project
  pt <- vect(data.frame(x = test_lon, y = test_lat), 
             geom = c("x", "y"), crs = "EPSG:4326")
  pt_proj <- project(pt, crs(r))
  
  proj_coords <- crds(pt_proj)
  cat("  Projected coords:", proj_coords[1,1], ",", proj_coords[1,2], "\n")
  
  # Extract
  val <- extract(r, pt_proj)
  cat("  Extracted NDVI:", val[[2]], "\n")
  
  if (!is.na(val[[2]])) {
    cat("\n  ✓ SUCCESS! Extraction working!\n")
  } else {
    cat("\n  ✗ Still NA - point may be outside raster or in NoData area\n")
  }
}

cat("\n")
cat(strrep("=", 60), "\n")
cat("  DONE\n")
cat(strrep("=", 60), "\n\n")
