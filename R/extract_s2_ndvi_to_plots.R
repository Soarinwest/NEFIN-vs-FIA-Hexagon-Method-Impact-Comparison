#!/usr/bin/env Rscript
# R/extract_s2_ndvi_to_plots.R
# Re-extract Sentinel-2 NDVI to NEFIN and FIA plots
# Ensures temporal alignment: 2020-2024 plots get 2020-2024 S2 NDVI
# Author: Soren Donisvitch
# Date: December 2024

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(sf)
  library(terra)
  library(fs)
})

extract_s2_to_plots <- function(
    nefin_path = "data/processed/nefin_processed.csv",
    fia_cov_path = "data/processed/fia_covariates.csv",
    s2_dir = "data/processed/ndvi/s2",
    output_dir = "data/processed/ndvi_at_plots"
) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  Extract Sentinel-2 NDVI to Plot Locations                                ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  fs::dir_create(output_dir, recurse = TRUE)
  
  # =========================================================================
  # 1. CHECK AVAILABLE S2 RASTERS
  # =========================================================================
  
  cat("═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 1: Checking Available S2 Rasters\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  s2_files <- list.files(s2_dir, pattern = "\\.tif$", full.names = TRUE)
  
  cat("  S2 directory:", s2_dir, "\n")
  cat("  Files found:\n")
  for (f in s2_files) {
    r <- try(rast(f), silent = TRUE)
    if (!inherits(r, "try-error")) {
      cat("    •", basename(f), "- res:", round(res(r)[1], 1), "m\n")
    } else {
      cat("    •", basename(f), "- (could not read)\n")
    }
  }
  
  # Find the 2020-2024 raster
  s2_2020_files <- s2_files[grepl("2020|2025", s2_files)]
  
  if (length(s2_2020_files) == 0) {
    cat("\n  ⚠ No 2020-2024 S2 raster found!\n")
    cat("  Looking for any S2 raster to use...\n")
    if (length(s2_files) > 0) {
      s2_path <- s2_files[1]
      cat("  Using:", basename(s2_path), "\n")
    } else {
      stop("No S2 rasters found in ", s2_dir)
    }
  } else {
    s2_path <- s2_2020_files[1]
    cat("\n  Using 2020-2024 raster:", basename(s2_path), "\n")
  }
  
  # Also check for 2016-2019 if doing temporal comparison
  s2_2016_files <- s2_files[grepl("2016|2019", s2_files) & !grepl("2020", s2_files)]
  s2_path_2016 <- if (length(s2_2016_files) > 0) s2_2016_files[1] else NULL
  
  if (!is.null(s2_path_2016)) {
    cat("  Also found 2016-2019 raster:", basename(s2_path_2016), "\n")
  }
  
  # Load primary raster
  r_s2 <- rast(s2_path)
  cat("\n  S2 raster info:\n")
  cat("    Resolution:", round(res(r_s2)[1], 1), "m\n")
  cat("    CRS:", crs(r_s2, describe = TRUE)$name, "\n")
  cat("    Extent:", round(ext(r_s2)[1:4], 0), "\n")
  
  # =========================================================================
  # 2. LOAD PLOT DATA
  # =========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 2: Loading Plot Data\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # NEFIN
  nefin <- read_csv(nefin_path, show_col_types = FALSE)
  lon_col <- names(nefin)[grepl("lon", names(nefin), ignore.case = TRUE)][1]
  lat_col <- names(nefin)[grepl("lat", names(nefin), ignore.case = TRUE)][1]
  
  cat("  NEFIN plots:", nrow(nefin), "\n")
  cat("    Coordinate columns:", lon_col, ",", lat_col, "\n")
  
  # FIA (from covariates file which has original coords)
  fia <- read_csv(fia_cov_path, show_col_types = FALSE)
  cat("  FIA plots:", nrow(fia), "\n")
  
  # Check for coordinate columns in FIA
  fia_lon_col <- names(fia)[grepl("lon", names(fia), ignore.case = TRUE)][1]
  fia_lat_col <- names(fia)[grepl("lat", names(fia), ignore.case = TRUE)][1]
  
  if (is.na(fia_lon_col)) {
    cat("  ⚠ FIA covariates missing coordinates - loading from plot_hex_assignments\n")
    fia_plots <- read_csv("data/processed/plot_hex_assignments.csv", show_col_types = FALSE)
    fia <- fia %>%
      left_join(fia_plots %>% select(CN, lon_original, lat_original), by = "CN")
    fia_lon_col <- "lon_original"
    fia_lat_col <- "lat_original"
  }
  
  cat("    FIA coordinate columns:", fia_lon_col, ",", fia_lat_col, "\n")
  
  # =========================================================================
  # 3. EXTRACT S2 NDVI TO NEFIN
  # =========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 3: Extracting S2 NDVI to NEFIN Plots\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Convert to sf
  nefin_sf <- st_as_sf(nefin, coords = c(lon_col, lat_col), crs = 4326)
  nefin_reproj <- st_transform(nefin_sf, crs(r_s2))
  
  cat("  Extracting...\n")
  s2_vals <- terra::extract(r_s2, vect(nefin_reproj))
  
  nefin$ndvi_s2 <- s2_vals[[2]]
  
  # Summary
  valid_s2 <- sum(!is.na(nefin$ndvi_s2))
  cat("  Valid extractions:", valid_s2, "/", nrow(nefin), 
      "(", round(100 * valid_s2 / nrow(nefin), 1), "%)\n")
  cat("  S2 NDVI range:", round(min(nefin$ndvi_s2, na.rm = TRUE), 3), "to",
      round(max(nefin$ndvi_s2, na.rm = TRUE), 3), "\n")
  cat("  S2 NDVI mean:", round(mean(nefin$ndvi_s2, na.rm = TRUE), 3), "\n")
  
  # Save NEFIN with S2
  nefin_out <- nefin %>% select(CN, ndvi_s2)
  
  # Check if there's an existing nefin_ndvi file to update
  nefin_ndvi_path <- fs::path(output_dir, "nefin_ndvi.csv")
  if (fs::file_exists(nefin_ndvi_path)) {
    existing <- read_csv(nefin_ndvi_path, show_col_types = FALSE)
    cat("\n  Updating existing nefin_ndvi.csv...\n")
    
    # Remove old ndvi_s2 if exists
    if ("ndvi_s2" %in% names(existing)) {
      existing <- existing %>% select(-ndvi_s2)
    }
    
    # Join new S2 values
    existing <- existing %>%
      left_join(nefin_out, by = "CN")
    
    write_csv(existing, nefin_ndvi_path)
    cat("  ✓ Updated:", nefin_ndvi_path, "\n")
  } else {
    write_csv(nefin_out, nefin_ndvi_path)
    cat("  ✓ Saved:", nefin_ndvi_path, "\n")
  }
  
  # =========================================================================
  # 4. EXTRACT S2 NDVI TO FIA
  # =========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 4: Extracting S2 NDVI to FIA Plots\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Filter to plots with coordinates
  fia_with_coords <- fia %>%
    filter(!is.na(.data[[fia_lon_col]]), !is.na(.data[[fia_lat_col]]))
  
  cat("  FIA plots with coordinates:", nrow(fia_with_coords), "\n")
  
  # Convert to sf
  fia_sf <- st_as_sf(fia_with_coords, 
                     coords = c(fia_lon_col, fia_lat_col), 
                     crs = 4326)
  fia_reproj <- st_transform(fia_sf, crs(r_s2))
  
  cat("  Extracting...\n")
  s2_vals_fia <- terra::extract(r_s2, vect(fia_reproj))
  
  fia_with_coords$ndvi_s2 <- s2_vals_fia[[2]]
  
  # Summary
  valid_s2_fia <- sum(!is.na(fia_with_coords$ndvi_s2))
  cat("  Valid extractions:", valid_s2_fia, "/", nrow(fia_with_coords),
      "(", round(100 * valid_s2_fia / nrow(fia_with_coords), 1), "%)\n")
  cat("  S2 NDVI range:", round(min(fia_with_coords$ndvi_s2, na.rm = TRUE), 3), "to",
      round(max(fia_with_coords$ndvi_s2, na.rm = TRUE), 3), "\n")
  
  # Update FIA covariates file
  cat("\n  Updating FIA covariates file...\n")
  
  # Remove old ndvi_s2 if exists
  if ("ndvi_s2" %in% names(fia)) {
    fia <- fia %>% select(-ndvi_s2)
  }
  
  # Join new S2 values
  fia_updated <- fia %>%
    left_join(fia_with_coords %>% select(CN, ndvi_s2), by = "CN")
  
  write_csv(fia_updated, fia_cov_path)
  cat("  ✓ Updated:", fia_cov_path, "\n")
  
  # =========================================================================
  # 5. SUMMARY
  # =========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  EXTRACTION COMPLETE                                                      ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  cat("Summary:\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  cat("  S2 raster used:", basename(s2_path), "\n")
  cat("  NEFIN: ", valid_s2, "/", nrow(nefin), "plots with S2 NDVI\n")
  cat("  FIA:   ", valid_s2_fia, "/", nrow(fia_with_coords), "plots with S2 NDVI\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  
  cat("\nNext steps:\n")
  cat("  1. Run: Rscript R/sensor_resolution_comparison.R\n")
  cat("  2. Run: Rscript R/spatial_model_comparison.R\n")
  
  invisible(list(
    nefin_s2 = nefin %>% select(CN, ndvi_s2),
    fia_s2 = fia_with_coords %>% select(CN, ndvi_s2)
  ))
}

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  extract_s2_to_plots()
}
