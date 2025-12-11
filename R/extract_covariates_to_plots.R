#!/usr/bin/env Rscript
# ==============================================================================
# extract_covariates_to_plots.R
# ==============================================================================
# Purpose:
#   Extract PRISM climate (tmean, ppt) and NDVI (MODIS, Sentinel-2) values
#   to plot locations for FIA (jittered) and NEFIN (true coords).
#
# Fixed: Uses terra for all spatial operations to avoid CRS transformation issues
#
# Usage:
#   Rscript R/extract_covariates_to_plots.R [--overwrite] [--max-reps=N]
#
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(fs)
  library(yaml)
  library(glue)
  library(purrr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==============================================================================
# Hardcoded Paths (matching your actual structure)
# ==============================================================================

PATHS <- list(
  # PRISM climate
  prism_dir = "data/processed/prism",
  prism_tmean_pattern = "prism_tmean_ne_{start}_{end}.tif",
  prism_ppt_pattern = "prism_ppt_ne_{start}_{end}.tif",
  
  # MODIS NDVI
  modis_dir = "data/processed/ndvi/modis",
  modis_pattern = "MODIS_NDVI_5yr_blocked_{start}_{end}.tif",
  
  # Sentinel-2 NDVI
  s2_dir = "data/processed/ndvi/s2",
  s2_pattern = "ndvi_{start}_{end}_raw.tif",
  
  # Plot data
  jitter_dir = "data/processed/mc_jitter_library/replicates",
  nefin_file = "data/processed/nefin_processed.csv",
  fia_plot_file = "data/interim/fia_region/plot.csv",
  
  # Outputs
  ndvi_out_dir = "data/processed/ndvi_at_plots",
  climate_out_dir = "data/processed/climate_at_plots"
)

# Temporal windows (5-year blocks)
WINDOWS <- list(
  list(name = "2000_2004", start = 2000, end = 2004, 
       has_modis = TRUE, has_s2 = FALSE, has_prism = TRUE),
  list(name = "2005_2009", start = 2005, end = 2009, 
       has_modis = TRUE, has_s2 = FALSE, has_prism = TRUE),
  list(name = "2010_2014", start = 2010, end = 2014, 
       has_modis = TRUE, has_s2 = FALSE, has_prism = TRUE),
  list(name = "2015_2019", start = 2015, end = 2019, 
       has_modis = TRUE, has_s2 = FALSE, has_prism = TRUE),
  list(name = "2016_2019", start = 2016, end = 2019, 
       has_modis = FALSE, has_s2 = TRUE, has_prism = FALSE),
  list(name = "2020_2024", start = 2020, end = 2024, 
       has_modis = TRUE, has_s2 = TRUE, has_prism = TRUE)
)

# ==============================================================================
# Raster Discovery
# ==============================================================================

discover_rasters <- function() {
  
  rasters <- list()
  
  for (w in WINDOWS) {
    # PRISM tmean
    if (w$has_prism) {
      path <- fs::path(PATHS$prism_dir, 
                       glue(PATHS$prism_tmean_pattern, start = w$start, end = w$end))
      rasters[[length(rasters) + 1]] <- list(
        product = "prism", variable = "tmean",
        window = w$name, start = w$start, end = w$end,
        path = as.character(path), exists = fs::file_exists(path)
      )
      
      # PRISM ppt
      path <- fs::path(PATHS$prism_dir,
                       glue(PATHS$prism_ppt_pattern, start = w$start, end = w$end))
      rasters[[length(rasters) + 1]] <- list(
        product = "prism", variable = "ppt",
        window = w$name, start = w$start, end = w$end,
        path = as.character(path), exists = fs::file_exists(path)
      )
    }
    
    # MODIS NDVI
    if (w$has_modis) {
      path <- fs::path(PATHS$modis_dir,
                       glue(PATHS$modis_pattern, start = w$start, end = w$end))
      rasters[[length(rasters) + 1]] <- list(
        product = "modis", variable = "ndvi",
        window = w$name, start = w$start, end = w$end,
        path = as.character(path), exists = fs::file_exists(path)
      )
    }
    
    # Sentinel-2 NDVI
    if (w$has_s2) {
      path <- fs::path(PATHS$s2_dir,
                       glue(PATHS$s2_pattern, start = w$start, end = w$end))
      rasters[[length(rasters) + 1]] <- list(
        product = "s2", variable = "ndvi",
        window = w$name, start = w$start, end = w$end,
        path = as.character(path), exists = fs::file_exists(path)
      )
    }
  }
  
  # Convert to data frame
  bind_rows(lapply(rasters, as.data.frame))
}

# ==============================================================================
# Year-to-Window Mapping
# ==============================================================================

map_year_to_window <- function(year) {
  if (is.na(year)) return(NA_character_)
  
  # Standard 5-year blocks for MODIS/PRISM
  if (year >= 2000 && year <= 2004) return("2000_2004")
  if (year >= 2005 && year <= 2009) return("2005_2009")
  if (year >= 2010 && year <= 2014) return("2010_2014")
  if (year >= 2015 && year <= 2019) return("2015_2019")
  if (year >= 2020 && year <= 2024) return("2020_2024")
  
 # Outside known windows - find closest
  if (year < 2000) return("2000_2004")
  if (year > 2024) return("2020_2024")
  
  NA_character_
}

# For S2, need special handling (only 2016-2019 and 2020-2024)
map_year_to_s2_window <- function(year) {
  if (is.na(year)) return(NA_character_)
  if (year >= 2016 && year <= 2019) return("2016_2019")
  if (year >= 2020 && year <= 2024) return("2020_2024")
  NA_character_
}

# ==============================================================================
# Extraction Function (using terra only - avoids sf CRS issues)
# ==============================================================================

extract_raster_values_terra <- function(lon, lat, raster_path) {
  # Extract values using terra only (no sf transformation)
  
  if (!fs::file_exists(raster_path)) {
    return(rep(NA_real_, length(lon)))
  }
  
  r <- tryCatch({
    terra::rast(raster_path)
  }, error = function(e) {
    warning("Failed to load: ", raster_path, " - ", e$message)
    return(NULL)
  })
  
  if (is.null(r)) return(rep(NA_real_, length(lon)))
  
  # Create points in WGS84
  pts <- terra::vect(
    data.frame(x = lon, y = lat),
    geom = c("x", "y"),
    crs = "EPSG:4326"
  )
  
  # Project points to raster CRS using terra (more robust)
  raster_crs <- terra::crs(r)
  
  pts_proj <- tryCatch({
    terra::project(pts, raster_crs)
  }, error = function(e) {
    # If projection fails, try using EPSG:5070 (Albers) as fallback
    message("    CRS projection issue, trying EPSG:5070 fallback...")
    tryCatch({
      terra::project(pts, "EPSG:5070")
    }, error = function(e2) {
      warning("Could not project points: ", e2$message)
      return(NULL)
    })
  })
  
  if (is.null(pts_proj)) return(rep(NA_real_, length(lon)))
  
  # Extract values
  vals <- tryCatch({
    terra::extract(r, pts_proj)
  }, error = function(e) {
    warning("Extraction failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(vals) || ncol(vals) < 2) {
    return(rep(NA_real_, length(lon)))
  }
  
  as.numeric(vals[[2]])
}

# ==============================================================================
# Batch extraction for a single raster (more efficient)
# ==============================================================================

extract_single_raster <- function(lon, lat, raster_path, chunk_size = 10000) {
  # For large datasets, extract in chunks to avoid memory issues
  
  n <- length(lon)
  if (n == 0) return(numeric(0))
  
  if (n <= chunk_size) {
    return(extract_raster_values_terra(lon, lat, raster_path))
  }
  
  # Process in chunks
  result <- numeric(n)
  n_chunks <- ceiling(n / chunk_size)
  
  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n)
    
    result[start_idx:end_idx] <- extract_raster_values_terra(
      lon[start_idx:end_idx],
      lat[start_idx:end_idx],
      raster_path
    )
  }
  
  result
}

# ==============================================================================
# FIA Extraction
# ==============================================================================

extract_fia <- function(raster_inv, max_reps = NULL) {
  
  message("\n", strrep("=", 60))
  message("EXTRACTING FOR FIA PLOTS (JITTERED)")
  message(strrep("=", 60), "\n")
  
  # Get replicate files
  rep_files <- list.files(PATHS$jitter_dir, pattern = "^rep_\\d+\\.csv$", 
                          full.names = TRUE)
  rep_files <- sort(rep_files)
  message("Found ", length(rep_files), " jitter replicates")
  
  if (!is.null(max_reps) && max_reps < length(rep_files)) {
    rep_files <- rep_files[1:max_reps]
    message("Processing first ", max_reps, " only")
  }
  
  # Load FIA plot metadata for MEASYEAR
  if (fs::file_exists(PATHS$fia_plot_file)) {
    plot_meta <- read_csv(PATHS$fia_plot_file, show_col_types = FALSE) %>%
      select(CN, STATECD, MEASYEAR) %>%
      distinct()
    message("Plot metadata: ", nrow(plot_meta), " plots")
  } else {
    stop("FIA plot.csv not found: ", PATHS$fia_plot_file)
  }
  
  # Filter available rasters
  available <- raster_inv %>% filter(exists)
  message("\nAvailable rasters: ", nrow(available))
  
  # Pre-load raster paths by window for efficiency
  raster_lookup <- list()
  for (i in 1:nrow(available)) {
    key <- paste(available$product[i], available$variable[i], available$window[i], sep = "_")
    raster_lookup[[key]] <- available$path[i]
  }
  
  # Storage
  all_results <- list()
  
  start_time <- Sys.time()
  
  for (rep_idx in seq_along(rep_files)) {
    rep_file <- rep_files[rep_idx]
    rep_id <- as.integer(gsub(".*rep_(\\d+)\\.csv", "\\1", basename(rep_file)))
    
    # Progress
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    if (rep_idx > 1) {
      rate <- (rep_idx - 1) / elapsed
      eta <- (length(rep_files) - rep_idx) / rate / 60
      message(sprintf("  Rep %d/%d (%.0f%%) - ETA: %.1f min",
                      rep_idx, length(rep_files),
                      100 * rep_idx / length(rep_files), eta))
    } else {
      message(sprintf("  Rep %d/%d", rep_idx, length(rep_files)))
    }
    
    # Load replicate
    rep_data <- read_csv(rep_file, show_col_types = FALSE)
    
    # Add replicate_id if not present
    if (!"replicate_id" %in% names(rep_data)) {
      rep_data$replicate_id <- rep_id
    }
    
    # Join MEASYEAR if not present
    if (!"MEASYEAR" %in% names(rep_data)) {
      rep_data <- rep_data %>%
        left_join(plot_meta %>% select(CN, MEASYEAR), by = "CN")
    }
    
    # Map to windows
    rep_data <- rep_data %>%
      mutate(
        window_modis = sapply(MEASYEAR, map_year_to_window),
        window_s2 = sapply(MEASYEAR, map_year_to_s2_window)
      )
    
    # Filter valid coords
    valid_idx <- !is.na(rep_data$lon_jittered) & !is.na(rep_data$lat_jittered)
    if (sum(valid_idx) == 0) {
      message("    No valid coordinates, skipping")
      next
    }
    
    rep_valid <- rep_data[valid_idx, ]
    
    # Initialize result columns
    rep_valid$tmean <- NA_real_
    rep_valid$ppt <- NA_real_
    rep_valid$ndvi_modis <- NA_real_
    rep_valid$ndvi_s2 <- NA_real_
    
    # Extract by window (MODIS/PRISM use same windows)
    for (window in unique(na.omit(rep_valid$window_modis))) {
      idx <- which(rep_valid$window_modis == window & !is.na(rep_valid$window_modis))
      if (length(idx) == 0) next
      
      lon <- rep_valid$lon_jittered[idx]
      lat <- rep_valid$lat_jittered[idx]
      
      # PRISM tmean
      key <- paste("prism", "tmean", window, sep = "_")
      if (!is.null(raster_lookup[[key]])) {
        rep_valid$tmean[idx] <- extract_single_raster(lon, lat, raster_lookup[[key]])
      }
      
      # PRISM ppt
      key <- paste("prism", "ppt", window, sep = "_")
      if (!is.null(raster_lookup[[key]])) {
        rep_valid$ppt[idx] <- extract_single_raster(lon, lat, raster_lookup[[key]])
      }
      
      # MODIS NDVI
      key <- paste("modis", "ndvi", window, sep = "_")
      if (!is.null(raster_lookup[[key]])) {
        rep_valid$ndvi_modis[idx] <- extract_single_raster(lon, lat, raster_lookup[[key]])
      }
    }
    
    # S2 extraction (different windows)
    for (window in unique(na.omit(rep_valid$window_s2))) {
      idx <- which(rep_valid$window_s2 == window & !is.na(rep_valid$window_s2))
      if (length(idx) == 0) next
      
      lon <- rep_valid$lon_jittered[idx]
      lat <- rep_valid$lat_jittered[idx]
      
      key <- paste("s2", "ndvi", window, sep = "_")
      if (!is.null(raster_lookup[[key]])) {
        rep_valid$ndvi_s2[idx] <- extract_single_raster(lon, lat, raster_lookup[[key]])
      }
    }
    
    # Store results
    all_results[[rep_idx]] <- rep_valid
  }
  
  # Combine all results
  combined <- bind_rows(all_results)
  
  message("\nFIA extraction complete:")
  message("  Total records: ", format(nrow(combined), big.mark = ","))
  message("  With NDVI (MODIS): ", format(sum(!is.na(combined$ndvi_modis)), big.mark = ","))
  message("  With climate: ", format(sum(!is.na(combined$tmean)), big.mark = ","))
  
  # Split into climate and NDVI outputs
  base_cols <- c("CN", "STATECD", "MEASYEAR", "replicate_id", 
                 "lon_jittered", "lat_jittered", "window_modis")
  base_cols <- intersect(base_cols, names(combined))
  
  # Add hex columns if present
  hex_cols <- names(combined)[grepl("^hex_id_", names(combined))]
  base_cols <- c(base_cols, hex_cols)
  
  climate_df <- combined %>% select(all_of(base_cols), tmean, ppt)
  ndvi_df <- combined %>% select(all_of(base_cols), ndvi_modis, ndvi_s2)
  
  list(climate = climate_df, ndvi = ndvi_df)
}

# ==============================================================================
# NEFIN Extraction
# ==============================================================================

extract_nefin <- function(raster_inv) {
  
  message("\n", strrep("=", 60))
  message("EXTRACTING FOR NEFIN PLOTS (TRUE COORDS)")
  message(strrep("=", 60), "\n")
  
  if (!fs::file_exists(PATHS$nefin_file)) {
    message("NEFIN file not found: ", PATHS$nefin_file)
    return(list(climate = NULL, ndvi = NULL))
  }
  
  nefin <- read_csv(PATHS$nefin_file, show_col_types = FALSE)
  message("NEFIN plots: ", nrow(nefin))
  
  # Identify coordinate columns
  if ("lon_public" %in% names(nefin)) {
    lon_col <- "lon_public"
    lat_col <- "lat_public"
  } else if ("lon" %in% names(nefin)) {
    lon_col <- "lon"
    lat_col <- "lat"
  } else {
    stop("Cannot find coordinate columns in NEFIN data")
  }
  message("Using coordinates: ", lon_col, ", ", lat_col)
  
  # Map years to windows
  nefin <- nefin %>%
    mutate(
      window_modis = sapply(MEASYEAR, map_year_to_window),
      window_s2 = sapply(MEASYEAR, map_year_to_s2_window)
    )
  
  # Filter valid coords
  valid_idx <- !is.na(nefin[[lon_col]]) & !is.na(nefin[[lat_col]])
  nefin <- nefin[valid_idx, ]
  message("Valid coordinates: ", nrow(nefin))
  
  # Initialize result columns
  nefin$tmean <- NA_real_
  nefin$ppt <- NA_real_
  nefin$ndvi_modis <- NA_real_
  nefin$ndvi_s2 <- NA_real_
  
  # Pre-load raster paths
  available <- raster_inv %>% filter(exists)
  raster_lookup <- list()
  for (i in 1:nrow(available)) {
    key <- paste(available$product[i], available$variable[i], available$window[i], sep = "_")
    raster_lookup[[key]] <- available$path[i]
  }
  
  # Extract by window (MODIS/PRISM)
  for (window in unique(na.omit(nefin$window_modis))) {
    message("  Processing window: ", window)
    idx <- which(nefin$window_modis == window & !is.na(nefin$window_modis))
    if (length(idx) == 0) next
    
    lon <- nefin[[lon_col]][idx]
    lat <- nefin[[lat_col]][idx]
    
    # PRISM tmean
    key <- paste("prism", "tmean", window, sep = "_")
    if (!is.null(raster_lookup[[key]])) {
      nefin$tmean[idx] <- extract_single_raster(lon, lat, raster_lookup[[key]])
    }
    
    # PRISM ppt
    key <- paste("prism", "ppt", window, sep = "_")
    if (!is.null(raster_lookup[[key]])) {
      nefin$ppt[idx] <- extract_single_raster(lon, lat, raster_lookup[[key]])
    }
    
    # MODIS NDVI
    key <- paste("modis", "ndvi", window, sep = "_")
    if (!is.null(raster_lookup[[key]])) {
      nefin$ndvi_modis[idx] <- extract_single_raster(lon, lat, raster_lookup[[key]])
    }
  }
  
  # S2 extraction
  for (window in unique(na.omit(nefin$window_s2))) {
    message("  Processing S2 window: ", window)
    idx <- which(nefin$window_s2 == window & !is.na(nefin$window_s2))
    if (length(idx) == 0) next
    
    lon <- nefin[[lon_col]][idx]
    lat <- nefin[[lat_col]][idx]
    
    key <- paste("s2", "ndvi", window, sep = "_")
    if (!is.null(raster_lookup[[key]])) {
      nefin$ndvi_s2[idx] <- extract_single_raster(lon, lat, raster_lookup[[key]])
    }
  }
  
  message("\nNEFIN extraction complete:")
  message("  With NDVI (MODIS): ", sum(!is.na(nefin$ndvi_modis)))
  message("  With climate: ", sum(!is.na(nefin$tmean)))
  
  # Build output - keep relevant columns
  base_cols <- c("CN", "STATECD", "MEASYEAR", lon_col, lat_col, 
                 "aglb_Mg_per_ha", "source", "window_modis")
  base_cols <- intersect(base_cols, names(nefin))
  
  hex_cols <- names(nefin)[grepl("^hex_id_", names(nefin))]
  base_cols <- c(base_cols, hex_cols)
  
  climate_df <- nefin %>% select(all_of(base_cols), tmean, ppt)
  ndvi_df <- nefin %>% select(all_of(base_cols), ndvi_modis, ndvi_s2)
  
  list(climate = climate_df, ndvi = ndvi_df)
}

# ==============================================================================
# Main
# ==============================================================================

main <- function() {
  
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("  COVARIATE EXTRACTION: PRISM + NDVI to Plots\n")
  cat("  (Using terra for robust CRS handling)\n")
  cat(strrep("=", 60), "\n\n")
  
  # Parse args
  args <- commandArgs(trailingOnly = TRUE)
  overwrite <- "--overwrite" %in% args
  
  max_reps <- NULL
  max_arg <- args[grepl("^--max-reps=", args)]
  if (length(max_arg) > 0) {
    max_reps <- as.integer(gsub("^--max-reps=", "", max_arg[1]))
  }
  
  # Create output dirs
  fs::dir_create(PATHS$ndvi_out_dir, recurse = TRUE)
  fs::dir_create(PATHS$climate_out_dir, recurse = TRUE)
  
  # Discover rasters
  message("Discovering rasters...")
  raster_inv <- discover_rasters()
  
  n_exist <- sum(raster_inv$exists)
  message("Found ", n_exist, "/", nrow(raster_inv), " rasters\n")
  
  # Show what's available
  message("Available rasters:")
  available <- raster_inv %>% filter(exists)
  for (i in 1:nrow(available)) {
    message("  [", available$window[i], "] ", 
            available$product[i], "/", available$variable[i], 
            ": ", basename(available$path[i]))
  }
  
  if (n_exist == 0) {
    stop("No rasters found! Check paths.")
  }
  
  # Check a raster's CRS
  test_raster <- available$path[1]
  r <- terra::rast(test_raster)
  message("\nRaster CRS: ", terra::crs(r, describe = TRUE)$name)
  
  # Extract FIA
  fia <- extract_fia(raster_inv, max_reps)
  
  if (!is.null(fia$climate) && nrow(fia$climate) > 0) {
    path <- fs::path(PATHS$climate_out_dir, "fia_climate.csv")
    write_csv(fia$climate, path)
    message("\nSaved: ", path, " (", nrow(fia$climate), " rows)")
  }
  
  if (!is.null(fia$ndvi) && nrow(fia$ndvi) > 0) {
    path <- fs::path(PATHS$ndvi_out_dir, "fia_ndvi.csv")
    write_csv(fia$ndvi, path)
    message("Saved: ", path, " (", nrow(fia$ndvi), " rows)")
  }
  
  # Extract NEFIN
  nefin <- extract_nefin(raster_inv)
  
  if (!is.null(nefin$climate) && nrow(nefin$climate) > 0) {
    path <- fs::path(PATHS$climate_out_dir, "nefin_climate.csv")
    write_csv(nefin$climate, path)
    message("\nSaved: ", path, " (", nrow(nefin$climate), " rows)")
  }
  
  if (!is.null(nefin$ndvi) && nrow(nefin$ndvi) > 0) {
    path <- fs::path(PATHS$ndvi_out_dir, "nefin_ndvi.csv")
    write_csv(nefin$ndvi, path)
    message("Saved: ", path, " (", nrow(nefin$ndvi), " rows)")
  }
  
  # Manifest
  manifest <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    rasters_found = n_exist,
    rasters_expected = nrow(raster_inv),
    fia_replicates_processed = if (!is.null(max_reps)) max_reps else length(list.files(PATHS$jitter_dir, pattern = "\\.csv$")),
    fia_records = if (!is.null(fia$ndvi)) nrow(fia$ndvi) else 0,
    nefin_records = if (!is.null(nefin$ndvi)) nrow(nefin$ndvi) else 0,
    outputs = list(
      fia_climate = as.character(fs::path(PATHS$climate_out_dir, "fia_climate.csv")),
      fia_ndvi = as.character(fs::path(PATHS$ndvi_out_dir, "fia_ndvi.csv")),
      nefin_climate = as.character(fs::path(PATHS$climate_out_dir, "nefin_climate.csv")),
      nefin_ndvi = as.character(fs::path(PATHS$ndvi_out_dir, "nefin_ndvi.csv"))
    )
  )
  
  write_yaml(manifest, "data/processed/extraction_manifest.yml")
  message("\nManifest saved: data/processed/extraction_manifest.yml")
  
  cat("\n", strrep("=", 60), "\n")
  cat("  EXTRACTION COMPLETE\n")
  cat(strrep("=", 60), "\n\n")
}

if (!interactive()) {
  main()
}
