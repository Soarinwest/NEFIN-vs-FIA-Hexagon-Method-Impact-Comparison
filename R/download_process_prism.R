# R/download_process_prism.R
# Download and process PRISM climate data for the study region

suppressPackageStartupMessages({
  library(prism); library(sf); library(dplyr); library(raster); library(fs)
  library(terra); library(exactextractr); library(readr); library(yaml)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

#' Download and process PRISM climate data
#' 
#' @param variables Character vector of PRISM variables to download
#'   Options: "tmean" (mean temperature), "ppt" (precipitation), 
#'            "tmax", "tmin", "tdmean", "vpdmax", "vpdmin"
#' @param years Numeric vector of years to download
#' @param temporal_period "annual" or "monthly"
#' @param study_region_path Path to study region boundary (hex grid or states)
#' @param out_dir Output directory for processed data
#' @param prism_dir Directory for raw PRISM downloads
#' @param overwrite Overwrite existing processed files
download_process_prism <- function(variables = c("tmean", "ppt"),
                                   years = 2015:2020,
                                   temporal_period = "annual",
                                   study_region_path = "data/hex/hex_grid.geojson",
                                   out_dir = "data/processed/prism",
                                   prism_dir = "data/raw/prism",
                                   overwrite = FALSE) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  PRISM Climate Data Download & Processing                 ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Check if already processed
  meta_file <- fs::path(out_dir, "prism_metadata.yml")
  if (fs::file_exists(meta_file) && !overwrite) {
    meta <- yaml::read_yaml(meta_file)
    cat("✓ PRISM data already processed\n")
    cat("  Variables:", paste(meta$variables, collapse = ", "), "\n")
    cat("  Years:", paste(range(meta$years), collapse = "-"), "\n")
    cat("  Use overwrite=TRUE to reprocess\n")
    return(invisible(out_dir))
  }
  
  # Create directories
  fs::dir_create(prism_dir, recurse = TRUE)
  fs::dir_create(out_dir, recurse = TRUE)
  
  # Set PRISM path
  prism::prism_set_dl_dir(prism_dir)
  
  cat("→ Loading study region:", study_region_path, "\n")
  study_region <- sf::st_read(study_region_path, quiet = TRUE)
  study_bbox <- sf::st_bbox(study_region)
  
  # Transform to PRISM projection (NAD83)
  study_region_4269 <- sf::st_transform(study_region, 4269)
  
  # Download PRISM data
  cat("\n→ Downloading PRISM data...\n")
  cat("  Variables:", paste(variables, collapse = ", "), "\n")
  cat("  Years:", paste(range(years), collapse = "-"), "\n")
  cat("  Period:", temporal_period, "\n")
  
  for (var in variables) {
    cat("\n  Downloading", var, "...\n")
    
    if (temporal_period == "annual") {
      for (yr in years) {
        cat("    Year", yr, "...")
        tryCatch({
          prism::get_prism_annual(type = var, years = yr, keepZip = FALSE)
          cat(" ✓\n")
        }, error = function(e) {
          cat(" ✗ Error:", e$message, "\n")
        })
      }
    } else {
      # Monthly data
      for (yr in years) {
        cat("    Year", yr, "...")
        tryCatch({
          prism::get_prism_monthlys(type = var, year = yr, mon = 1:12, keepZip = FALSE)
          cat(" ✓\n")
        }, error = function(e) {
          cat(" ✗ Error:", e$message, "\n")
        })
      }
    }
  }
  
  # Process data to create multi-year summaries
  cat("\n→ Processing PRISM data...\n")
  
  results <- list()
  
  for (var in variables) {
    cat("  Processing", var, "...\n")
    
    # Get file paths
    if (temporal_period == "annual") {
      pattern <- paste0("PRISM_", var, "_.*_bil.bil$")
    } else {
      pattern <- paste0("PRISM_", var, "_.*_bil.bil$")
    }
    
    prism_files <- list.files(prism_dir, pattern = pattern, 
                              full.names = TRUE, recursive = TRUE)
    
    # Filter to requested years
    year_pattern <- paste0("(", paste(years, collapse = "|"), ")")
    prism_files <- prism_files[grep(year_pattern, prism_files)]
    
    if (length(prism_files) == 0) {
      cat("    ⚠ No files found for", var, "\n")
      next
    }
    
    cat("    Found", length(prism_files), "files\n")
    
    # Stack rasters using terra (faster than raster)
    rast_stack <- terra::rast(prism_files)
    
    # Calculate statistics across years
    cat("    Calculating multi-year statistics...\n")
    
    # Mean across all years
    mean_rast <- terra::mean(rast_stack, na.rm = TRUE)
    names(mean_rast) <- paste0(var, "_mean")
    
    # SD across years (temporal variability)
    sd_rast <- terra::stdev(rast_stack, na.rm = TRUE)
    names(sd_rast) <- paste0(var, "_sd")
    
    # CV (coefficient of variation)
    cv_rast <- sd_rast / mean_rast * 100
    names(cv_rast) <- paste0(var, "_cv")
    
    # Min and max
    min_rast <- terra::min(rast_stack, na.rm = TRUE)
    names(min_rast) <- paste0(var, "_min")
    
    max_rast <- terra::max(rast_stack, na.rm = TRUE)
    names(max_rast) <- paste0(var, "_max")
    
    # Combine
    var_stack <- c(mean_rast, sd_rast, cv_rast, min_rast, max_rast)
    
    # Crop to study region
    cat("    Cropping to study region...\n")
    var_crop <- terra::crop(var_stack, terra::vect(study_region_4269))
    
    # Save processed rasters
    out_file <- fs::path(out_dir, paste0(var, "_stats.tif"))
    terra::writeRaster(var_crop, out_file, overwrite = TRUE)
    cat("    ✓ Saved:", out_file, "\n")
    
    results[[var]] <- out_file
  }
  
  # Create metadata
  metadata <- list(
    created = as.character(Sys.time()),
    variables = variables,
    years = years,
    temporal_period = temporal_period,
    study_region = study_region_path,
    statistics = c("mean", "sd", "cv", "min", "max"),
    files = results,
    crs = "EPSG:4269 (NAD83)",
    units = list(
      tmean = "degrees C",
      tmin = "degrees C", 
      tmax = "degrees C",
      ppt = "mm",
      tdmean = "degrees C",
      vpdmin = "hPa",
      vpdmax = "hPa"
    )
  )
  
  yaml::write_yaml(metadata, meta_file)
  
  cat("\n✓ PRISM processing complete!\n")
  cat("  Output directory:", out_dir, "\n")
  cat("  Metadata:", meta_file, "\n")
  
  invisible(list(
    out_dir = out_dir,
    metadata = metadata
  ))
}

# CLI interface
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  get_arg <- function(flag, default = NULL) {
    hit <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
  }
  
  # Load config
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  # Parse arguments
  vars <- strsplit(get_arg("--vars", "tmean,ppt"), ",")[[1]]
  years_str <- get_arg("--years", "2015:2020")
  
  if (grepl(":", years_str)) {
    year_parts <- as.integer(strsplit(years_str, ":")[[1]])
    years <- seq(year_parts[1], year_parts[2])
  } else {
    years <- as.integer(strsplit(years_str, ",")[[1]])
  }
  
  download_process_prism(
    variables = vars,
    years = years,
    temporal_period = get_arg("--period", "annual"),
    study_region_path = get_arg("--region", cfg$hex_path %||% "data/hex/hex_grid.geojson"),
    out_dir = get_arg("--out", "data/processed/prism"),
    prism_dir = get_arg("--prism-dir", "data/raw/prism"),
    overwrite = "--overwrite" %in% args
  )
}