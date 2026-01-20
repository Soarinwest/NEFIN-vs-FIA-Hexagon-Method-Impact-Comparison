# R/extract_prism_to_hex.R
# Extract PRISM climate statistics to all hex grid scales

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(readr)
  library(exactextractr); library(fs); library(yaml)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

extract_prism_to_hex <- function(prism_dir = "data/processed/prism",
                                 hex_grids = NULL,
                                 out_file = "data/processed/hex_prism_values.csv",
                                 overwrite = FALSE) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Extract PRISM Climate Data to Hex Grids                 ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  if (fs::file_exists(out_file) && !overwrite) {
    cat("✓ Hex PRISM values already exist:", out_file, "\n")
    cat("  Use overwrite=TRUE to regenerate\n")
    return(invisible(out_file))
  }
  
  # Load config to get hex grids
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  if (is.null(hex_grids)) {
    hex_grids <- cfg$hex_grids
    if (is.null(hex_grids)) {
      stop("No hex_grids specified. Set in function call or configs/process.yml")
    }
  }
  
  # Load PRISM metadata
  prism_meta_file <- fs::path(prism_dir, "prism_metadata.yml")
  if (!fs::file_exists(prism_meta_file)) {
    stop("PRISM metadata not found. Run download_process_prism.R first")
  }
  
  prism_meta <- yaml::read_yaml(prism_meta_file)
  
  cat("→ PRISM variables:", paste(prism_meta$variables, collapse = ", "), "\n")
  cat("→ Processing", length(hex_grids), "hex grid scales\n\n")
  
  all_results <- list()
  
  # Process each hex grid
  for (i in seq_along(hex_grids)) {
    grid_info <- hex_grids[[i]]
    grid_name <- grid_info$name
    
    cat("Processing", grid_name, "grid...\n")
    
    # Load hex grid
    hex_sf <- sf::st_read(grid_info$path, quiet = TRUE)
    
    # Ensure hex_id column
    if (!("hex_id" %in% names(hex_sf))) {
      if ("ID" %in% names(hex_sf)) hex_sf <- dplyr::rename(hex_sf, hex_id = ID)
      else if ("OBJECTID" %in% names(hex_sf)) hex_sf <- dplyr::rename(hex_sf, hex_id = OBJECTID)
      else hex_sf$hex_id <- sprintf("H%06d", seq_len(nrow(hex_sf)))
    }
    
    hex_sf$hex_id <- as.character(hex_sf$hex_id)
    
    # Transform to PRISM CRS (NAD83)
    hex_4269 <- sf::st_transform(hex_sf, 4269)
    
    # Process each PRISM variable
    hex_data <- hex_4269 |> sf::st_drop_geometry() |> dplyr::select(hex_id)
    
    for (var in prism_meta$variables) {
      raster_file <- prism_meta$files[[var]]
      
      if (!fs::file_exists(raster_file)) {
        cat("  ⚠ Missing:", raster_file, "\n")
        next
      }
      
      cat("  Extracting", var, "...")
      
      # Load raster
      r <- terra::rast(raster_file)
      
      # Extract values for each statistic
      stats <- c("mean", "sd", "cv", "min", "max")
      
      for (stat in stats) {
        layer_name <- paste0(var, "_", stat)
        
        if (layer_name %in% names(r)) {
          # Extract mean value per hex
          extracted <- exactextractr::exact_extract(
            r[[layer_name]], 
            hex_4269, 
            'mean',
            progress = FALSE
          )
          
          col_name <- paste0(var, "_", stat, "_", grid_name)
          hex_data[[col_name]] <- extracted
          
          # Also extract spatial SD within hex (heterogeneity)
          extracted_sd <- exactextractr::exact_extract(
            r[[layer_name]], 
            hex_4269, 
            'stdev',
            progress = FALSE
          )
          
          col_name_sd <- paste0(var, "_", stat, "_spatial_sd_", grid_name)
          hex_data[[col_name_sd]] <- extracted_sd
        }
      }
      
      cat(" ✓\n")
    }
    
    # Add grid scale identifier
    hex_data$grid_scale <- grid_name
    
    # Pivot to long format for easier joining
    hex_long <- hex_data |>
      tidyr::pivot_longer(
        cols = -c(hex_id, grid_scale),
        names_to = "variable",
        values_to = "value"
      ) |>
      dplyr::filter(!is.na(value))
    
    all_results[[grid_name]] <- hex_long
  }
  
  # Combine all scales
  cat("\n→ Combining results...\n")
  combined <- dplyr::bind_rows(all_results)
  
  # Pivot back to wide format with all scales
  final_wide <- combined |>
    tidyr::pivot_wider(
      names_from = c(variable, grid_scale),
      values_from = value,
      names_sep = "_"
    )
  
  # Write output
  fs::dir_create(fs::path_dir(out_file), recurse = TRUE)
  readr::write_csv(final_wide, out_file)
  
  cat("\n✓ Wrote:", out_file, "\n")
  
  # Summary statistics
  cat("\nSummary by grid scale:\n")
  summary_stats <- combined |>
    dplyr::group_by(grid_scale, variable) |>
    dplyr::summarise(
      n_hexes = dplyr::n_distinct(hex_id),
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(as.data.frame(summary_stats), row.names = FALSE)
  
  # Create summary file
  summary_file <- fs::path(fs::path_dir(out_file), "hex_prism_summary.csv")
  readr::write_csv(summary_stats, summary_file)
  
  invisible(list(
    output = out_file,
    summary = summary_file,
    n_hexes = dplyr::n_distinct(combined$hex_id),
    variables = unique(combined$variable)
  ))
}

# CLI interface
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  extract_prism_to_hex(
    prism_dir = get_arg("--prism", "data/processed/prism"),
    hex_grids = cfg$hex_grids,
    out_file = get_arg("--out", "data/processed/hex_prism_values.csv"),
    overwrite = "--overwrite" %in% args
  )
}