# check_and_create_missing_grids.R
suppressPackageStartupMessages({
  library(fs); library(yaml); library(dplyr)
})

# Load config
cfg <- yaml::read_yaml("configs/process.yml")

# Check which grids exist
cat("\n════════════════════════════════════════════════════════════\n")
cat("HEX GRID STATUS CHECK\n")
cat("════════════════════════════════════════════════════════════\n\n")

existing <- c()
missing <- c()

for (i in seq_along(cfg$hex_grids)) {
  grid <- cfg$hex_grids[[i]]
  exists <- fs::file_exists(grid$path)
  
  if (exists) {
    existing <- c(existing, grid$name)
    cat("✓ ", sprintf("%-10s", grid$name), ": ", grid$path, "\n", sep = "")
  } else {
    missing <- c(missing, grid$name)
    cat("✗ ", sprintf("%-10s", grid$name), ": ", grid$path, " [MISSING]\n", sep = "")
  }
}

cat("\n")
cat("Summary:\n")
cat("  Existing: ", length(existing), " grids\n")
cat("  Missing:  ", length(missing), " grids\n")

if (length(missing) > 0) {
  cat("\nMissing grids: ", paste(missing, collapse = ", "), "\n")
  
  # Find the hex_sizes entries for missing grids
  missing_specs <- list()
  for (spec in cfg$hex_sizes) {
    # Extract grid name from output path
    grid_name <- gsub(".*hex_grid_(.+)\\.geojson", "\\1", spec$out_path)
    grid_name <- gsub("ac$", "", grid_name)  # Remove 'ac' suffix
    
    # Check if this matches any missing grid
    for (m in missing) {
      if (grepl(grid_name, m, ignore.case = TRUE) || 
          grepl(gsub("k$", "", m), grid_name, ignore.case = TRUE)) {
        missing_specs[[length(missing_specs) + 1]] <- spec
        break
      }
    }
  }
  
  if (length(missing_specs) > 0) {
    cat("\nWould you like to create the ", length(missing_specs), " missing grids? (y/n): ")
    answer <- readline()
    
    if (tolower(answer) == "y") {
      source("R/create_hex_grid.R")
      
      cat("\nCreating missing grids...\n")
      
      for (i in seq_along(missing_specs)) {
        spec <- missing_specs[[i]]
        cat("\n[", i, "/", length(missing_specs), "] Creating ", spec$out_path, "\n", sep = "")
        
        create_hex_grid(
          area_acres = spec$area_acres,
          clip_to = spec$clip_to,
          exclude_water = spec$exclude_water,
          buffer_miles = spec$buffer_miles %||% NULL,
          out_path = spec$out_path,
          overwrite = FALSE
        )
      }
    }
  }
}

cat("\n")