# R/create_all_hex_scales.R
# Generate full range of hex grids for scale threshold analysis

suppressPackageStartupMessages({
  library(fs); library(yaml)
})

source("R/create_hex_grid.R")

create_all_hex_scales <- function(overwrite = FALSE) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Creating Hex Grids for Scale Threshold Analysis         ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Define all hex sizes
  hex_specs <- list(
    #list(area_acres = 532, label = "532ac"),      # 1.5 km
    list(area_acres = 946, label = "946ac"),      # 2 km
    list(area_acres = 2128, label = "2.1kac"),    # 3 km
    list(area_acres = 3785, label = "3.8kac"),    # 4 km
    list(area_acres = 5914, label = "5.9kac"),    # 5 km
    list(area_acres = 8514, label = "8.5kac"),    # 6 km
    list(area_acres = 15138, label = "15.1kac"),  # 8 km
    list(area_acres = 23668, label = "23.7kac"),  # 10 km
    list(area_acres = 53254, label = "53.3kac"),  # 15 km
    list(area_acres = 94671, label = "94.7kac"),  # 20 km
    list(area_acres = 213008, label = "213kac"),  # 30 km
    list(area_acres = 378684, label = "379kac"),  # 40 km
    list(area_acres = 591537, label = "592kac"),  # 50 km
    list(area_acres = 1514736, label = "1515kac"), # 80 km
    list(area_acres = 2366855, label = "2367kac"),  # 100 km
    list(area_acres = 6059008,  label = "6059kac"),   # 160 km
    list(area_acres = 24236032, label = "24236kac"),  # 320 km
    list(area_acres = 54531072, label = "54531kac")  # 480 km (optional)
  )
  
  # Get clipping region (use existing FIA grid as boundary)
  clip_to <- "data/hex/hex_grid.geojson"
  if (!fs::file_exists(clip_to)) {
    stop("Need existing hex grid as clipping boundary: ", clip_to)
  }
  
  # Water bodies for exclusion
  water_path <- "data/hex/Waterbodies_FeaturesToJSON.geojson"
  exclude_water <- if (fs::file_exists(water_path)) water_path else NULL
  
  created <- 0
  skipped <- 0
  
  for (spec in hex_specs) {
    out_path <- fs::path("data/hex", paste0("hex_grid_", spec$label, ".geojson"))
    
    if (fs::file_exists(out_path) && !overwrite) {
      cat("  ✓ Exists: ", spec$label, " (", spec$area_acres, " acres)\n", sep = "")
      skipped <- skipped + 1
    } else {
      cat("  → Creating: ", spec$label, " (", spec$area_acres, " acres)...\n", sep = "")
      
      tryCatch({
        create_hex_grid(
          area_acres = spec$area_acres,
          clip_to = clip_to,
          exclude_water = exclude_water,
          buffer_miles = 5,
          out_path = out_path,
          overwrite = TRUE
        )
        created <- created + 1
      }, error = function(e) {
        cat("    ✗ Error: ", e$message, "\n")
      })
    }
  }
  
  cat("\n")
  cat("Summary:\n")
  cat("  Created: ", created, " grids\n")
  cat("  Skipped: ", skipped, " grids (already exist)\n")
  cat("  Total:   ", length(hex_specs), " grids\n")
  
  invisible(TRUE)
}

# Run if called directly
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  overwrite <- "--overwrite" %in% args
  
  create_all_hex_scales(overwrite = overwrite)
}
