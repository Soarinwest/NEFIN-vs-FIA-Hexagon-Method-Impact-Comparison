# create_all_grids_sequential.R
suppressPackageStartupMessages({
  library(fs); library(yaml)
})

source("R/create_hex_grid.R")

# Load config
cfg <- yaml::read_yaml("configs/process.yml")

cat("\n════════════════════════════════════════════════════════════\n")
cat("CREATING ALL HEX GRIDS SEQUENTIALLY\n")
cat("════════════════════════════════════════════════════════════\n\n")

# Get command line args
args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
skip_existing <- "--skip-existing" %in% args

# Count existing
n_existing <- 0
for (spec in cfg$hex_sizes) {
  if (fs::file_exists(spec$out_path)) {
    n_existing <- n_existing + 1
  }
}

cat("Total grids to create: ", length(cfg$hex_sizes), "\n")
cat("Already exist: ", n_existing, "\n")
cat("Overwrite mode: ", if(overwrite) "YES" else "NO", "\n")
cat("Skip existing: ", if(skip_existing) "YES" else "NO", "\n\n")

if (n_existing == length(cfg$hex_sizes) && !overwrite) {
  cat("All grids already exist! Use --overwrite to recreate.\n")
  quit(status = 0)
}

start_time <- Sys.time()
created <- 0
skipped <- 0
errors <- 0

for (i in seq_along(cfg$hex_sizes)) {
  spec <- cfg$hex_sizes[[i]]
  
  cat("\n════════════════════════════════════════════════════════════\n")
  cat("[", i, "/", length(cfg$hex_sizes), "] Grid: ", spec$area_acres, " acres\n", sep = "")
  cat("════════════════════════════════════════════════════════════\n")
  
  if (fs::file_exists(spec$out_path) && skip_existing) {
    cat("✓ Already exists: ", spec$out_path, "\n")
    cat("  Skipping...\n")
    skipped <- skipped + 1
    next
  }
  
  result <- tryCatch({
    create_hex_grid(
      area_acres = spec$area_acres,
      clip_to = spec$clip_to,
      exclude_water = spec$exclude_water,
      buffer_miles = spec$buffer_miles,  # This will be NULL from your config
      out_path = spec$out_path,
      overwrite = overwrite
    )
    created <- created + 1
    TRUE
  }, error = function(e) {
    cat("\n✗ ERROR: ", e$message, "\n")
    errors <- errors + 1
    FALSE
  })
}

elapsed <- difftime(Sys.time(), start_time, units = "mins")

cat("\n\n════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("════════════════════════════════════════════════════════════\n")
cat("  Total time: ", round(elapsed, 1), " minutes\n")
cat("  Created: ", created, " grids\n")
cat("  Skipped: ", skipped, " grids\n")
cat("  Errors: ", errors, " grids\n")
cat("  Success rate: ", round(100 * created / (created + errors), 1), "%\n")
cat("════════════════════════════════════════════════════════════\n\n")