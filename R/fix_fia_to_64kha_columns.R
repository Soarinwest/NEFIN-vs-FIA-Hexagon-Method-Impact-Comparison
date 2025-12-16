#!/usr/bin/env Rscript
# R/fix_fia_to_64kha_columns.R
# Rename hex_id_fia column to hex_id_64kha in all relevant data files

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(fs)
})

fix_fia_column_names <- function() {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Fix Column Names: hex_id_fia → hex_id_64kha             ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  files_to_check <- c(
    "data/processed/plot_hex_assignments.csv",
    "data/processed/nefin_hex_assignments.csv",
    "data/processed/fia_covariates.csv",
    "data/processed/climate_at_plots/fia_climate.csv",
    "data/processed/climate_at_plots/nefin_climate.csv",
    "data/processed/ndvi_at_plots/fia_ndvi.csv",
    "data/processed/ndvi_at_plots/nefin_ndvi.csv"
  )
  
 # Also check jitter library replicates
  jitter_dir <- "data/processed/mc_jitter_library/replicates"
  if (fs::dir_exists(jitter_dir)) {
    jitter_files <- list.files(jitter_dir, pattern = "rep_.*\\.csv$", full.names = TRUE)
    files_to_check <- c(files_to_check, jitter_files)
  }
  
  files_updated <- 0
  
  for (f in files_to_check) {
    if (!fs::file_exists(f)) {
      next
    }
    
    # Read header only first to check columns
    header <- names(readr::read_csv(f, n_max = 0, show_col_types = FALSE))
    
    if ("hex_id_fia" %in% header) {
      cat("  Updating:", f, "\n")
      
      # Read full file
      df <- readr::read_csv(f, show_col_types = FALSE)
      
      # Rename column
      df <- df %>% rename(hex_id_64kha = hex_id_fia)
      
      # Write back
      readr::write_csv(df, f)
      
      files_updated <- files_updated + 1
    }
  }
  
  cat("\n")
  cat("══════════════════════════════════════════════════════════\n")
  cat("Summary:\n")
  cat("  Files checked:", length(files_to_check), "\n")
  cat("  Files updated:", files_updated, "\n")
  cat("══════════════════════════════════════════════════════════\n")
  
  if (files_updated > 0) {
    cat("\n✓ Column renaming complete!\n")
    cat("  You can now run: Rscript run_pipeline.R --all-scales\n\n")
  } else {
    cat("\n✓ No files needed updating (already using hex_id_64kha)\n\n")
  }
}

# Run
fix_fia_column_names()
