#!/usr/bin/env Rscript
# R/batch_rename_fia_to_64kha.R
# One-time script to convert all existing "fia" scale names to "64kha"
# Run this once after the pipeline to update all CSVs

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(fs)
})

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  Batch Rename: fia → 64kha in all CSVs                   ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n")
cat("\n")

# Directories to process
dirs_to_process <- c("runs", "data/processed")

# Columns that might contain scale names
scale_columns <- c("grid_scale", "scale_name", "scale", "grid")

# Track changes
files_updated <- character(0)
files_checked <- 0

for (dir in dirs_to_process) {
  if (!fs::dir_exists(dir)) {
    cat("⚠ Directory not found:", dir, "\n")
    next
  }
  
  csv_files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  cat("Checking", length(csv_files), "CSV files in", dir, "...\n")
  
  for (f in csv_files) {
    files_checked <- files_checked + 1
    
    # Read file
    df <- tryCatch({
      readr::read_csv(f, show_col_types = FALSE)
    }, error = function(e) {
      cat("  ⚠ Could not read:", basename(f), "\n")
      return(NULL)
    })
    
    if (is.null(df)) next
    
    # Check each potential scale column
    needs_update <- FALSE
    for (col in scale_columns) {
      if (col %in% names(df)) {
        # Check if "fia" exists (case insensitive)
        if (any(tolower(df[[col]]) == "fia", na.rm = TRUE)) {
          # Replace fia with 64kha
          df[[col]][tolower(df[[col]]) == "fia"] <- "64kha"
          needs_update <- TRUE
        }
      }
    }
    
    if (needs_update) {
      # Write back
      readr::write_csv(df, f)
      cat("  ✓ Updated:", f, "\n")
      files_updated <- c(files_updated, f)
    }
  }
}

cat("\n")
cat("══════════════════════════════════════════════════════════\n")
cat("Summary:\n")
cat("  Files checked:", files_checked, "\n")
cat("  Files updated:", length(files_updated), "\n")
cat("══════════════════════════════════════════════════════════\n")

if (length(files_updated) > 0) {
  cat("\nUpdated files:\n")
  for (f in files_updated) {
    cat("  -", f, "\n")
  }
}

cat("\n✓ Done! All 'fia' scale names converted to '64kha'\n\n")

# Also rename directories if they contain "_fia_" in the name
cat("Checking for directories to rename...\n")
run_dirs <- list.dirs("runs", recursive = FALSE)
fia_dirs <- run_dirs[grepl("_fia_", run_dirs)]

if (length(fia_dirs) > 0) {
  cat("Found", length(fia_dirs), "directories with '_fia_' in name:\n")
  for (d in fia_dirs) {
    new_name <- gsub("_fia_", "_64kha_", d)
    cat("  ", basename(d), " → ", basename(new_name), "\n")
    
    # Actually rename
    if (!fs::dir_exists(new_name)) {
      fs::dir_copy(d, new_name)
      fs::dir_delete(d)
      cat("    ✓ Renamed\n")
    } else {
      cat("    ⚠ Target already exists, skipping\n")
    }
  }
} else {
  cat("  No directories need renaming\n")
}

cat("\n✓ Batch conversion complete!\n")
