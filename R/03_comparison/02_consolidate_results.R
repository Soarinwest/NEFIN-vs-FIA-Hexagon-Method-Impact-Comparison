#!/usr/bin/env Rscript
# R/consolidate_nefin_comparison.R
# Combines all per-scale NEFIN comparison results into consolidated directory

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(fs)
})

# Source scale utilities if available
if (file.exists("R/utils_scale_names.R")) {
  source("R/utils_scale_names.R")
} else {
  standardize_scale_name <- function(x) {
    x <- as.character(x)
    x[tolower(x) == "fia"] <- "64kha"
    x
  }
}

consolidate_nefin_results <- function(runs_dir = "runs", 
                                       consolidated_dir = NULL) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Consolidating NEFIN Comparison Results                  ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  

  # Find consolidated directory if not specified
  if (is.null(consolidated_dir)) {
    all_dirs <- list.dirs(runs_dir, recursive = FALSE)
    consol_dirs <- all_dirs[grepl("^.*consolidated_", all_dirs)]
    
    if (length(consol_dirs) == 0) {
      stop("No consolidated directory found. Run --viz first to create one.")
    }
    
    # Pick most recent
    consol_dirs <- consol_dirs[order(file.mtime(consol_dirs), decreasing = TRUE)]
    consolidated_dir <- consol_dirs[1]
  }
  
  cat("Target consolidated directory:\n")
  cat("  ", consolidated_dir, "\n\n")
  
  # Find all NEFIN comparison directories
  all_dirs <- list.dirs(runs_dir, recursive = FALSE)
  nefin_dirs <- all_dirs[grepl("fia_nefin_comparison_", all_dirs)]
  
  if (length(nefin_dirs) == 0) {
    stop("No NEFIN comparison directories found. Run --compare-nefin first.")
  }
  
  cat("Found", length(nefin_dirs), "NEFIN comparison directories:\n")
  for (d in nefin_dirs) {
    cat("  -", basename(d), "\n")
  }
  cat("\n")
  
  # Collect all summary_5y CSVs
  all_data <- list()
  
  for (d in nefin_dirs) {
    # Extract scale name from directory
    scale_name <- gsub(".*fia_nefin_comparison_", "", basename(d))
    
    # Find summary_5y CSV
    csv_file <- fs::path(d, paste0("summary_5y_", scale_name, ".csv"))
    
    if (fs::file_exists(csv_file)) {
      cat("  Loading:", basename(csv_file), "\n")
      
      df <- readr::read_csv(csv_file, show_col_types = FALSE)
      
      # Ensure grid_scale column exists and is standardized
      if (!"grid_scale" %in% names(df)) {
        df$grid_scale <- scale_name
      }
      df$grid_scale <- standardize_scale_name(df$grid_scale)
      
      # Also standardize scale_name column if it exists
      if ("scale_name" %in% names(df)) {
        df$scale_name <- standardize_scale_name(df$scale_name)
      }
      
      all_data[[scale_name]] <- df
    } else {
      cat("  ⚠ Not found:", csv_file, "\n")
    }
  }
  
  if (length(all_data) == 0) {
    stop("No NEFIN comparison data files found.")
  }
  
  # Combine all scales
  combined <- dplyr::bind_rows(all_data)
  
  cat("\n")
  cat("Combined data:\n")
  cat("  Total rows:", nrow(combined), "\n")
  cat("  Scales:", paste(unique(combined$grid_scale), collapse = ", "), "\n")
  
  # Standardize column names for compatibility with advanced_analysis.R
  # The comparison file expects these columns:
  # hex_id, grid_scale, year_label, mean_fia, mean_nefin, diff, 
  # se_fia, se_nefin, n_plots_fia, n_plots_nefin, has_both
  
  combined <- combined %>%
    mutate(
      # Create coverage flags
      has_fia = !is.na(fia_only_mean),
      has_nefin = !is.na(nefin_mean_5y),
      has_both = has_fia & has_nefin,
      
      # Rename columns to match expected format
      mean_fia = fia_only_mean,
      mean_nefin = nefin_mean_5y,
      se_fia = fia_only_se,
      se_nefin = nefin_se_5y,
      n_plots_fia = fia_only_n,
      n_plots_nefin = nefin_n_5y,
      diff = nefin_mean_5y - fia_only_mean,
      
      # Percent difference (relative to FIA)
      pct_diff = dplyr::case_when(
        is.na(fia_only_mean) | is.na(nefin_mean_5y) ~ NA_real_,
        fia_only_mean == 0 ~ NA_real_,
        TRUE ~ 100 * (nefin_mean_5y - fia_only_mean) / fia_only_mean
      ),
      
      # Add year_label if missing (use middle of window)
      year_label = if ("year_label" %in% names(.)) year_label else 2022
    )
  
  # Write to consolidated directory
  output_file <- fs::path(consolidated_dir, "fia_nefin_comparison_all_scales.csv")
  readr::write_csv(combined, output_file)
  
  cat("\n")
  cat("✓ Wrote consolidated NEFIN comparison:\n")
  cat("  ", output_file, "\n")
  cat("  Rows:", nrow(combined), "\n")
  cat("  Columns:", ncol(combined), "\n")
  
  # Summary stats
  cat("\n")
  cat("Summary by scale:\n")
  summary_stats <- combined %>%
    group_by(grid_scale) %>%
    summarise(
      n_hexes = n(),
      n_both = sum(has_both, na.rm = TRUE),
      pct_both = round(100 * n_both / n_hexes, 1),
      .groups = "drop"
    )
  
  print(as.data.frame(summary_stats), row.names = FALSE)
  
  cat("\n✓ Consolidation complete!\n\n")
  
  invisible(output_file)
}

# Run if called directly
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  consolidate_nefin_results()
}
