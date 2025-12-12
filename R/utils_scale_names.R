# R/utils_scale_names.R
# Utility functions for standardizing grid scale names
# Renames "fia" → "64kha" for proper ordering in figures

suppressPackageStartupMessages({
  library(dplyr)
})

# =============================================================================
# SCALE NAME STANDARDIZATION
# =============================================================================

#' Standardize grid scale names
#' Converts "fia" to "64kha" and ensures consistent naming
#' 
#' @param x Character vector of scale names
#' @return Character vector with standardized names
standardize_scale_name <- function(x) {
  x <- as.character(x)
  # Replace "fia" with "64kha" (case insensitive)
  x[tolower(x) == "fia"] <- "64kha"
  x
}

#' Standardize grid_scale column in a data frame
#' 
#' @param df Data frame with grid_scale column
#' @param scale_col Name of the scale column (default "grid_scale")
#' @return Data frame with standardized scale names
standardize_scale_df <- function(df, scale_col = "grid_scale") {
  if (scale_col %in% names(df)) {
    df[[scale_col]] <- standardize_scale_name(df[[scale_col]])
  }
  df
}

#' Get correct ordering of grid scales (smallest to largest)
#' 
#' @return Character vector of scale names in size order
get_scale_order <- function() {
  c("100ha", "500ha", "1kha", "5kha", "10kha", "50kha", "64kha", "100kha")
}

#' Convert grid_scale to ordered factor with correct size ordering
#' 
#' @param df Data frame with grid_scale column
#' @param scale_col Name of the scale column (default "grid_scale")
#' @return Data frame with grid_scale as ordered factor
order_scales <- function(df, scale_col = "grid_scale") {
  if (scale_col %in% names(df)) {
    # First standardize names
    df[[scale_col]] <- standardize_scale_name(df[[scale_col]])
    # Then convert to ordered factor
    scale_order <- get_scale_order()
    # Only include levels that exist in data
    existing_scales <- unique(df[[scale_col]])
    ordered_levels <- scale_order[scale_order %in% existing_scales]
    # Add any scales not in our predefined order (shouldn't happen, but safe)
    extra_scales <- setdiff(existing_scales, scale_order)
    if (length(extra_scales) > 0) {
      ordered_levels <- c(ordered_levels, sort(extra_scales))
    }
    df[[scale_col]] <- factor(df[[scale_col]], levels = ordered_levels, ordered = TRUE)
  }
  df
}

#' Get approximate area in hectares for each scale
#' 
#' @param scale_name Character scale name
#' @return Numeric area in hectares
scale_to_area_ha <- function(scale_name) {
  scale_name <- standardize_scale_name(scale_name)
  dplyr::case_when(
    scale_name == "100ha" ~ 100,
    scale_name == "500ha" ~ 500,
    scale_name == "1kha" ~ 1000,
    scale_name == "5kha" ~ 5000,
    scale_name == "10kha" ~ 10000,
    scale_name == "50kha" ~ 50000,
    scale_name == "64kha" ~ 64000,
    scale_name == "100kha" ~ 100000,
    TRUE ~ NA_real_
  )
}

#' Add area column to data frame based on grid_scale
#' 
#' @param df Data frame with grid_scale column
#' @param scale_col Name of the scale column
#' @param area_col Name for new area column
#' @return Data frame with added area column
add_scale_area <- function(df, scale_col = "grid_scale", area_col = "grid_area_ha") {
  if (scale_col %in% names(df)) {
    df[[area_col]] <- scale_to_area_ha(df[[scale_col]])
  }
  df
}

# =============================================================================
# FILE/PATH UTILITIES
# =============================================================================

#' Standardize scale name in file path
#' Converts paths like "runs/2025-12-11_aglb_fia_W5y" to use "64kha"
#' 
#' @param path File path string
#' @return Path with "fia" replaced by "64kha" in scale context
standardize_path_scale <- function(path) {
  # Replace _fia_ with _64kha_ (scale in directory names)
  path <- gsub("_fia_", "_64kha_", path)
  # Replace _fia/ with _64kha/ (end of directory name)
  path <- gsub("_fia/", "_64kha/", path)
  # Replace _fia$ with _64kha (end of string)
  path <- gsub("_fia$", "_64kha", path)
  path
}

# =============================================================================
# CSV I/O WITH STANDARDIZATION
# =============================================================================

#' Read CSV and standardize grid_scale names
#' 
#' @param path Path to CSV file
#' @param ... Additional arguments passed to readr::read_csv
#' @return Data frame with standardized scale names
read_csv_standardized <- function(path, ...) {
  df <- readr::read_csv(path, show_col_types = FALSE, ...)
  
  # Standardize common scale columns
  scale_cols <- c("grid_scale", "scale_name", "scale")
  for (col in scale_cols) {
    if (col %in% names(df)) {
      df[[col]] <- standardize_scale_name(df[[col]])
    }
  }
  
  df
}

#' Write CSV with standardized scale names
#' 
#' @param df Data frame to write
#' @param path Output path
#' @param ... Additional arguments passed to readr::write_csv
write_csv_standardized <- function(df, path, ...) {
  # Standardize before writing
  scale_cols <- c("grid_scale", "scale_name", "scale")
  for (col in scale_cols) {
    if (col %in% names(df)) {
      df[[col]] <- standardize_scale_name(df[[col]])
    }
  }
  
  readr::write_csv(df, path, ...)
}

# =============================================================================
# BATCH CONVERSION UTILITY
# =============================================================================

#' Convert all CSVs in a directory to use standardized scale names
#' 
#' @param dir Directory containing CSV files
#' @param recursive Search recursively?
#' @param dry_run If TRUE, only report what would be changed
#' @return Invisible list of modified files
batch_standardize_csvs <- function(dir, recursive = TRUE, dry_run = FALSE) {
  csv_files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE, recursive = recursive)
  
  modified <- character(0)
  
  for (f in csv_files) {
    df <- readr::read_csv(f, show_col_types = FALSE)
    
    # Check if any scale columns contain "fia"
    scale_cols <- c("grid_scale", "scale_name", "scale")
    needs_update <- FALSE
    
    for (col in scale_cols) {
      if (col %in% names(df)) {
        if (any(tolower(df[[col]]) == "fia", na.rm = TRUE)) {
          needs_update <- TRUE
          break
        }
      }
    }
    
    if (needs_update) {
      if (dry_run) {
        message("Would update: ", f)
      } else {
        # Standardize and rewrite
        for (col in scale_cols) {
          if (col %in% names(df)) {
            df[[col]] <- standardize_scale_name(df[[col]])
          }
        }
        readr::write_csv(df, f)
        message("Updated: ", f)
      }
      modified <- c(modified, f)
    }
  }
  
  if (length(modified) == 0) {
    message("No files needed updating")
  } else {
    message("\n", length(modified), " file(s) ", 
            if(dry_run) "would be " else "", "updated")
  }
  
  invisible(modified)
}

# =============================================================================
# GGPLOT SCALE HELPERS
# =============================================================================

#' Custom x-axis scale for grid scales with proper ordering
#' 
#' @param ... Additional arguments passed to scale_x_discrete
#' @return ggplot2 scale object
scale_x_grid_scale <- function(...) {
  ggplot2::scale_x_discrete(limits = get_scale_order(), ...)
}

#' Custom fill scale for grid scales
#' 
#' @param ... Additional arguments passed to scale_fill_viridis_d
#' @return ggplot2 scale object
scale_fill_grid_scale <- function(...) {
  ggplot2::scale_fill_viridis_d(option = "D", ...)
}

#' Custom color scale for grid scales  
#' 
#' @param ... Additional arguments passed to scale_color_viridis_d
#' @return ggplot2 scale object
scale_color_grid_scale <- function(...) {
  ggplot2::scale_color_viridis_d(option = "D", ...)
}
