#!/usr/bin/env Rscript
# =============================================================================
# enrich_fia_complete.R
# Joins climate and NDVI covariates into fia_complete.csv for dashboard use
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  Enriching fia_complete.csv with Covariates                         ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# File paths
fia_path <- "data/processed/fia_complete.csv"
climate_path <- "data/processed/climate_at_plots/fia_climate.csv"
ndvi_path <- "data/processed/ndvi_at_plots/fia_ndvi.csv"
backup_path <- "data/processed/fia_complete_backup.csv"

# Check files exist
if (!file.exists(fia_path)) {
  stop("ERROR: fia_complete.csv not found at: ", fia_path)
}

# Load FIA data
cat("Loading fia_complete.csv...\n")
fia <- read_csv(fia_path, show_col_types = FALSE)
cat("  Loaded", format(nrow(fia), big.mark = ","), "rows,", ncol(fia), "columns\n")
cat("  Existing columns:", paste(names(fia), collapse = ", "), "\n\n")

# Check if covariates already present
existing_covars <- intersect(c("ndvi_modis", "tmean", "ppt"), names(fia))
if (length(existing_covars) == 3) {
  cat("✓ All covariates already present! No action needed.\n\n")
  quit(save = "no", status = 0)
}

# Create backup
cat("Creating backup at:", backup_path, "\n")
write_csv(fia, backup_path)

# Find the join column (CN is standard)
join_col <- if ("CN" %in% names(fia)) "CN" else {
  # Try to find an alternative
  cn_like <- names(fia)[grepl("^CN$|^cn$|PLOT_CN", names(fia), ignore.case = TRUE)][1]
  if (is.na(cn_like)) stop("Cannot find CN column for joining")
  cn_like
}
cat("Using join column:", join_col, "\n\n")

# Join climate data
if (file.exists(climate_path)) {
  cat("Loading climate data...\n")
  climate <- read_csv(climate_path, show_col_types = FALSE)
  cat("  Loaded", nrow(climate), "rows\n")
  cat("  Columns:", paste(names(climate), collapse = ", "), "\n")
  
  # Find matching join column in climate
  climate_join <- if (join_col %in% names(climate)) join_col else {
    names(climate)[grepl("^CN$|^cn$|PLOT_CN", names(climate), ignore.case = TRUE)][1]
  }
  
  if (!is.na(climate_join)) {
    # Rename if needed
    if (climate_join != join_col) {
      names(climate)[names(climate) == climate_join] <- join_col
    }
    
    # Select only the columns we need (avoid duplicates)
    climate_cols <- c(join_col, intersect(c("tmean", "ppt", "tmax", "tmin"), names(climate)))
    climate <- climate %>% select(all_of(climate_cols))
    
    # Remove duplicates
    climate <- climate %>% distinct(across(all_of(join_col)), .keep_all = TRUE)
    
    n_before <- nrow(fia)
    fia <- fia %>% left_join(climate, by = join_col)
    cat("  Joined climate:", sum(!is.na(fia$tmean)), "plots with tmean\n")
  } else {
    cat("  ⚠ Could not find matching join column in climate data\n")
  }
} else {
  cat("⚠ Climate file not found:", climate_path, "\n")
}

# Join NDVI data
if (file.exists(ndvi_path)) {
  cat("\nLoading NDVI data...\n")
  ndvi <- read_csv(ndvi_path, show_col_types = FALSE)
  cat("  Loaded", nrow(ndvi), "rows\n")
  cat("  Columns:", paste(names(ndvi), collapse = ", "), "\n")
  
  # Find matching join column
  ndvi_join <- if (join_col %in% names(ndvi)) join_col else {
    names(ndvi)[grepl("^CN$|^cn$|PLOT_CN", names(ndvi), ignore.case = TRUE)][1]
  }
  
  if (!is.na(ndvi_join)) {
    if (ndvi_join != join_col) {
      names(ndvi)[names(ndvi) == ndvi_join] <- join_col
    }
    
    # Select NDVI columns (look for various naming conventions)
    ndvi_col_options <- c("ndvi_modis", "ndvi_s2", "NDVI", "ndvi", "modis_ndvi", "s2_ndvi")
    available_ndvi <- intersect(ndvi_col_options, names(ndvi))
    
    if (length(available_ndvi) > 0) {
      ndvi_cols <- c(join_col, available_ndvi)
      ndvi <- ndvi %>% select(all_of(ndvi_cols))
      
      # Standardize column name if needed
      if ("NDVI" %in% names(ndvi) && !"ndvi_modis" %in% names(ndvi)) {
        names(ndvi)[names(ndvi) == "NDVI"] <- "ndvi_modis"
      }
      if ("ndvi" %in% names(ndvi) && !"ndvi_modis" %in% names(ndvi)) {
        names(ndvi)[names(ndvi) == "ndvi"] <- "ndvi_modis"
      }
      if ("modis_ndvi" %in% names(ndvi) && !"ndvi_modis" %in% names(ndvi)) {
        names(ndvi)[names(ndvi) == "modis_ndvi"] <- "ndvi_modis"
      }
      
      # Remove duplicates
      ndvi <- ndvi %>% distinct(across(all_of(join_col)), .keep_all = TRUE)
      
      fia <- fia %>% left_join(ndvi, by = join_col)
      
      if ("ndvi_modis" %in% names(fia)) {
        cat("  Joined NDVI:", sum(!is.na(fia$ndvi_modis)), "plots with ndvi_modis\n")
      }
    } else {
      cat("  ⚠ No recognized NDVI columns found\n")
    }
  } else {
    cat("  ⚠ Could not find matching join column in NDVI data\n")
  }
} else {
  cat("⚠ NDVI file not found:", ndvi_path, "\n")
}

# Summary
cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("RESULT\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

final_covars <- intersect(c("ndvi_modis", "ndvi_s2", "tmean", "ppt"), names(fia))
cat("Final columns:", ncol(fia), "\n")
cat("Covariates added:", paste(final_covars, collapse = ", "), "\n")

# Coverage stats
if ("tmean" %in% names(fia)) {
  cat("  tmean coverage:", round(100 * sum(!is.na(fia$tmean)) / nrow(fia), 1), "%\n")
}
if ("ppt" %in% names(fia)) {
  cat("  ppt coverage:", round(100 * sum(!is.na(fia$ppt)) / nrow(fia), 1), "%\n")
}
if ("ndvi_modis" %in% names(fia)) {
  cat("  ndvi_modis coverage:", round(100 * sum(!is.na(fia$ndvi_modis)) / nrow(fia), 1), "%\n")
}

# Save enriched file
cat("\nSaving enriched fia_complete.csv...\n")
write_csv(fia, fia_path)
cat("✓ Done! Original backed up to:", backup_path, "\n\n")
