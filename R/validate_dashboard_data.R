#!/usr/bin/env Rscript
# =============================================================================
# validate_dashboard_data.R
# Run this script to verify all data files are ready for the dashboard
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(sf)
})

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  NEFIN vs FIA Dashboard Data Validation                             ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Track issues
issues <- character()
warnings <- character()

# =============================================================================
# 1. Core Plot Data
# =============================================================================
cat("1. CORE PLOT DATA\n")
cat("─────────────────────────────────────────────────────────────────────\n")

# FIA
fia_path <- "data/processed/fia_complete.csv"
if (file.exists(fia_path)) {
  fia <- read_csv(fia_path, show_col_types = FALSE)
  cat("  ✓ fia_complete.csv:", format(nrow(fia), big.mark = ","), "rows\n")
  
  # Check required columns
  required_fia <- c("CN", "aglb_Mg_per_ha")
  lon_cols <- names(fia)[grepl("^lon|longitude", names(fia), ignore.case = TRUE)]
  lat_cols <- names(fia)[grepl("^lat|latitude", names(fia), ignore.case = TRUE)]
  
  if (length(lon_cols) > 0) {
    cat("    - Longitude column:", lon_cols[1], "\n")
  } else {
    issues <- c(issues, "FIA: No longitude column found")
    cat("    ✗ No longitude column found!\n")
  }
  
  if (length(lat_cols) > 0) {
    cat("    - Latitude column:", lat_cols[1], "\n")
  } else {
    issues <- c(issues, "FIA: No latitude column found")
    cat("    ✗ No latitude column found!\n")
  }
  
  if ("aglb_Mg_per_ha" %in% names(fia)) {
    cat("    - Biomass: ✓ (mean:", round(mean(fia$aglb_Mg_per_ha, na.rm = TRUE), 1), "Mg/ha)\n")
  } else {
    issues <- c(issues, "FIA: Missing aglb_Mg_per_ha column")
    cat("    ✗ aglb_Mg_per_ha missing!\n")
  }
  
  if ("MEASYEAR" %in% names(fia)) {
    cat("    - MEASYEAR: ✓ (", min(fia$MEASYEAR, na.rm = TRUE), "-", 
        max(fia$MEASYEAR, na.rm = TRUE), ")\n")
  } else {
    warnings <- c(warnings, "FIA: Missing MEASYEAR (temporal plots won't work)")
    cat("    ⚠ MEASYEAR missing (temporal plots disabled)\n")
  }
  
  # Check for covariates
  covars <- c("ndvi_modis", "tmean", "ppt")
  present_covars <- covars[covars %in% names(fia)]
  missing_covars <- setdiff(covars, names(fia))
  
  if (length(present_covars) > 0) {
    cat("    - Covariates present:", paste(present_covars, collapse = ", "), "\n")
  }
  if (length(missing_covars) > 0) {
    warnings <- c(warnings, paste("FIA: Missing covariates:", paste(missing_covars, collapse = ", ")))
    cat("    ⚠ Missing covariates:", paste(missing_covars, collapse = ", "), "\n")
  }
  
} else {
  issues <- c(issues, "FIA: fia_complete.csv not found")
  cat("  ✗ fia_complete.csv NOT FOUND\n")
}

# NEFIN
nefin_path <- "data/processed/nefin_processed.csv"
if (file.exists(nefin_path)) {
  nefin <- read_csv(nefin_path, show_col_types = FALSE)
  cat("  ✓ nefin_processed.csv:", format(nrow(nefin), big.mark = ","), "rows\n")
  
  lon_cols <- names(nefin)[grepl("^lon|longitude", names(nefin), ignore.case = TRUE)]
  lat_cols <- names(nefin)[grepl("^lat|latitude", names(nefin), ignore.case = TRUE)]
  
  if (length(lon_cols) > 0) {
    cat("    - Longitude column:", lon_cols[1], "\n")
  } else {
    issues <- c(issues, "NEFIN: No longitude column found")
    cat("    ✗ No longitude column found!\n")
  }
  
  if (length(lat_cols) > 0) {
    cat("    - Latitude column:", lat_cols[1], "\n")
  } else {
    issues <- c(issues, "NEFIN: No latitude column found")
    cat("    ✗ No latitude column found!\n")
  }
  
  if ("aglb_Mg_per_ha" %in% names(nefin)) {
    cat("    - Biomass: ✓ (mean:", round(mean(nefin$aglb_Mg_per_ha, na.rm = TRUE), 1), "Mg/ha)\n")
  } else {
    issues <- c(issues, "NEFIN: Missing aglb_Mg_per_ha column")
  }
  
} else {
  issues <- c(issues, "NEFIN: nefin_processed.csv not found")
  cat("  ✗ nefin_processed.csv NOT FOUND\n")
}

cat("\n")

# =============================================================================
# 2. Hex Grids and Results
# =============================================================================
cat("2. HEX GRIDS AND RESULTS\n")
cat("─────────────────────────────────────────────────────────────────────\n")

scales <- c("100ha", "500ha", "1kha", "5kha", "10kha", "50kha", "64kha", "100kha")

# Handle 64kha special case (file is named hex_grid_6kac.geojson)
grid_files <- c(
  "100ha" = "data/hex/hex_grid_100ha.geojson",
  "500ha" = "data/hex/hex_grid_500ha.geojson",
  "1kha" = "data/hex/hex_grid_1kha.geojson",
  "5kha" = "data/hex/hex_grid_5kha.geojson",
  "10kha" = "data/hex/hex_grid_10kha.geojson",
  "50kha" = "data/hex/hex_grid_50kha.geojson",
  "64kha" = "data/hex/hex_grid_6kac.geojson",  # Note: different name
  "100kha" = "data/hex/hex_grid_100kha.geojson"
)

result_pattern <- "runs/2025-12-15_aglb_%s_W5y/hex_aglb_results.csv"

for (scale in scales) {
  grid_path <- grid_files[scale]
  result_path <- sprintf(result_pattern, scale)
  
  grid_ok <- file.exists(grid_path)
  result_ok <- file.exists(result_path)
  
  status <- if (grid_ok && result_ok) "✓" else if (!grid_ok && !result_ok) "✗" else "⚠"
  
  cat(sprintf("  %s %s: Grid=%s Results=%s\n", 
              status, scale,
              if (grid_ok) "✓" else "✗",
              if (result_ok) "✓" else "✗"))
  
  if (!grid_ok) issues <- c(issues, paste("Hex grid missing:", grid_path))
  if (!result_ok) issues <- c(issues, paste("Hex results missing:", result_path))
  
  # Check hex_id column in both
  if (grid_ok) {
    hex_sf <- st_read(grid_path, quiet = TRUE)
    hex_id_cols <- intersect(names(hex_sf), c("hex_id", "HEX_ID", "id", "ID", "hexid"))
    if (length(hex_id_cols) == 0) {
      warnings <- c(warnings, paste(scale, "grid: No hex_id column"))
    }
  }
}

cat("\n")

# =============================================================================
# 3. Model Comparison Results
# =============================================================================
cat("3. MODEL COMPARISON RESULTS\n")
cat("─────────────────────────────────────────────────────────────────────\n")

model_comp_path <- "runs/spatial_model_comparison/model_comparison_summary.csv"
holdout_path <- "runs/spatial_model_comparison/holdout_prediction_results.csv"

if (file.exists(model_comp_path)) {
  model_comp <- read_csv(model_comp_path, show_col_types = FALSE)
  cat("  ✓ model_comparison_summary.csv:", nrow(model_comp), "rows\n")
  cat("    Columns:", paste(names(model_comp), collapse = ", "), "\n")
} else {
  warnings <- c(warnings, "Model comparison: model_comparison_summary.csv missing")
  cat("  ⚠ model_comparison_summary.csv NOT FOUND (will use demo data)\n")
}

if (file.exists(holdout_path)) {
  holdout <- read_csv(holdout_path, show_col_types = FALSE)
  cat("  ✓ holdout_prediction_results.csv:", nrow(holdout), "rows\n")
} else {
  warnings <- c(warnings, "Holdout: holdout_prediction_results.csv missing")
  cat("  ⚠ holdout_prediction_results.csv NOT FOUND\n")
}

cat("\n")

# =============================================================================
# 4. Fuzzing Analysis Results
# =============================================================================
cat("4. FUZZING ANALYSIS RESULTS\n")
cat("─────────────────────────────────────────────────────────────────────\n")

covar_unc_path <- "runs/fuzzing_effect_analysis/covariate_uncertainty_summary.csv"
pred_unc_path <- "runs/fuzzing_effect_analysis/prediction_uncertainty_by_plot.csv"

if (file.exists(covar_unc_path)) {
  covar_unc <- read_csv(covar_unc_path, show_col_types = FALSE)
  cat("  ✓ covariate_uncertainty_summary.csv:", nrow(covar_unc), "rows\n")
  cat("    Columns:", paste(names(covar_unc), collapse = ", "), "\n")
} else {
  warnings <- c(warnings, "Fuzzing: covariate_uncertainty_summary.csv missing")
  cat("  ⚠ covariate_uncertainty_summary.csv NOT FOUND\n")
}

if (file.exists(pred_unc_path)) {
  pred_unc <- read_csv(pred_unc_path, show_col_types = FALSE)
  cat("  ✓ prediction_uncertainty_by_plot.csv:", nrow(pred_unc), "rows\n")
} else {
  warnings <- c(warnings, "Fuzzing: prediction_uncertainty_by_plot.csv missing")
  cat("  ⚠ prediction_uncertainty_by_plot.csv NOT FOUND\n")
}

cat("\n")

# =============================================================================
# 5. Prediction Rasters
# =============================================================================
cat("5. PREDICTION RASTERS\n")
cat("─────────────────────────────────────────────────────────────────────\n")

rasters <- c(
  "modis_biomass_fia.tif" = "runs/prediction_maps/modis_biomass_fia.tif",
  "modis_biomass_nefin.tif" = "runs/prediction_maps/modis_biomass_nefin.tif",
  "modis_prediction_difference.tif" = "runs/prediction_maps/modis_prediction_difference.tif"
)

for (name in names(rasters)) {
  path <- rasters[name]
  if (file.exists(path)) {
    cat("  ✓", name, "\n")
  } else {
    warnings <- c(warnings, paste("Raster:", name, "missing"))
    cat("  ⚠", name, "NOT FOUND\n")
  }
}

cat("\n")

# =============================================================================
# 6. State Boundaries
# =============================================================================
cat("6. STATE BOUNDARIES\n")
cat("─────────────────────────────────────────────────────────────────────\n")

states_path <- "data/boundaries/states_5070.geojson"
if (file.exists(states_path)) {
  states <- st_read(states_path, quiet = TRUE)
  cat("  ✓ states_5070.geojson:", nrow(states), "states\n")
} else {
  warnings <- c(warnings, "States: states_5070.geojson missing")
  cat("  ⚠ states_5070.geojson NOT FOUND\n")
}

cat("\n")

# =============================================================================
# Summary
# =============================================================================
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("VALIDATION SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

if (length(issues) == 0 && length(warnings) == 0) {
  cat("\n✅ ALL CHECKS PASSED - Dashboard is ready to launch!\n")
} else {
  if (length(issues) > 0) {
    cat("\n❌ CRITICAL ISSUES (", length(issues), "):\n")
    for (issue in issues) {
      cat("  •", issue, "\n")
    }
  }
  
  if (length(warnings) > 0) {
    cat("\n⚠️  WARNINGS (", length(warnings), "):\n")
    for (warning in warnings) {
      cat("  •", warning, "\n")
    }
  }
}

cat("\n")
cat("To run the dashboard:\n")
cat("  cd dashboard\n")
cat("  Rscript -e \"shiny::runApp('.')\"\n")
cat("\n")
