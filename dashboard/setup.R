#!/usr/bin/env Rscript
# =============================================================================
# setup.R - Install dependencies and verify dashboard setup
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  NEFIN vs FIA Dashboard Setup                                        ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Required packages
required_packages <- c(
  # Shiny ecosystem
  "shiny",
  "bslib",
  "shinycssloaders",
  "htmltools",
  
  # Mapping
  "leaflet",
  "leaflet.extras",
  "sf",
  "terra",
  
  # Plotting
  "plotly",
  "ggplot2",
  "viridis",
  "scales",
  "patchwork",
  
  # Data manipulation
  "dplyr",
  "tidyr",
  "readr",
  
  # Tables
  "DT"
)

cat("Checking required packages...\n\n")

# Check and install missing packages
missing <- c()
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing <- c(missing, pkg)
    cat("  ✗", pkg, "- NOT INSTALLED\n")
  } else {
    cat("  ✓", pkg, "\n")
  }
}

if (length(missing) > 0) {
  cat("\n")
  cat("Missing packages:", paste(missing, collapse = ", "), "\n")
  cat("\nInstall missing packages? (y/n): ")
  
  response <- readline()
  
  if (tolower(response) == "y") {
    cat("\nInstalling packages...\n")
    install.packages(missing, repos = "https://cloud.r-project.org")
    cat("\n✓ Installation complete!\n")
  } else {
    cat("\nSkipping installation. Install manually with:\n")
    cat('  install.packages(c("', paste(missing, collapse = '", "'), '"))\n')
  }
} else {
  cat("\n✓ All packages installed!\n")
}

# Verify directory structure
cat("\nVerifying dashboard structure...\n")

required_files <- c(
  "app.R",
  "R/utils.R",
  "R/mod_overview.R",
  "R/mod_data.R",
  "R/mod_fuzzing.R",
  "R/mod_models.R",
  "R/mod_maps.R",
  "www/styles.css"
)

all_present <- TRUE
for (f in required_files) {
  if (file.exists(f)) {
    cat("  ✓", f, "\n")
  } else {
    cat("  ✗", f, "- MISSING\n")
    all_present <- FALSE
  }
}

if (!all_present) {
  cat("\n⚠ Some files are missing. Dashboard may not work correctly.\n")
} else {
  cat("\n✓ All dashboard files present!\n")
}

# Check for data files
cat("\nChecking for data files...\n")

data_paths <- list(
  fia = "data/processed/fia_complete.csv",
  nefin = "data/processed/nefin_processed.csv",
  states = "data/boundaries/states_5070.geojson",
  hex = "data/hex/hex_grid.geojson"
)

has_data <- FALSE
for (name in names(data_paths)) {
  path <- data_paths[[name]]
  if (file.exists(path)) {
    cat("  ✓", name, ":", path, "\n")
    has_data <- TRUE
  } else {
    cat("  ⚠", name, ":", path, "- not found (will use demo data)\n")
  }
}

if (!has_data) {
  cat("\n")
  cat("No data files found. The dashboard will run in DEMO MODE with\n")
  cat("synthetic data. To use real data, update paths in app.R config.\n")
}

# Summary
cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("SETUP COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("\n")
cat("To run the dashboard:\n")
cat("  shiny::runApp('.')\n")
cat("\n")
cat("Or from command line:\n")
cat("  Rscript -e \"shiny::runApp('.')\"\n")
cat("\n")
