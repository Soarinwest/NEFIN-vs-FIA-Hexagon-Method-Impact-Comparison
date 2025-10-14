# R/process_nefin_data.R
# Standardize NEFIN data to per-hectare basis for comparison with FIA
# Joins TREE_PLOT_DATA.csv with NEFIN_plots.csv for coordinates

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(sf); library(fs); library(yaml)
})

source("R/utils_spatial.R")

`%||%` <- function(a,b) if (!is.null(a)) a else b

#' Process NEFIN tree/plot data to FIA-comparable format
#' 
#' @param tree_csv Path to NEFIN TREE_PLOT_DATA.csv (biomass data)
#' @param plot_csv Path to NEFIN_plots.csv (location data)
#' @param out_dir Output directory for processed data
#' @param overwrite Overwrite existing output
#' 
#' @return Path to processed NEFIN data file
process_nefin_data <- function(tree_csv = "TREE_PLOT_DATA.csv",
                               plot_csv = "NEFIN_plots.csv",
                               out_dir = "data/processed",
                               overwrite = FALSE) {
  
  message("\n═══════════════════════════════════════════════════════════")
  message("NEFIN DATA STANDARDIZATION")
  message("═══════════════════════════════════════════════════════════\n")
  
  out_file <- fs::path(out_dir, "nefin_processed.csv")
  
  if (fs::file_exists(out_file) && !overwrite) {
    message("✓ Processed NEFIN data already exists: ", out_file)
    message("  Use overwrite=TRUE to regenerate")
    return(invisible(out_file))
  }
  
  fs::dir_create(out_dir, recurse = TRUE)
  
  # Read NEFIN tree data
  message("→ Reading NEFIN tree data: ", tree_csv)
  if (!fs::file_exists(tree_csv)) {
    stop("Tree data not found: ", tree_csv)
  }
  tree_data <- readr::read_csv(tree_csv, show_col_types = FALSE)
  message("  Rows: ", format(nrow(tree_data), big.mark = ","))
  message("  Columns: ", ncol(tree_data))
  
  # Read NEFIN plot location data
  message("\n→ Reading NEFIN plot locations: ", plot_csv)
  if (!fs::file_exists(plot_csv)) {
    stop("Plot location data not found: ", plot_csv)
  }
  plot_data <- readr::read_csv(plot_csv, show_col_types = FALSE)
  message("  Plots: ", format(nrow(plot_data), big.mark = ","))
  message("  Columns: ", ncol(plot_data))
  
  # Examine plot data structure
  message("\n→ Examining plot data columns...")
  plot_cols <- names(plot_data)
  message("  Available columns: ", paste(head(plot_cols, 10), collapse = ", "))
  
  # Look for coordinate columns (common variations)
  coord_candidates <- list(
    lat = c("lat", "latitude", "Latitude", "LAT", "y", "Y", "northing"),
    lon = c("lon", "long", "longitude", "Longitude", "LON", "LONG", "x", "X", "easting")
  )
  
  lat_col <- NULL
  lon_col <- NULL
  
  for (cand in coord_candidates$lat) {
    if (cand %in% plot_cols) {
      lat_col <- cand
      break
    }
  }
  
  for (cand in coord_candidates$lon) {
    if (cand %in% plot_cols) {
      lon_col <- cand
      break
    }
  }
  
  if (is.null(lat_col) || is.null(lon_col)) {
    message("  Available columns in NEFIN_plots.csv:")
    print(plot_cols)
    stop("Could not find lat/lon columns in plot data.\n",
         "  Expected names: lat, latitude, lon, long, longitude, etc.\n",
         "  Please specify correct column names or rename columns.")
  }
  
  message("  ✓ Found latitude column: ", lat_col)
  message("  ✓ Found longitude column: ", lon_col)
  
  # Look for plot ID column to join on
  plot_id_candidates <- c("_nefin_plotID", "plotID", "plot_id", "PLOTID", "PlotID")
  plot_id_col <- NULL
  
  for (cand in plot_id_candidates) {
    if (cand %in% names(tree_data) && cand %in% names(plot_data)) {
      plot_id_col <- cand
      break
    }
  }
  
  if (is.null(plot_id_col)) {
    # Try to find matching columns
    tree_id_cols <- names(tree_data)[grepl("plot|PLOT|Plot", names(tree_data), ignore.case = TRUE)]
    plot_id_cols <- names(plot_data)[grepl("plot|PLOT|Plot", names(plot_data), ignore.case = TRUE)]
    
    message("\n  Plot ID columns in tree data: ", paste(tree_id_cols, collapse = ", "))
    message("  Plot ID columns in plot data: ", paste(plot_id_cols, collapse = ", "))
    
    stop("Could not find matching plot ID column for joining.\n",
         "  Please ensure both files have a common plot identifier column.")
  }
  
  message("  ✓ Found plot ID column: ", plot_id_col)
  
  # Examine tree data structure
  message("\n→ Examining tree data structure...")
  key_cols <- c("treeSampleYear", plot_id_col, "AGB_kgPH", "BAPH", "TPH", "QMD", "STATUS")
  for (col in key_cols) {
    if (col %in% names(tree_data)) {
      message("    ✓ ", col)
    } else {
      message("    ⚠ ", col, " (MISSING - may affect processing)")
    }
  }
  
  # Join tree data with plot locations
  message("\n→ Joining tree data with plot locations...")
  
  nefin_joined <- tree_data |>
    dplyr::left_join(
      plot_data |> dplyr::select(!!rlang::sym(plot_id_col), 
                                 lat = !!rlang::sym(lat_col), 
                                 lon = !!rlang::sym(lon_col)),
      by = plot_id_col
    )
  
  # Check join success
  n_with_coords <- sum(!is.na(nefin_joined$lat) & !is.na(nefin_joined$lon))
  n_without_coords <- nrow(nefin_joined) - n_with_coords
  
  message("  Records with coordinates: ", format(n_with_coords, big.mark = ","), 
          " (", round(100 * n_with_coords / nrow(nefin_joined), 1), "%)")
  
  if (n_without_coords > 0) {
    message("  ⚠ Records without coordinates: ", format(n_without_coords, big.mark = ","))
    message("    These will be excluded from spatial analysis")
  }
  
  # Data validation
  message("\n→ Data validation...")
  
  # Check for NA values in critical columns
  na_checks <- list(
    Year = sum(is.na(nefin_joined$treeSampleYear)),
    PlotID = sum(is.na(nefin_joined[[plot_id_col]])),
    Biomass = sum(is.na(nefin_joined$AGB_kgPH)),
    Lat = sum(is.na(nefin_joined$lat)),
    Lon = sum(is.na(nefin_joined$lon))
  )
  
  for (check in names(na_checks)) {
    n_na <- na_checks[[check]]
    if (n_na > 0) {
      message("  ⚠ ", check, ": ", format(n_na, big.mark = ","), " missing values (", 
              round(100 * n_na / nrow(nefin_joined), 1), "%)")
    } else {
      message("  ✓ ", check, ": No missing values")
    }
  }
  
  # Year range
  year_range <- range(nefin_joined$treeSampleYear, na.rm = TRUE)
  message("  Year range: ", year_range[1], " - ", year_range[2])
  
  # Unique plots
  n_plots <- dplyr::n_distinct(nefin_joined[[plot_id_col]])
  message("  Unique plots: ", format(n_plots, big.mark = ","))
  
  # Convert NEFIN to FIA-comparable format
  message("\n→ Converting to FIA-comparable format...")
  
  # NEFIN appears to already have per-hectare values
  # Need to aggregate to plot level (currently may be tree-level records)
  
  nefin_processed <- nefin_joined |>
    dplyr::filter(
      !is.na(treeSampleYear),
      !is.na(!!rlang::sym(plot_id_col)),
      !is.na(lat), !is.na(lon),
      is.finite(lat), is.finite(lon)
    ) |>
    # Aggregate to plot level (in case data is tree-level)
    dplyr::group_by(
      treeSampleYear,
      !!rlang::sym(plot_id_col),
      lat, lon
    ) |>
    dplyr::summarise(
      # For per-hectare values that are already calculated, take mean
      # (if already plot-level, this will just keep the values)
      BAPH = mean(BAPH, na.rm = TRUE),  # Basal area per hectare
      TPH = mean(TPH, na.rm = TRUE),    # Trees per hectare
      QMD = mean(QMD, na.rm = TRUE),    # Quadratic mean diameter
      AGB_kgPH = mean(AGB_kgPH, na.rm = TRUE),  # Biomass kg per hectare
      
      # Count records (to detect if aggregating)
      n_records = dplyr::n(),
      .groups = "drop"
    ) |>
    # Convert kg/ha to Mg/ha (divide by 1000)
    dplyr::mutate(
      aglb_Mg_per_ha = AGB_kgPH / 1000,
      
      # Rename for consistency with FIA
      CN = as.character(!!rlang::sym(plot_id_col)),
      MEASYEAR = treeSampleYear,
      lat_public = lat,
      lon_public = lon,
      STATECD = "NEFIN",  # Or extract from data if available
      
      # Data source flag
      source = "NEFIN"
    ) |>
    # Validate coordinates
    dplyr::filter(
      lat_public >= -90, lat_public <= 90,
      lon_public >= -180, lon_public <= 180,
      is.finite(aglb_Mg_per_ha),
      aglb_Mg_per_ha >= 0  # Biomass can't be negative
    ) |>
    # Select final columns (matching FIA structure)
    dplyr::select(
      CN, STATECD, MEASYEAR, 
      lat_public, lon_public,
      aglb_Mg_per_ha,
      BAPH, TPH, QMD,
      n_records, source
    )
  
  # Summary statistics
  message("\n→ Processing summary...")
  message("  Input tree records: ", format(nrow(tree_data), big.mark = ","))
  message("  After join with locations: ", format(nrow(nefin_joined), big.mark = ","))
  message("  Output plots: ", format(nrow(nefin_processed), big.mark = ","))
  
  if (any(nefin_processed$n_records > 1)) {
    message("  ⚠ Some plots had multiple records (aggregated by mean)")
    message("    Max records per plot: ", max(nefin_processed$n_records))
  }
  
  # Biomass statistics
  biomass_stats <- nefin_processed |>
    dplyr::summarise(
      mean_aglb = mean(aglb_Mg_per_ha, na.rm = TRUE),
      median_aglb = median(aglb_Mg_per_ha, na.rm = TRUE),
      min_aglb = min(aglb_Mg_per_ha, na.rm = TRUE),
      max_aglb = max(aglb_Mg_per_ha, na.rm = TRUE),
      sd_aglb = sd(aglb_Mg_per_ha, na.rm = TRUE)
    )
  
  message("\n  Biomass statistics (Mg/ha):")
  message("    Mean: ", round(biomass_stats$mean_aglb, 2))
  message("    Median: ", round(biomass_stats$median_aglb, 2))
  message("    SD: ", round(biomass_stats$sd_aglb, 2))
  message("    Range: ", round(biomass_stats$min_aglb, 2), " - ", 
          round(biomass_stats$max_aglb, 2))
  
  # Coordinate range (sanity check)
  message("\n  Coordinate ranges:")
  message("    Latitude: ", round(min(nefin_processed$lat_public), 4), " to ",
          round(max(nefin_processed$lat_public), 4))
  message("    Longitude: ", round(min(nefin_processed$lon_public), 4), " to ",
          round(max(nefin_processed$lon_public), 4))
  
  # Write output
  message("\n→ Writing processed data...")
  readr::write_csv(nefin_processed, out_file)
  message("  ✓ Wrote: ", out_file)
  
  # Create summary report
  summary_file <- fs::path(out_dir, "nefin_processing_summary.txt")
  summary_lines <- c(
    "NEFIN Data Processing Summary",
    paste("Generated:", Sys.time()),
    "",
    paste("Tree data file:", tree_csv),
    paste("Plot location file:", plot_csv),
    paste("Output file:", out_file),
    "",
    "=== Data Overview ===",
    paste("Input tree records:", format(nrow(tree_data), big.mark = ",")),
    paste("Input plot locations:", format(nrow(plot_data), big.mark = ",")),
    paste("After joining:", format(nrow(nefin_joined), big.mark = ",")),
    paste("Output plots:", format(nrow(nefin_processed), big.mark = ",")),
    paste("Year range:", year_range[1], "-", year_range[2]),
    "",
    "=== Join Details ===",
    paste("Plot ID column:", plot_id_col),
    paste("Latitude column:", lat_col),
    paste("Longitude column:", lon_col),
    paste("Records with coordinates:", format(n_with_coords, big.mark = ",")),
    paste("Records without coordinates:", format(n_without_coords, big.mark = ",")),
    "",
    "=== Biomass Statistics (Mg/ha) ===",
    paste("Mean:", round(biomass_stats$mean_aglb, 2)),
    paste("Median:", round(biomass_stats$median_aglb, 2)),
    paste("SD:", round(biomass_stats$sd_aglb, 2)),
    paste("Min:", round(biomass_stats$min_aglb, 2)),
    paste("Max:", round(biomass_stats$max_aglb, 2)),
    "",
    "=== Column Mapping ===",
    "NEFIN → FIA equivalent:",
    paste0("  ", plot_id_col, " → CN (plot identifier)"),
    "  treeSampleYear → MEASYEAR (measurement year)",
    "  AGB_kgPH / 1000 → aglb_Mg_per_ha (biomass in Mg/ha)",
    paste0("  ", lat_col, " → lat_public"),
    paste0("  ", lon_col, " → lon_public"),
    "",
    "=== Notes ===",
    "- NEFIN data already in per-hectare units (BAPH, TPH, AGB_kgPH)",
    "- Converted kg/ha to Mg/ha by dividing by 1000",
    "- Aggregated by plot (if multiple tree records per plot)",
    "- Joined tree data with plot locations",
    "- Ready for direct comparison with FIA data"
  )
  
  writeLines(summary_lines, summary_file)
  message("  ✓ Wrote summary: ", summary_file)
  
  message("\n═══════════════════════════════════════════════════════════")
  message("NEFIN DATA READY")
  message("═══════════════════════════════════════════════════════════\n")
  message("Next steps:")
  message("  1. Assign NEFIN plots to hex grids:")
  message("     Rscript R/assign_nefin_to_hex.R")
  message("")
  message("  2. Compare FIA vs NEFIN:")
  message("     Rscript R/compare_fia_nefin.R")
  message("")
  
  invisible(list(
    output = out_file,
    summary = summary_file,
    n_plots = nrow(nefin_processed),
    biomass_stats = biomass_stats
  ))
}

# CLI
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  get_arg <- function(flag, default = NULL) {
    hit <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
  }
  
  process_nefin_data(
    tree_csv = get_arg("--tree", "TREE_PLOT_DATA.csv"),
    plot_csv = get_arg("--plots", "NEFIN_plots.csv"),
    out_dir = get_arg("--output", "data/processed"),
    overwrite = "--overwrite" %in% args
  )
}


message("\n═══════════════════════════════════════════════════════════")
message("NEFIN DATA STANDARDIZATION")
message("═══════════════════════════════════════════════════════════\n")

out_file <- fs::path(out_dir, "nefin_processed.csv")

if (fs::file_exists(out_file) && !overwrite) {
  message("✓ Processed NEFIN data already exists: ", out_file)
  message("  Use overwrite=TRUE to regenerate")
  return(invisible(out_file))
}

fs::dir_create(out_dir, recurse = TRUE)

# Read NEFIN data
message("→ Reading NEFIN data: ", nefin_csv)
nefin <- readr::read_csv(nefin_csv, show_col_types = FALSE)

message("  Rows: ", format(nrow(nefin), big.mark = ","))
message("  Columns: ", ncol(nefin))

# Examine column names
message("\n→ Examining data structure...")
message("  Key columns found:")
key_cols <- c("treeSampleYear", "_nefin_plotID", "AGB_kgPH", "BAPH", "TPH", 
              "QMD", "lat", "long", "STATUS")
for (col in key_cols) {
  if (col %in% names(nefin)) {
    message("    ✓ ", col)
  } else {
    message("    ✗ ", col, " (MISSING)")
  }
}

# Check for plot area information
# NEFIN typically stores: BAPH (basal area per hectare), TPH (trees per hectare)
# AGB_kgPH is already per hectare if "PH" = per hectare

message("\n→ Understanding NEFIN metrics...")
message("  BAPH: Basal Area per Hectare (m²/ha)")
message("  TPH: Trees per Hectare (stems/ha)")
message("  AGB_kgPH: Aboveground Biomass (kg/ha)")
message("  QMD: Quadratic Mean Diameter (cm)")

# Data validation
message("\n→ Data validation...")

# Check for NA values in critical columns
na_checks <- list(
  Year = sum(is.na(nefin$treeSampleYear)),
  PlotID = sum(is.na(nefin$`_nefin_plotID`)),
  Biomass = sum(is.na(nefin$AGB_kgPH)),
  Lat = sum(is.na(nefin$lat)),
  Lon = sum(is.na(nefin$long))
)

for (check in names(na_checks)) {
  n_na <- na_checks[[check]]
  if (n_na > 0) {
    message("  ⚠ ", check, ": ", n_na, " missing values (", 
            round(100 * n_na / nrow(nefin), 1), "%)")
  } else {
    message("  ✓ ", check, ": No missing values")
  }
}

# Year range
year_range <- range(nefin$treeSampleYear, na.rm = TRUE)
message("  Year range: ", year_range[1], " - ", year_range[2])

# Unique plots
n_plots <- dplyr::n_distinct(nefin$`_nefin_plotID`)
message("  Unique plots: ", format(n_plots, big.mark = ","))

# Convert NEFIN to FIA-comparable format
message("\n→ Converting to FIA-comparable format...")

# NEFIN appears to already have per-hectare values
# Need to aggregate to plot level (currently may be tree-level records)

nefin_processed <- nefin |>
  dplyr::filter(
    !is.na(treeSampleYear),
    !is.na(`_nefin_plotID`),
    !is.na(lat), !is.na(long),
    is.finite(lat), is.finite(long)
  ) |>
  # Aggregate to plot level (in case data is tree-level)
  dplyr::group_by(
    treeSampleYear,
    `_nefin_plotID`,
    lat, long,
    `_nefin_state`
  ) |>
  dplyr::summarise(
    # Sum per-hectare values if multiple records per plot
    # (if already plot-level, this will just keep the values)
    BAPH = mean(BAPH, na.rm = TRUE),  # Basal area per hectare
    TPH = mean(TPH, na.rm = TRUE),    # Trees per hectare
    QMD = mean(QMD, na.rm = TRUE),    # Quadratic mean diameter
    AGB_kgPH = mean(AGB_kgPH, na.rm = TRUE),  # Biomass kg per hectare
    
    # Count records (to detect if aggregating)
    n_records = dplyr::n(),
    .groups = "drop"
  ) |>
  # Convert kg/ha to Mg/ha (divide by 1000)
  dplyr::mutate(
    aglb_Mg_per_ha = AGB_kgPH / 1000,
    
    # Rename for consistency
    CN = `_nefin_plotID`,
    MEASYEAR = treeSampleYear,
    lat_public = lat,
    lon_public = long,
    STATECD = `_nefin_state`,
    
    # Data source flag
    source = "NEFIN"
  ) |>
  # Validate coordinates
  dplyr::filter(
    lat_public >= -90, lat_public <= 90,
    lon_public >= -180, lon_public <= 180,
    is.finite(aglb_Mg_per_ha)
  ) |>
  # Select final columns (matching FIA structure)
  dplyr::select(
    CN, STATECD, MEASYEAR, 
    lat_public, lon_public,
    aglb_Mg_per_ha,
    BAPH, TPH, QMD,
    n_records, source
  )

# Summary statistics
message("\n→ Processing summary...")
message("  Input rows: ", format(nrow(nefin), big.mark = ","))
message("  Output plots: ", format(nrow(nefin_processed), big.mark = ","))

if (any(nefin_processed$n_records > 1)) {
  message("  ⚠ Some plots had multiple records (aggregated by mean)")
  message("    Max records per plot: ", max(nefin_processed$n_records))
}

# Biomass statistics
biomass_stats <- nefin_processed |>
  dplyr::summarise(
    mean_aglb = mean(aglb_Mg_per_ha, na.rm = TRUE),
    median_aglb = median(aglb_Mg_per_ha, na.rm = TRUE),
    min_aglb = min(aglb_Mg_per_ha, na.rm = TRUE),
    max_aglb = max(aglb_Mg_per_ha, na.rm = TRUE),
    sd_aglb = sd(aglb_Mg_per_ha, na.rm = TRUE)
  )

message("\n  Biomass statistics (Mg/ha):")
message("    Mean: ", round(biomass_stats$mean_aglb, 2))
message("    Median: ", round(biomass_stats$median_aglb, 2))
message("    SD: ", round(biomass_stats$sd_aglb, 2))
message("    Range: ", round(biomass_stats$min_aglb, 2), " - ", 
        round(biomass_stats$max_aglb, 2))

# Write output
message("\n→ Writing processed data...")
readr::write_csv(nefin_processed, out_file)
message("  ✓ Wrote: ", out_file)

# Create summary report
summary_file <- fs::path(out_dir, "nefin_processing_summary.txt")
summary_lines <- c(
  "NEFIN Data Processing Summary",
  paste("Generated:", Sys.time()),
  "",
  paste("Input file:", nefin_csv),
  paste("Output file:", out_file),
  "",
  "=== Data Overview ===",
  paste("Input rows:", format(nrow(nefin), big.mark = ",")),
  paste("Output plots:", format(nrow(nefin_processed), big.mark = ",")),
  paste("Year range:", year_range[1], "-", year_range[2]),
  paste("States:", paste(unique(nefin_processed$STATECD), collapse = ", ")),
  "",
  "=== Biomass Statistics (Mg/ha) ===",
  paste("Mean:", round(biomass_stats$mean_aglb, 2)),
  paste("Median:", round(biomass_stats$median_aglb, 2)),
  paste("SD:", round(biomass_stats$sd_aglb, 2)),
  paste("Min:", round(biomass_stats$min_aglb, 2)),
  paste("Max:", round(biomass_stats$max_aglb, 2)),
  "",
  "=== Column Mapping ===",
  "NEFIN → FIA equivalent:",
  "  _nefin_plotID → CN (plot identifier)",
  "  treeSampleYear → MEASYEAR (measurement year)",
  "  AGB_kgPH / 1000 → aglb_Mg_per_ha (biomass in Mg/ha)",
  "  lat → lat_public",
  "  long → lon_public",
  "",
  "=== Notes ===",
  "- NEFIN data already in per-hectare units (BAPH, TPH, AGB_kgPH)",
  "- Converted kg/ha to Mg/ha by dividing by 1000",
  "- Aggregated by plot (if multiple tree records per plot)",
  "- Ready for direct comparison with FIA data"
)

writeLines(summary_lines, summary_file)
message("  ✓ Wrote summary: ", summary_file)

message("\n═══════════════════════════════════════════════════════════")
message("NEFIN DATA READY")
message("═══════════════════════════════════════════════════════════\n")
message("Next steps:")
message("  1. Assign NEFIN plots to hex grids:")
message("     Rscript R/assign_nefin_to_hex.R")
message("")
message("  2. Compare FIA vs NEFIN:")
message("     Rscript R/compare_fia_nefin.R")
message("")

invisible(list(
  output = out_file,
  summary = summary_file,
  n_plots = nrow(nefin_processed),
  biomass_stats = biomass_stats
))
}

# CLI
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  get_arg <- function(flag, default = NULL) {
    hit <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
  }
  
  process_nefin_data(
    nefin_csv = get_arg("--input", "TREE_PLOT_DATA.csv"),
    out_dir = get_arg("--output", "data/processed"),
    overwrite = "--overwrite" %in% args
  )
}