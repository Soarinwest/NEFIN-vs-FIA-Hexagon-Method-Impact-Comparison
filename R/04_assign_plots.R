
suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(fs); library(yaml)
})

source("R/utils_spatial.R")
source("R/utils_metrics.R")

#' Stage 2: Assign plots to hexes and cache
#' This runs ONCE and stores the plot→hex mapping with original coords
stage2_assign_plots <- function(project_dir = ".",
                                hex_path = "data/hex/hex_grid.geojson",
                                hex_layer = NULL,
                                overwrite = FALSE) {
  
  # Paths
  fia_root <- if (fs::dir_exists(fs::path(project_dir, "data", "interim", "fia_region"))) {
    fs::path(project_dir, "data", "interim", "fia_region")
  } else {
    fs::path(project_dir, "data", "interim", "fia")
  }
  
  out_dir <- fs::path(project_dir, "data", "processed")
  out_file <- fs::path(out_dir, "plot_hex_assignments.csv")
  
  # Check if already exists
  if (fs::file_exists(out_file) && !overwrite) {
    message("✓ Plot assignments already exist: ", out_file)
    message("  Use overwrite=TRUE to regenerate")
    return(invisible(out_file))
  }
  
  fs::dir_create(out_dir, recurse = TRUE)
  
  message("→ Reading FIA PLOT data from: ", fia_root)
  plot_csv <- fs::path(fia_root, "plot.csv")
  tree_csv <- fs::path(fia_root, "tree.csv")
  
  if (!fs::file_exists(plot_csv)) stop("PLOT not found: ", plot_csv)
  if (!fs::file_exists(tree_csv)) stop("TREE not found: ", tree_csv)
  
  pl <- suppressMessages(readr::read_csv(plot_csv, guess_max = 1e6, show_col_types = FALSE))
  tr <- suppressMessages(readr::read_csv(tree_csv, guess_max = 1e6, show_col_types = FALSE))
  
  message("  PLOT rows: ", nrow(pl))
  message("  TREE rows: ", nrow(tr))
  
  # Normalize coordinates
  pl <- normalize_plot_coords(pl)
  
  # Store ORIGINAL coordinates explicitly
  pl <- pl |>
    dplyr::mutate(
      lat_original = lat_public,
      lon_original = lon_public
    )
  
  # Assign to hex (this uses the spatial join logic you already have)
  message("→ Assigning plots to hexes...")
  plot_hex <- assign_plots_to_hex(pl, hex_path, hex_layer)
  
  # Verify we kept original coords
  if (!all(c("lat_original", "lon_original") %in% names(plot_hex))) {
    stop("Original coordinates not preserved!")
  }
  
  # Select essential columns
  assignments <- plot_hex |>
    dplyr::select(
      CN, STATECD, MEASYEAR, UNITCD, COUNTYCD, PLOT,
      lat_original, lon_original,  # ORIGINAL FIA-provided fuzzy coords
      hex_id
    )
  
  # Add plot count summary
  summary_stats <- assignments |>
    dplyr::summarise(
      total_plots = dplyr::n(),
      plots_with_hex = sum(!is.na(hex_id)),
      plots_missing_hex = sum(is.na(hex_id)),
      unique_hexes = dplyr::n_distinct(hex_id, na.rm = TRUE),
      year_range = paste(range(MEASYEAR, na.rm = TRUE), collapse = "-")
    )
  
  message("\n=== Assignment Summary ===")
  message("  Total plots:        ", summary_stats$total_plots)
  message("  Assigned to hex:    ", summary_stats$plots_with_hex, 
          " (", round(100 * summary_stats$plots_with_hex / summary_stats$total_plots, 1), "%)")
  message("  Missing hex:        ", summary_stats$plots_missing_hex)
  message("  Unique hexes used:  ", summary_stats$unique_hexes)
  message("  Year range:         ", summary_stats$year_range)
  
  # Write
  readr::write_csv(assignments, out_file)
  message("\n✓ Wrote: ", out_file)
  
  # Also write summary
  summary_file <- fs::path(out_dir, "plot_hex_summary.txt")
  writeLines(c(
    "Plot-to-Hex Assignment Summary",
    paste("Generated:", Sys.time()),
    "",
    paste("Total plots:", summary_stats$total_plots),
    paste("Assigned to hex:", summary_stats$plots_with_hex),
    paste("Missing hex:", summary_stats$plots_missing_hex),
    paste("Unique hexes:", summary_stats$unique_hexes),
    paste("Year range:", summary_stats$year_range),
    "",
    "Files:",
    paste("  Assignments:", out_file),
    paste("  Hex grid:", hex_path)
  ), summary_file)
  
  invisible(list(
    assignments = out_file,
    summary = summary_file,
    stats = summary_stats
  ))
}

# CLI
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  overwrite <- "--overwrite" %in% args
  
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  stage2_assign_plots(
    project_dir = cfg$project_dir %||% ".",
    hex_path = cfg$hex_path %||% "data/hex/hex_grid.geojson",
    hex_layer = cfg$hex_layer %||% NULL,
    overwrite = overwrite
  )
}