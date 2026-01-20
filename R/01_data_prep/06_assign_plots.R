# R/04_assign_plots.R (MULTI-SCALE)
# Assign plots to ALL hex grids simultaneously

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(fs); library(yaml)
})

source("R/utils_spatial.R")
source("R/utils_metrics.R")

`%||%` <- function(a,b) if (!is.null(a)) a else b

#' Stage 2: Assign plots to ALL hex grids and cache
#' This runs ONCE and stores plot→hex mappings for ALL grids with original coords
stage2_assign_plots <- function(project_dir = ".",
                                hex_grids = NULL,
                                overwrite = FALSE) {
  
  # Load config to get hex grids if not provided
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  if (is.null(hex_grids)) {
    hex_grids <- cfg$hex_grids
    if (is.null(hex_grids)) {
      stop("No hex_grids specified. Set in function call or configs/process.yml")
    }
  }
  
  # Validate hex grids
  for (i in seq_along(hex_grids)) {
    grid <- hex_grids[[i]]
    if (is.null(grid$name)) stop("hex_grids[[", i, "]] missing 'name'")
    if (is.null(grid$path)) stop("hex_grids[[", i, "]] missing 'path'")
    if (!fs::file_exists(grid$path)) {
      stop("Hex grid not found: ", grid$path)
    }
  }
  
  grid_names <- sapply(hex_grids, function(x) x$name)
  message("→ Multi-scale plot assignment with ", length(hex_grids), " hex grids:")
  for (grid in hex_grids) {
    message("    ", grid$name, ": ", grid$path)
  }
  
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
    # Check if it has all grid assignments
    existing <- suppressMessages(readr::read_csv(out_file, n_max = 5, show_col_types = FALSE))
    expected_cols <- paste0("hex_id_", grid_names)
    has_cols <- all(expected_cols %in% names(existing))
    
    if (has_cols) {
      message("✓ Multi-scale plot assignments already exist: ", out_file)
      message("  Use overwrite=TRUE to regenerate")
      return(invisible(out_file))
    } else {
      message("⚠ Existing assignments missing some grids. Regenerating...")
    }
  }
  
  fs::dir_create(out_dir, recurse = TRUE)
  
  message("\n→ Reading FIA PLOT data from: ", fia_root)
  plot_csv <- fs::path(fia_root, "plot.csv")
  
  if (!fs::file_exists(plot_csv)) stop("PLOT not found: ", plot_csv)
  
  pl <- suppressMessages(readr::read_csv(plot_csv, guess_max = 1e6, show_col_types = FALSE))
  
  message("  PLOT rows: ", nrow(pl))
  
  # Normalize coordinates
  pl <- normalize_plot_coords(pl)
  
  # Store ORIGINAL coordinates explicitly
  pl <- pl |>
    dplyr::mutate(
      lat_original = lat_public,
      lon_original = lon_public
    )
  
  # Filter valid coordinates
  valid_plots <- pl |>
    dplyr::filter(
      is.finite(lat_original), is.finite(lon_original),
      lat_original >= -90, lat_original <= 90,
      lon_original >= -180, lon_original <= 180
    )
  
  message("  Valid plots with coordinates: ", nrow(valid_plots))
  
  if (!nrow(valid_plots)) stop("No valid plots to assign!")
  
  # Turn off spherical geometry
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  crs_ref <- sf::st_crs(5070)
  
  # Load ALL hex grids for assignment
  message("\n→ Loading all hex grids for multi-scale assignment...")
  hex_grids_sf <- lapply(hex_grids, function(grid) {
    message("    Loading ", grid$name, "...")
    hx <- if (is.null(grid$layer)) {
      sf::st_read(grid$path, quiet = TRUE)
    } else {
      sf::st_read(grid$path, layer = grid$layer, quiet = TRUE)
    }
    
    if (!("hex_id" %in% names(hx))) {
      if ("ID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = ID)
      else if ("OBJECTID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = OBJECTID)
      else hx$hex_id <- seq_len(nrow(hx))
    }
    hx$hex_id <- as.character(hx$hex_id)
    hx <- sf::st_make_valid(sf::st_zm(hx, drop = TRUE, what = "ZM"))
    hx_5070 <- sf::st_transform(hx, crs_ref)
    sf::st_crs(hx_5070) <- crs_ref
    
    list(name = grid$name, sf = hx_5070)
  })
  
  # Convert plots to spatial points
  message("\n→ Converting plots to spatial points...")
  pts <- sf::st_as_sf(valid_plots, 
                      coords = c("lon_original", "lat_original"), 
                      crs = 4326, remove = FALSE)
  pts_5070 <- sf::st_transform(pts, crs_ref)
  sf::st_crs(pts_5070) <- crs_ref
  
  # Assign plots to ALL hex grids
  message("\n→ Assigning plots to all ", length(hex_grids_sf), " hex grids...")
  
  assignments <- sf::st_drop_geometry(pts_5070) |>
    dplyr::select(CN, STATECD, MEASYEAR, UNITCD, COUNTYCD, PLOT, 
                  lat_original, lon_original)
  
  for (grid_info in hex_grids_sf) {
    grid_name <- grid_info$name
    grid_sf <- grid_info$sf
    
    message("  → ", grid_name, " (", nrow(grid_sf), " hexes)")
    
    # Spatial join for this grid
    joined <- sf::st_join(pts_5070, grid_sf["hex_id"], left = TRUE, join = sf::st_intersects)
    
    # Extract hex_id for this grid
    hex_col_name <- paste0("hex_id_", grid_name)
    
    if ("hex_id.y" %in% names(joined)) {
      assignments[[hex_col_name]] <- joined$hex_id.y
    } else if ("hex_id.x" %in% names(joined)) {
      assignments[[hex_col_name]] <- joined$hex_id.x
    } else {
      assignments[[hex_col_name]] <- joined$hex_id
    }
    
    # Report assignment stats
    n_assigned <- sum(!is.na(assignments[[hex_col_name]]))
    n_missing <- sum(is.na(assignments[[hex_col_name]]))
    pct <- round(100 * n_assigned / nrow(assignments), 1)
    
    message("    Assigned: ", format(n_assigned, big.mark = ","), " plots (", pct, "%)")
    if (n_missing > 0) {
      message("    Missing: ", format(n_missing, big.mark = ","), " plots")
    }
  }
  
  # Summary statistics
  message("\n=== Assignment Summary ===")
  message("  Total plots: ", nrow(assignments))
  message("  Years: ", paste(range(assignments$MEASYEAR, na.rm = TRUE), collapse = "-"))
  
  for (grid_name in grid_names) {
    hex_col <- paste0("hex_id_", grid_name)
    n_assigned <- sum(!is.na(assignments[[hex_col]]))
    n_unique_hexes <- dplyr::n_distinct(assignments[[hex_col]], na.rm = TRUE)
    message("  ", grid_name, ": ", n_assigned, " plots → ", n_unique_hexes, " hexes")
  }
  
  # Write
  readr::write_csv(assignments, out_file)
  message("\n✓ Wrote: ", out_file)
  
  # Summary file
  summary_file <- fs::path(out_dir, "plot_hex_summary.txt")
  summary_lines <- c(
    "Multi-Scale Plot-to-Hex Assignment Summary",
    paste("Generated:", Sys.time()),
    "",
    paste("Total plots:", nrow(assignments)),
    paste("Year range:", paste(range(assignments$MEASYEAR, na.rm = TRUE), collapse = "-")),
    "",
    "Hex Grids:",
    sapply(seq_along(hex_grids), function(i) {
      grid <- hex_grids[[i]]
      hex_col <- paste0("hex_id_", grid$name)
      n_assigned <- sum(!is.na(assignments[[hex_col]]))
      n_hexes <- dplyr::n_distinct(assignments[[hex_col]], na.rm = TRUE)
      paste0("  ", grid$name, ": ", grid$path)
    }),
    "",
    "Assignments per grid:",
    sapply(grid_names, function(gn) {
      hex_col <- paste0("hex_id_", gn)
      n_assigned <- sum(!is.na(assignments[[hex_col]]))
      n_hexes <- dplyr::n_distinct(assignments[[hex_col]], na.rm = TRUE)
      pct <- round(100 * n_assigned / nrow(assignments), 1)
      paste0("  ", gn, ": ", n_assigned, " plots (", pct, "%) → ", n_hexes, " hexes")
    }),
    "",
    "Files:",
    paste("  Assignments:", out_file)
  )
  
  writeLines(summary_lines, summary_file)
  
  invisible(list(
    assignments = out_file,
    summary = summary_file,
    n_plots = nrow(assignments),
    n_grids = length(hex_grids)
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
    hex_grids = cfg$hex_grids,
    overwrite = overwrite
  )
}