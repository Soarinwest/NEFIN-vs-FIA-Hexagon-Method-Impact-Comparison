# R/assign_nefin_to_hex.R
# Assign NEFIN plots to all hex grids (matching FIA assignment structure)

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(fs); library(yaml)
})

source("R/utils_spatial.R")

`%||%` <- function(a,b) if (!is.null(a)) a else b

assign_nefin_to_hex <- function(nefin_csv = "data/processed/nefin_processed.csv",
                                hex_grids = NULL,
                                out_file = "data/processed/nefin_hex_assignments.csv",
                                overwrite = FALSE) {
  
  message("\n═══════════════════════════════════════════════════════════")
  message("NEFIN PLOT-TO-HEX ASSIGNMENT")
  message("═══════════════════════════════════════════════════════════\n")
  
  if (fs::file_exists(out_file) && !overwrite) {
    message("✓ NEFIN hex assignments already exist: ", out_file)
    message("  Use overwrite=TRUE to regenerate")
    return(invisible(out_file))
  }
  
  # Load config to get hex grids
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
    if (!fs::file_exists(grid$path)) {
      stop("Hex grid not found: ", grid$path)
    }
  }
  
  grid_names <- sapply(hex_grids, function(x) x$name)
  message("→ Assigning NEFIN plots to ", length(hex_grids), " hex grids:")
  for (grid in hex_grids) {
    message("    ", grid$name, ": ", grid$path)
  }
  
  # Read NEFIN data
  message("\n→ Reading NEFIN data: ", nefin_csv)
  if (!fs::file_exists(nefin_csv)) {
    stop("NEFIN processed data not found: ", nefin_csv, 
         "\n  Run: Rscript R/process_nefin_data.R")
  }
  
  nefin <- readr::read_csv(nefin_csv, show_col_types = FALSE)
  message("  Plots: ", format(nrow(nefin), big.mark = ","))
  
  # Filter valid coordinates
  nefin_valid <- nefin |>
    dplyr::filter(
      is.finite(lat_public), is.finite(lon_public),
      lat_public >= -90, lat_public <= 90,
      lon_public >= -180, lon_public <= 180
    )
  
  message("  Valid coordinates: ", format(nrow(nefin_valid), big.mark = ","))
  
  if (!nrow(nefin_valid)) stop("No valid NEFIN plots to assign!")
  
  # Turn off spherical geometry
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  crs_ref <- sf::st_crs(5070)
  
  # Load ALL hex grids
  message("\n→ Loading hex grids...")
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
  
  # Convert NEFIN plots to spatial points
  message("\n→ Converting NEFIN plots to spatial points...")
  pts <- sf::st_as_sf(nefin_valid, 
                      coords = c("lon_public", "lat_public"), 
                      crs = 4326, remove = FALSE)
  pts_5070 <- sf::st_transform(pts, crs_ref)
  sf::st_crs(pts_5070) <- crs_ref
  
  # Assign to ALL hex grids
  message("\n→ Assigning NEFIN plots to all hex grids...")
  
  assignments <- sf::st_drop_geometry(pts_5070) |>
    dplyr::select(CN, STATECD, MEASYEAR, lat_public, lon_public, 
                  aglb_Mg_per_ha, source)
  
  for (grid_info in hex_grids_sf) {
    grid_name <- grid_info$name
    grid_sf <- grid_info$sf
    
    message("  → ", grid_name, " (", nrow(grid_sf), " hexes)")
    
    # Spatial join
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
    
    # Report stats
    n_assigned <- sum(!is.na(assignments[[hex_col_name]]))
    n_missing <- sum(is.na(assignments[[hex_col_name]]))
    pct <- round(100 * n_assigned / nrow(assignments), 1)
    
    message("    Assigned: ", format(n_assigned, big.mark = ","), " plots (", pct, "%)")
    if (n_missing > 0) {
      message("    Missing: ", format(n_missing, big.mark = ","), " plots")
    }
  }
  
  # Summary
  message("\n=== Assignment Summary ===")
  message("  Total NEFIN plots: ", nrow(assignments))
  message("  Years: ", paste(range(assignments$MEASYEAR, na.rm = TRUE), collapse = "-"))
  
  for (grid_name in grid_names) {
    hex_col <- paste0("hex_id_", grid_name)
    n_assigned <- sum(!is.na(assignments[[hex_col]]))
    n_unique_hexes <- dplyr::n_distinct(assignments[[hex_col]], na.rm = TRUE)
    message("  ", grid_name, ": ", n_assigned, " plots → ", n_unique_hexes, " hexes")
  }
  
  # Write
  fs::dir_create(fs::path_dir(out_file), recurse = TRUE)
  readr::write_csv(assignments, out_file)
  message("\n✓ Wrote: ", out_file)
  
  message("\n═══════════════════════════════════════════════════════════")
  message("NEFIN ASSIGNMENTS READY")
  message("═══════════════════════════════════════════════════════════\n")
  message("Next step: Compare FIA vs NEFIN by hex")
  message("  Rscript R/compare_fia_nefin.R")
  message("")
  
  invisible(list(
    output = out_file,
    n_plots = nrow(assignments),
    n_grids = length(hex_grids)
  ))
}

# CLI
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  get_arg <- function(flag, default = NULL) {
    hit <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
  }
  
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  assign_nefin_to_hex(
    nefin_csv = get_arg("--input", "data/processed/nefin_processed.csv"),
    hex_grids = cfg$hex_grids,
    out_file = get_arg("--output", "data/processed/nefin_hex_assignments.csv"),
    overwrite = "--overwrite" %in% args
  )
}