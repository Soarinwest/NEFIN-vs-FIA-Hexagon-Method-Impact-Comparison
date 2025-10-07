# R/stage3_build_jitter_library.R
# Pre-generate N jittered coordinate sets for reuse

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(fs); library(yaml); library(glue)
})

source("R/utils_spatial.R")
source("R/utils_metrics.R")

`%||%` <- function(a,b) if (!is.null(a)) a else b

stage3_build_jitter_library <- function(project_dir = ".",
                                        hex_path = "data/hex/hex_grid.geojson",
                                        hex_layer = NULL,
                                        n_replicates = 100,
                                        radius_m = 1609.34,
                                        use_constraints = TRUE,
                                        overwrite = FALSE,
                                        format = "parquet") {
  
  assignments_file <- fs::path(project_dir, "data", "processed", "plot_hex_assignments.csv")
  if (!fs::file_exists(assignments_file)) {
    stop("Run stage2_assign_plots.R first. Missing: ", assignments_file)
  }
  
  out_dir <- fs::path(project_dir, "data", "processed", "mc_jitter_library")
  manifest_file <- fs::path(out_dir, "manifest.yml")
  
  if (fs::dir_exists(out_dir) && fs::file_exists(manifest_file) && !overwrite) {
    manifest <- yaml::read_yaml(manifest_file)
    message("✓ Jitter library already exists:")
    message("    Created: ", manifest$created)
    message("    Replicates: ", manifest$n_replicates)
    message("    Plots: ", manifest$n_plots)
    message("  Use overwrite=TRUE to regenerate")
    return(invisible(out_dir))
  }
  
  fs::dir_create(out_dir, recurse = TRUE)
  
  message("→ Reading plot assignments...")
  assignments <- readr::read_csv(assignments_file, show_col_types = FALSE)
  
  valid_plots <- assignments |>
    dplyr::filter(!is.na(hex_id), 
                  is.finite(lat_original), is.finite(lon_original))
  
  message("  Total plots with hex assignment: ", nrow(valid_plots))
  
  if (!nrow(valid_plots)) stop("No valid plots to jitter!")
  
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  use_hex_union <- use_constraints && isTRUE(cfg$mask$use_hex_union)
  use_state <- use_constraints && isTRUE(cfg$mask$use_state_constraint)
  
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  crs_ref <- sf::st_crs(5070)
  
  hex_union_5070 <- if (use_hex_union) {
    message("  Building hex union for jitter constraints...")
    hu <- build_hex_union_5070(hex_path, hex_layer)
    sf::st_crs(hu) <- crs_ref
    hu
  } else NULL
  
  state_polys_5070 <- if (use_state) {
    message("  Loading state boundaries for jitter constraints...")
    sp <- build_state_polys_5070(cfg$mask$state_geo_path %||% "", 
                                 cfg$mask$state_field %||% "STATEFP")
    sf::st_crs(sp) <- crs_ref
    sp
  } else NULL
  
  message("→ Converting to spatial points...")
  pts <- sf::st_as_sf(valid_plots, 
                      coords = c("lon_original", "lat_original"), 
                      crs = 4326, remove = FALSE)
  pts_5070 <- sf::st_transform(pts, crs_ref)
  sf::st_crs(pts_5070) <- crs_ref
  
  message("→ Generating ", n_replicates, " jittered coordinate sets...")
  message("  This may take a while...")
  
  jitter_results <- vector("list", n_replicates)
  
  report_every <- max(1, floor(n_replicates / 10))
  start_time <- Sys.time()
  
  for (r in seq_len(n_replicates)) {
    if (r %% report_every == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      rate <- r / elapsed
      remaining <- (n_replicates - r) / rate
      message(sprintf("  Progress: %d/%d (%.1f%%) - ETA: %.1f min", 
                      r, n_replicates, 100 * r / n_replicates, remaining / 60))
    }
    
    if (!is.null(hex_union_5070) || !is.null(state_polys_5070)) {
      pts_j <- constrained_jitter_once(pts_5070, radius_m, hex_union_5070, 
                                       state_polys_5070, 
                                       max_reroll = as.integer(cfg$mask$max_reroll %||% 20))
    } else {
      u <- runif(nrow(pts_5070))
      rr <- sqrt(u) * radius_m
      theta <- runif(nrow(pts_5070), 0, 2 * pi)
      dx <- rr * cos(theta)
      dy <- rr * sin(theta)
      pts_j <- sf::st_set_geometry(pts_5070, sf::st_geometry(pts_5070) + cbind(dx, dy))
      sf::st_crs(pts_j) <- crs_ref
    }
    
    pts_j_4326 <- sf::st_transform(pts_j, 4326)
    coords_j <- sf::st_coordinates(pts_j_4326)
    
    jitter_results[[r]] <- valid_plots |>
      dplyr::select(CN, STATECD, hex_id) |>
      dplyr::mutate(
        replicate_id = r,
        lon_jittered = coords_j[, 1],
        lat_jittered = coords_j[, 2]
      )
  }
  
  message("✓ Jittering complete in ", 
          round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1), " minutes")
  
  message("→ Writing jitter library...")
  all_jitters <- dplyr::bind_rows(jitter_results)
  
  if (tolower(format) == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      message("  'arrow' package not available, falling back to CSV")
      format <- "csv"
    }
  }
  
  if (tolower(format) == "parquet") {
    out_file <- fs::path(out_dir, "jitter_library.parquet")
    arrow::write_parquet(all_jitters, out_file, compression = "zstd")
  } else {
    out_file <- fs::path(out_dir, "jitter_library.csv")
    readr::write_csv(all_jitters, out_file)
  }
  
  message("✓ Wrote: ", out_file)
  message("  Size: ", format(file.size(out_file) / 1024^2, digits = 1), " MB")
  
  manifest <- list(
    created = as.character(Sys.time()),
    n_replicates = n_replicates,
    n_plots = nrow(valid_plots),
    total_jitters = nrow(all_jitters),
    radius_m = radius_m,
    format = format,
    constrained = use_constraints,
    used_hex_union = !is.null(hex_union_5070),
    used_state_constraint = !is.null(state_polys_5070),
    file = basename(out_file)
  )
  
  yaml::write_yaml(manifest, manifest_file)
  message("✓ Wrote: ", manifest_file)
  
  message("\n=== Jitter Library Summary ===")
  message("  Replicates:     ", n_replicates)
  message("  Plots:          ", nrow(valid_plots))
  message("  Total records:  ", nrow(all_jitters))
  message("  Format:         ", format)
  message("  Constrained:    ", use_constraints)
  
  invisible(list(
    library_dir = out_dir,
    library_file = out_file,
    manifest = manifest_file
  ))
}

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  overwrite <- "--overwrite" %in% args
  
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  stage3_build_jitter_library(
    project_dir = cfg$project_dir %||% ".",
    hex_path = cfg$hex_path %||% "data/hex/hex_grid.geojson",
    hex_layer = cfg$hex_layer %||% NULL,
    n_replicates = cfg$mc_reps %||% 100,
    radius_m = cfg$jitter_radius_m %||% 1609.34,
    use_constraints = TRUE,
    overwrite = overwrite,
    format = "parquet"
  )
}