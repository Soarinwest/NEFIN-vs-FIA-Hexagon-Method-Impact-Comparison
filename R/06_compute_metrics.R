#R/stage4_compute_metrics.R
# Fast metric computation using pre-generated jitter library

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(fs); library(yaml); 
  library(glue); library(tidyr)
})

source("R/utils_spatial.R")
source("R/utils_metrics.R")

`%||%` <- function(a,b) if (!is.null(a)) a else b

stage4_compute_metrics <- function(project_dir = ".",
                                   hex_path = "data/hex/hex_grid.geojson",
                                   hex_layer = NULL,
                                   metric = "aglb",
                                   years = 2018:2020,
                                   level_window = 3,
                                   run_id = NULL,
                                   metric_params = list()) {
  
  assignments_file <- fs::path(project_dir, "data", "processed", "plot_hex_assignments.csv")
  jitter_dir <- fs::path(project_dir, "data", "processed", "mc_jitter_library")
  jitter_manifest <- fs::path(jitter_dir, "manifest.yml")
  
  if (!fs::file_exists(assignments_file)) {
    stop("Run stage2_assign_plots.R first. Missing: ", assignments_file)
  }
  if (!fs::dir_exists(jitter_dir) || !fs::file_exists(jitter_manifest)) {
    stop("Run stage3_build_jitter_library.R first. Missing jitter library in: ", jitter_dir)
  }
  
  jitter_meta <- yaml::read_yaml(jitter_manifest)
  message("→ Using jitter library:")
  message("    Created: ", jitter_meta$created)
  message("    Replicates: ", jitter_meta$n_replicates)
  message("    Format: ", jitter_meta$format)
  
  metric_col <- switch(tolower(metric),
                       "aglb" = "aglb_Mg_per_ha",
                       "carbon" = "carbon_Mg_per_ha",
                       "mortality" = "mortality_pct_per_yr",
                       "growth" = "growth_Mg_ha_yr",
                       "regeneration" = "regen_ratio",
                       stop("Unknown metric: ", metric)
  )
  
  message("→ Computing metric: ", metric)
  fia_root <- if (fs::dir_exists(fs::path(project_dir, "data", "interim", "fia_region"))) {
    fs::path(project_dir, "data", "interim", "fia_region")
  } else {
    fs::path(project_dir, "data", "interim", "fia")
  }
  
  # Build plot-level metric
  fia_plot <- build_plot_aglb(
    tree_csv = fs::path(fia_root, "tree.csv"),
    plot_csv = fs::path(fia_root, "plot.csv")
  )
  
  message("→ Joining with hex assignments...")
  assignments <- readr::read_csv(assignments_file, show_col_types = FALSE)
  
  fia_with_hex <- fia_plot |>
    dplyr::inner_join(assignments |> dplyr::select(CN, STATECD, hex_id), 
                      by = c("CN", "STATECD"))
  
  message("  Plots with metric + hex: ", nrow(fia_with_hex))
  
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  run_id <- run_id %||% cfg$run_id %||%
    paste0(format(Sys.Date(), "%Y-%m-%d"), "_", metric, "_W", level_window, "y")
  out_dir <- fs::path(project_dir, "runs", run_id)
  fs::dir_create(out_dir, recurse = TRUE)
  message("→ Output dir: ", out_dir)
  
  out_design <- vector("list", length(years))
  out_pos_sd <- vector("list", length(years))
  
  for (i in seq_along(years)) {
    yr <- years[i]
    message(glue("\n=== Year {yr} ==="))
    
    fia_hex_stats <- aggregate_hex(fia_with_hex, year_label = yr, 
                                   window_years = level_window, y_col = metric_col) |>
      dplyr::mutate(source = "fia")
    
    pos_sd <- compute_positional_sd_from_library(
      fia_with_hex = fia_with_hex,
      jitter_dir = jitter_dir,
      jitter_meta = jitter_meta,
      year_label = yr,
      window_years = level_window,
      metric_col = metric_col,
      hex_path = hex_path,
      hex_layer = hex_layer
    )
    
    out_design[[i]] <- fia_hex_stats |> dplyr::mutate(metric = metric)
    out_pos_sd[[i]] <- pos_sd |> dplyr::mutate(metric = metric)
  }
  
  design <- dplyr::bind_rows(out_design) |> dplyr::mutate(hex_id = as.character(hex_id))
  pos_sd <- dplyr::bind_rows(out_pos_sd) |> dplyr::mutate(hex_id = as.character(hex_id))
  
  joined <- design |>
    dplyr::left_join(pos_sd |> dplyr::select(hex_id, year_label, window, positional_sd),
                     by = c("hex_id", "year_label", "window")) |>
    dplyr::mutate(
      total_sd = sqrt(se^2 + dplyr::coalesce(positional_sd, 0)^2)
    )
  
  out_file <- fs::path(out_dir, glue("hex_{metric}_results.csv"))
  readr::write_csv(joined, out_file)
  message("\n✓ Wrote: ", out_file)
  
  invisible(list(
    results = out_file,
    run_dir = out_dir,
    run_id = run_id
  ))
}

compute_positional_sd_from_library <- function(fia_with_hex, jitter_dir, jitter_meta,
                                               year_label, window_years, metric_col,
                                               hex_path, hex_layer = NULL) {
  
  message("  Computing positional SD from jitter library...")
  
  jitter_file <- fs::path(jitter_dir, jitter_meta$file)
  if (jitter_meta$format == "parquet") {
    jitters <- arrow::read_parquet(jitter_file)
  } else {
    jitters <- readr::read_csv(jitter_file, show_col_types = FALSE)
  }
  
  message("    Loaded ", nrow(jitters), " jittered coordinates")
  
  years <- (year_label - floor(window_years/2)):(year_label + floor(window_years/2))
  fia_year <- fia_with_hex |> dplyr::filter(MEASYEAR %in% years)
  
  if (!nrow(fia_year)) {
    return(dplyr::tibble(hex_id = character(), positional_sd = numeric(), 
                         year_label = integer(), window = character()))
  }
  
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  hx <- if (is.null(hex_layer)) sf::st_read(hex_path, quiet = TRUE) else sf::st_read(hex_path, layer = hex_layer, quiet = TRUE)
  if (!("hex_id" %in% names(hx))) {
    if ("ID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = ID)
    else if ("OBJECTID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = OBJECTID)
    else hx$hex_id <- seq_len(nrow(hx))
  }
  hx$hex_id <- as.character(hx$hex_id)
  hx_5070 <- sf::st_transform(hx, 5070)
  
  n_reps <- jitter_meta$n_replicates
  rep_results <- vector("list", n_reps)
  
  message("    Processing ", n_reps, " replicates...")
  for (r in seq_len(n_reps)) {
    jitter_r <- jitters |> dplyr::filter(replicate_id == r)
    
    data_r <- fia_year |>
      dplyr::inner_join(jitter_r |> dplyr::select(CN, STATECD, lon_jittered, lat_jittered),
                        by = c("CN", "STATECD"))
    
    pts_r <- sf::st_as_sf(data_r, coords = c("lon_jittered", "lat_jittered"), 
                          crs = 4326, remove = FALSE)
    pts_r_5070 <- sf::st_transform(pts_r, 5070)
    
    joined_r <- sf::st_join(pts_r_5070, hx_5070["hex_id"], left = TRUE, join = sf::st_intersects) |>
      sf::st_drop_geometry()
    
    hex_mean_r <- joined_r |>
      dplyr::group_by(hex_id) |>
      dplyr::summarise(mean_rep = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(replicate_id = r)
    
    rep_results[[r]] <- hex_mean_r
  }
  
  all_reps <- dplyr::bind_rows(rep_results)
  
  pos_sd <- all_reps |>
    dplyr::group_by(hex_id) |>
    dplyr::summarise(
      positional_sd = sd(mean_rep, na.rm = TRUE),
      n_reps = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      year_label = year_label,
      window = paste0(window_years, "y")
    )
  
  message("    Computed positional SD for ", nrow(pos_sd), " hexes")
  
  pos_sd
}

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  years_vec <- if (!is.null(cfg$years)) {
    if (is.list(cfg$years) && length(cfg$years) == 2) {
      seq(cfg$years[[1]], cfg$years[[2]])
    } else {
      unlist(cfg$years)
    }
  } else {
    2018:2020
  }
  
  stage4_compute_metrics(
    project_dir = cfg$project_dir %||% ".",
    hex_path = cfg$hex_path %||% "data/hex/hex_grid.geojson",
    hex_layer = cfg$hex_layer %||% NULL,
    metric = cfg$metric %||% "aglb",
    years = years_vec,
    level_window = cfg$level_window %||% 3,
    run_id = cfg$run_id %||% NULL,
    metric_params = cfg$metric_params %||% list()
  )
}