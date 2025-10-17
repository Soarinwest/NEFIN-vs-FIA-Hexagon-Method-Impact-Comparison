# R/06_compute_metrics.R
# Fast metric computation using pre-generated jitter library

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(fs); library(yaml); 
  library(glue); library(tidyr)
})

source("R/utils_spatial.R")
source("R/utils_metrics.R")

`%||%` <- function(a,b) if (!is.null(a)) a else b

stage4_compute_metrics <- function(project_dir = ".",
                                   hex_grid_name = NULL,    # IMPORTANT: which grid to process
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
    stop("Run stage4_build_jitter_library.R first. Missing jitter library in: ", jitter_dir)
  }
  
  jitter_meta <- yaml::read_yaml(jitter_manifest)
  
  # Check how many replicates actually exist
  reps_dir <- fs::path(jitter_dir, "replicates")
  rep_files <- list.files(reps_dir, pattern = "^rep_\\d{4}\\.csv$")
  n_available <- length(rep_files)
  
  if (n_available == 0) {
    stop("No replicate files found in: ", reps_dir)
  }
  
  message("→ Using jitter library:")
  message("    Created: ", jitter_meta$created)
  message("    Expected replicates: ", jitter_meta$n_replicates)
  message("    Available replicates: ", n_available)
  if (n_available < jitter_meta$n_replicates) {
    message("    ⚠ Using partial library (", n_available, " of ", 
            jitter_meta$n_replicates, " replicates)")
  }
  
  metric_col <- switch(tolower(metric),
                       "aglb" = "aglb_Mg_per_ha",
                       "carbon" = "carbon_Mg_per_ha",
                       "mortality" = "mortality_pct_per_yr",
                       "growth" = "growth_Mg_ha_yr",
                       "regeneration" = "regen_ratio",
                       stop("Unknown metric: ", metric)
  )
  
  message("→ Computing metric: ", metric)
  if (!is.null(hex_grid_name)) {
    message("→ Grid scale: ", hex_grid_name)
  }
  
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
  
  # FIXED: Get the correct hex column for this grid scale
  if (!is.null(hex_grid_name)) {
    hex_col <- paste0("hex_id_", hex_grid_name)
    if (!(hex_col %in% names(assignments))) {
      # Fallback: hex_id might be used for fia grid
      if ("hex_id" %in% names(assignments) && hex_grid_name == "fia") {
        hex_col <- "hex_id"
      } else {
        available_hex_cols <- names(assignments)[grepl("^hex_id", names(assignments))]
        stop("Column not found for grid '", hex_grid_name, "': ", hex_col,
             "\n  Available columns: ", paste(available_hex_cols, collapse = ", "))
      }
    }
    message("  Using hex assignment column: ", hex_col)
  } else {
    # Auto-detect (fallback for backward compatibility)
    hex_col_candidates <- c("hex_id", "hex_id_fia", "hex_id_1.5k", "hex_id_3k", "hex_id_6k")
    hex_col <- NULL
    for (candidate in hex_col_candidates) {
      if (candidate %in% names(assignments)) {
        hex_col <- candidate
        break
      }
    }
    if (is.null(hex_col)) {
      stop("No hex_id column found in plot assignments!")
    }
    message("  Auto-detected hex column: ", hex_col)
  }
  
  fia_with_hex <- fia_plot |>
    dplyr::inner_join(
      assignments |> 
        dplyr::select(CN, STATECD, hex_id = !!rlang::sym(hex_col)), 
      by = c("CN", "STATECD")
    )
  
  message("  Plots with metric + hex: ", nrow(fia_with_hex))
  
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  # Include grid name in run_id
  grid_suffix <- if (!is.null(hex_grid_name)) paste0("_", hex_grid_name) else ""
  
  run_id <- run_id %||% cfg$run_id %||%
    paste0(format(Sys.Date(), "%Y-%m-%d"), "_", metric, grid_suffix, "_W", level_window, "y")
  
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
      hex_grid_name = hex_grid_name  # Pass grid name!
    )
    
    out_design[[i]] <- fia_hex_stats |> dplyr::mutate(metric = metric)
    out_pos_sd[[i]] <- pos_sd |> dplyr::mutate(metric = metric)
  }
  
  design <- dplyr::bind_rows(out_design) |> dplyr::mutate(hex_id = as.character(hex_id))
  pos_sd <- dplyr::bind_rows(out_pos_sd) |> dplyr::mutate(hex_id = as.character(hex_id))
  
  joined <- design |>
    dplyr::left_join(pos_sd |> dplyr::select(hex_id, year_label, window, positional_sd, n_reps),
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

add_prism_covariates <- function(results, hex_grid_name, project_dir = ".") {
  prism_file <- fs::path(project_dir, "data", "processed", "hex_prism_values.csv")
  
  if (!fs::file_exists(prism_file)) {
    message("  ℹ No PRISM covariates found. Run extract_prism_to_hex.R to add climate data.")
    return(results)
  }
  
  message("  → Adding PRISM climate covariates...")
  
  # Load PRISM data
  prism_data <- readr::read_csv(prism_file, show_col_types = FALSE)
  
  # Select columns for this grid scale
  grid_cols <- names(prism_data)[grepl(paste0("_", hex_grid_name, "$"), names(prism_data))]
  
  if (length(grid_cols) == 0) {
    message("    ⚠ No PRISM data found for grid: ", hex_grid_name)
    return(results)
  }
  
  # Rename columns to remove grid suffix for joining
  prism_subset <- prism_data |>
    dplyr::select(hex_id, all_of(grid_cols))
  
  names(prism_subset) <- gsub(paste0("_", hex_grid_name, "$"), "", names(prism_subset))
  
  # Join with results
  results_with_prism <- results |>
    dplyr::left_join(prism_subset, by = "hex_id")
  
  n_vars <- length(grid_cols)
  message("    ✓ Added ", n_vars, " PRISM variables")
  
  return(results_with_prism)
}

# Insert this after computing the joined results (around line 160)
joined <- add_prism_covariates(joined, hex_grid_name, project_dir)

# FIXED: Pass hex_grid_name to get correct column
compute_positional_sd_from_library <- function(fia_with_hex, jitter_dir, jitter_meta,
                                               year_label, window_years, metric_col,
                                               hex_grid_name = NULL) {
  
  message("  Computing positional SD from jitter library...")
  
  # Read ALL available replicate CSVs
  reps_dir <- fs::path(jitter_dir, "replicates")
  rep_files <- list.files(reps_dir, pattern = "^rep_\\d{4}\\.csv$", full.names = TRUE)
  
  if (!length(rep_files)) {
    stop("No replicate files found in: ", reps_dir)
  }
  
  message("    Loading ", length(rep_files), " replicate files...")
  jitters <- dplyr::bind_rows(lapply(rep_files, function(f) {
    suppressMessages(readr::read_csv(f, show_col_types = FALSE))
  }))
  
  message("    Loaded ", format(nrow(jitters), big.mark = ","), " jittered coordinates")
  
  # FIXED: Get the correct hex column for this grid scale
  if (!is.null(hex_grid_name)) {
    hex_col <- paste0("hex_id_", hex_grid_name)
    if (!(hex_col %in% names(jitters))) {
      # Fallback for fia grid
      if ("hex_id_jittered" %in% names(jitters) && hex_grid_name == "fia") {
        hex_col <- "hex_id_jittered"
      } else if ("hex_id" %in% names(jitters) && hex_grid_name == "fia") {
        hex_col <- "hex_id"
      } else {
        available_cols <- names(jitters)[grepl("^hex_id", names(jitters))]
        stop("Jitter library missing hex column for grid '", hex_grid_name, "': ", hex_col,
             "\n  Available: ", paste(available_cols, collapse = ", "))
      }
    }
    message("    Using jittered hex column: ", hex_col)
  } else {
    # Auto-detect fallback
    hex_col_candidates <- c("hex_id", "hex_id_jittered", "hex_id_fia")
    hex_col <- NULL
    for (candidate in hex_col_candidates) {
      if (candidate %in% names(jitters)) {
        hex_col <- candidate
        break
      }
    }
    if (is.null(hex_col)) {
      stop("No hex assignment column found in jitter library!")
    }
    message("    Auto-detected hex column: ", hex_col)
  }
  
  # Filter to analysis window
  years <- (year_label - floor(window_years/2)):(year_label + floor(window_years/2))
  fia_year <- fia_with_hex |> dplyr::filter(MEASYEAR %in% years)
  
  if (!nrow(fia_year)) {
    return(dplyr::tibble(hex_id = character(), positional_sd = numeric(), 
                         year_label = integer(), window = character(), n_reps = integer()))
  }
  
  # Get unique replicate IDs from actual data
  rep_ids <- sort(unique(jitters$replicate_id))
  n_reps <- length(rep_ids)
  
  message("    Processing ", n_reps, " replicates...")
  rep_results <- vector("list", n_reps)
  
  for (idx in seq_along(rep_ids)) {
    r <- rep_ids[idx]
    
    # Filter jitter data for this replicate
    jitter_r <- jitters |> dplyr::filter(replicate_id == r)
    
    # Join FIA data with jittered hex assignments
    # KEY: Use the detected hex column name
    data_r <- fia_year |>
      dplyr::inner_join(
        jitter_r |> dplyr::select(CN, STATECD, hex_id_jitter = !!rlang::sym(hex_col)),
        by = c("CN", "STATECD")
      )
    
    # Compute hex means using JITTERED hex assignments
    hex_mean_r <- data_r |>
      dplyr::group_by(hex_id_jitter) |>
      dplyr::summarise(mean_rep = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop") |>
      dplyr::rename(hex_id = hex_id_jitter) |>
      dplyr::mutate(replicate_id = r)
    
    rep_results[[idx]] <- hex_mean_r
  }
  
  all_reps <- dplyr::bind_rows(rep_results)
  
  # Compute positional SD across replicates
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
  
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  get_arg <- function(flag, default = NULL) {
    hit <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
  }
  
  # Get grid name from command line
  grid_name <- get_arg("--grid", NULL)
  
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
    hex_grid_name = grid_name,  # Use command line grid
    hex_path = cfg$hex_path %||% "data/hex/hex_grid.geojson",
    hex_layer = cfg$hex_layer %||% NULL,
    metric = cfg$metric %||% "aglb",
    years = years_vec,
    level_window = cfg$level_window %||% 3,
    run_id = cfg$run_id %||% NULL,
    metric_params = cfg$metric_params %||% list()
  )
}