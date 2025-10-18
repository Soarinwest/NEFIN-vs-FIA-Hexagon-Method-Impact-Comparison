### run_all.R — master runner.R ###
### Foreword ---------------------------------------------------------------------------------
### Title: run_all.R
### Author: Soren Donisvitch
### Date: 10/02/2025
### Dependents: R (>= 3.5)
### Description: run_all.R — master runner: init - hex filter - FIA pull (state-by-state or single) - hex processing
### Foreword: The use or application of these code without permission of the author is prohibited.The author is not ->
###           liable for the use, modification, or any other application of this or other provided scripts.
# run_pipeline.R
# Master pipeline controller - run stages individually or together

suppressPackageStartupMessages({
  library(fs); library(yaml)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  stages <- list(
    init = "--init" %in% args || "--all" %in% args,
    create_hex = "--create-hex" %in% args,
    hex_filter = "--hex" %in% args || "--all" %in% args,
    fia_pull = "--fia" %in% args || "--all" %in% args,
    assign_plots = "--assign" %in% args || "--all" %in% args || "--cached" %in% args,
    jitter_library = "--jitter" %in% args || "--all" %in% args || "--cached" %in% args,
    compute_metrics = "--compute" %in% args || "--all" %in% args || "--cached" %in% args,
    error_analysis = "--analyze" %in% args || "--all" %in% args
  )
  
  overwrite <- "--overwrite" %in% args
  
  # If no stages specified, default to compute + analyze
  if (!any(unlist(stages))) {
    stages$compute_metrics <- TRUE
    stages$error_analysis <- TRUE
  }
  
  list(stages = stages, overwrite = overwrite)
}

safe_source <- function(path) {
  if (!file.exists(path)) stop("Required script missing: ", path, call. = FALSE)
  source(path, local = FALSE)
}

run_pipeline <- function(stages, overwrite = FALSE) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║       FIA Multi-Stage Processing Pipeline                ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Load configurations
  cfg_hex <- if (file.exists("configs/hex_filter.yml")) {
    yaml::read_yaml("configs/hex_filter.yml")
  } else list()
  
  cfg_fia <- if (file.exists("configs/fia_pull.yml")) {
    yaml::read_yaml("configs/fia_pull.yml")
  } else list()
  
  cfg_states <- if (file.exists("configs/fia_states.yml")) {
    yaml::read_yaml("configs/fia_states.yml")
  } else list()
  
  cfg_proc <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else list()
  
  # STAGE 0: CREATE HEX GRIDS
  if (stages$create_hex) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 0: Hex Grid Creation\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/create_hex_grid.R")
    
    # Get hex creation settings from config
    hex_settings <- cfg_proc$hex_creation_settings %||% list()
    area_ha_list <- hex_settings$area_ha_list %||% c(100, 500, 1000, 5000, 10000, 50000, 100000)
    clip_to <- hex_settings$clip_to %||% "data/hex/hex_grid.geojson"
    exclude_water <- hex_settings$exclude_water %||% NULL
    
    # Create all grids
    create_multiple_hex_grids(
      area_ha_list = area_ha_list,
      clip_to = clip_to,
      exclude_water = exclude_water,
      overwrite = overwrite
    )
    cat("\n")
  }
  
  # STAGE 1: INIT
  if (stages$init) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 1: Project Initialization\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/01_init_project.R")
    paths <- init_project(".")
    cat("✓ Project structure initialized\n\n")
  }
  
  # STAGE 2: HEX FILTER
  if (stages$hex_filter) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 2: Hex Grid Filtering\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/02_filter_hex_grid.R")
    
    hex_out <- cfg_hex$out_path %||% "data/hex/hex_grid.geojson"
    
    if (!fs::file_exists(hex_out) || overwrite) {
      if (is.null(cfg_hex$national_hex_path)) {
        stop("configs/hex_filter.yml must set 'national_hex_path'")
      }
      
      filter_hex_grid(
        in_path = cfg_hex$national_hex_path,
        states = cfg_hex$states %||% c("VT","NH","ME","MA","CT","RI","NY"),
        out_path = hex_out,
        method = cfg_hex$method %||% "attribute"
      )
    } else {
      cat("✓ Using existing hex grid: ", hex_out, "\n")
    }
    cat("\n")
  }
  
  # STAGE 3: FIA PULL
  if (stages$fia_pull) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 3: FIA Data Pull\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/03_fia_pull.R")
    
    if (length(cfg_states) && isTRUE(cfg_states$use_states)) {
      if (!file.exists("R/03_fia_pull_states.R")) {
        stop("configs/fia_states.yml requests state-by-state pull, but R/03_fia_pull_states.R is missing.")
      }
      cat("→ Using state-by-state FIA pull\n")
      safe_source("R/03_fia_pull_states.R")
      fia_pull_states(project_dir = ".", cfg_path = "configs/fia_states.yml")
    } else {
      cat("→ Using single-bundle FIA pull\n")
      years_vec <- if (!is.null(cfg_fia$years)) {
        if (is.list(cfg_fia$years) && length(cfg_fia$years) == 2) {
          seq(cfg_fia$years[[1]], cfg_fia$years[[2]])
        } else unlist(cfg_fia$years)
      } else 2008:2022
      
      fia_pull(
        project_dir = ".",
        states2 = cfg_fia$states2 %||% c("VT","NH","ME","MA","CT","RI","NY"),
        years = years_vec,
        eval_type = cfg_fia$eval_type %||% "EXPN",
        db_zip_url = cfg_fia$db_zip_url %||% "https://apps.fs.usda.gov/fia/datamart/Databases/SQLite_FIADB_ENTIRE.zip",
        overwrite_zip = isTRUE(cfg_fia$overwrite_zip %||% FALSE),
        overwrite_unzip = isTRUE(cfg_fia$overwrite_unzip %||% FALSE),
        keep_zip = isTRUE(cfg_fia$keep_zip %||% FALSE)
      )
    }
    cat("\n")
  }
  
  # STAGE 4: PLOT ASSIGNMENT
  if (stages$assign_plots) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 4: Plot-to-Hex Assignment (Multi-Scale)\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/04_assign_plots.R")
    
    stage2_assign_plots(
      project_dir = cfg_proc$project_dir %||% ".",
      hex_grids = cfg_proc$hex_grids,
      overwrite = overwrite
    )
    cat("\n")
  }
  
  # STAGE 5: JITTER LIBRARY
  if (stages$jitter_library) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 5: Monte Carlo Jitter Library (Multi-Scale)\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/05_build_jitter_library.R")
    
    stage4_build_jitter_library(
      project_dir = cfg_proc$project_dir %||% ".",
      hex_grids = cfg_proc$hex_grids,
      n_replicates = cfg_proc$mc_reps %||% 100,
      radius_m = cfg_proc$jitter_radius_m %||% 1609.34,
      use_constraints = TRUE,
      overwrite = overwrite
    )
    cat("\n")
  }
  
  # STAGE 6: COMPUTE METRICS
  result_info <- NULL
  if (stages$compute_metrics) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 6: Metric Computation\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/06_compute_metrics.R")
    
    # Get grid name from command line args
    args <- commandArgs(trailingOnly = TRUE)
    get_arg <- function(flag, default = NULL) {
      hit <- grep(paste0("^", flag, "="), args, value = TRUE)
      if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
    }
    grid_name <- get_arg("--grid", NULL)
    
    # If no specific grid, use first grid in config
    if (is.null(grid_name) && !is.null(cfg_proc$hex_grids)) {
      grid_name <- cfg_proc$hex_grids[[1]]$name
      cat("  → No --grid specified, using: ", grid_name, "\n")
    }
    
    years_vec <- if (!is.null(cfg_proc$years)) {
      if (is.list(cfg_proc$years) && length(cfg_proc$years) == 2) {
        seq(cfg_proc$years[[1]], cfg_proc$years[[2]])
      } else unlist(cfg_proc$years)
    } else 2018:2020
    
    result_info <- stage4_compute_metrics(
      project_dir = cfg_proc$project_dir %||% ".",
      hex_grid_name = grid_name,
      hex_path = cfg_proc$hex_path %||% "data/hex/hex_grid.geojson",
      hex_layer = cfg_proc$hex_layer %||% NULL,
      metric = cfg_proc$metric %||% "aglb",
      years = years_vec,
      level_window = cfg_proc$level_window %||% 3,
      run_id = cfg_proc$run_id %||% NULL
    )
    cat("\n")
  }
  
  # STAGE 7: ERROR ANALYSIS
  if (stages$error_analysis) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 7: Error Analysis & Visualization\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/07_error_analysis.R")
    
    # If we just computed metrics, use that run_dir
    run_dir <- if (!is.null(result_info)) {
      result_info$run_dir
    } else {
      # Find most recent run
      runs_dir <- fs::path(cfg_proc$project_dir %||% ".", "runs")
      if (!fs::dir_exists(runs_dir)) {
        cat("⚠ No runs directory found. Skipping error analysis.\n")
        return(invisible(NULL))
      }
      
      run_dirs <- fs::dir_ls(runs_dir, type = "directory")
      if (!length(run_dirs)) {
        cat("⚠ No run directories found. Skipping error analysis.\n")
        return(invisible(NULL))
      }
      
      # Get most recent
      run_dirs <- run_dirs[order(file.mtime(run_dirs), decreasing = TRUE)]
      run_dirs[1]
    }
    
    cat("→ Analyzing run: ", run_dir, "\n\n")
    
    error_analysis(
      run_dir = run_dir,
      hex_path = cfg_proc$hex_path %||% "data/hex/hex_grid.geojson",
      hex_layer = cfg_proc$hex_layer %||% NULL,
      create_maps = TRUE,
      create_plots = TRUE
    )
    cat("\n")
  }
  
  cat("══════════════════════════════════════════════════════════\n")
  cat("✓ Pipeline Complete\n")
  cat("══════════════════════════════════════════════════════════\n")
}

# Execute
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- parse_args()
  
  if (any(unlist(args$stages))) {
    cat("\nStages to run:\n")
    for (stage in names(args$stages)) {
      if (args$stages[[stage]]) {
        cat("  ✓", stage, "\n")
      }
    }
    cat("\n")
    
    run_pipeline(args$stages, args$overwrite)
  } else {
    cat("\nUsage: Rscript run_pipeline.R [options]\n\n")
    cat("Options:\n")
    cat("  --all            Run all stages (first-time setup)\n")
    cat("  --cached         Run cached stages (4-6) + analysis\n")
    cat("  --init           Initialize project structure\n")
    cat("  --create-hex     Create new hex grid(s)\n")
    cat("  --hex            Filter hex grid\n")
    cat("  --fia            Pull FIA data\n")
    cat("  --assign         Assign plots to hexes\n")
    cat("  --jitter         Build jitter library\n")
    cat("  --compute        Compute metrics (default)\n")
    cat("  --analyze        Error analysis & visualization (default)\n")
    cat("  --overwrite      Force regeneration\n")
    cat("  --grid=NAME      Specify grid for compute stage\n")
    cat("\nExamples:\n")
    cat("  # First-time full setup:\n")
    cat("  Rscript run_pipeline.R --all\n\n")
    cat("  # Create hex grids from config:\n")
    cat("  Rscript run_pipeline.R --create-hex\n\n")
    cat("  # Build cache only (if FIA data already exists):\n")
    cat("  Rscript run_pipeline.R --cached\n\n")
    cat("  # Compute metrics + error analysis (fast, default):\n")
    cat("  Rscript run_pipeline.R\n\n")
    cat("  # Compute for specific grid:\n")
    cat("  Rscript run_pipeline.R --compute --grid=1.5k\n\n")
    cat("  # Rebuild jitter library:\n")
    cat("  Rscript run_pipeline.R --jitter --overwrite\n\n")
    cat("  # Just error analysis on existing run:\n")
    cat("  Rscript run_pipeline.R --analyze\n")
  }
}