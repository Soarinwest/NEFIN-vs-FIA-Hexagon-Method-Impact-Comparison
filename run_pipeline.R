### run_pipeline.R (v2.0) ###
### Foreword ---------------------------------------------------------------------------------
### Title: run_pipeline.R - Master Pipeline Controller
### Author: Soren Donisvitch
### Date: 12/05/2025
### Dependents: R (>= 3.5)
### Description: Master pipeline with FIA processing, NEFIN comparison, and climate analysis
### Foreword: The use or application of these code without permission of the author is prohibited.
###           The author is not liable for the use, modification, or any other application of 
###           this or other provided scripts.

suppressPackageStartupMessages({
  library(fs); library(yaml)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  stages <- list(
    # Core setup stages
    init = "--init" %in% args || "--all" %in% args,
    create_hex = "--create-hex" %in% args,
    hex_filter = "--hex" %in% args || "--all" %in% args,
    fia_pull = "--fia" %in% args || "--all" %in% args,
    
    # Plot processing
    assign_plots = "--assign" %in% args || "--all" %in% args || "--cached" %in% args,
    jitter_library = "--jitter" %in% args || "--all" %in% args || "--cached" %in% args,
    
    # NEFIN stages
    process_nefin = "--nefin" %in% args || "--all" %in% args || "--cached" %in% args,
    assign_nefin = "--nefin" %in% args || "--all" %in% args || "--cached" %in% args,
    
    # Climate stages
    download_prism = "--climate" %in% args || "--all" %in% args,
    extract_prism = "--climate" %in% args || "--all" %in% args || "--cached" %in% args,
    
    # Core metrics
    compute_metrics = "--compute" %in% args || "--all" %in% args || "--cached" %in% args,
    error_analysis = "--analyze" %in% args || "--all" %in% args,
    
    # Comparison & analysis stages
    compare_nefin = "--compare-nefin" %in% args || "--all" %in% args,
    climate_analysis = "--analyze-climate" %in% args,
    
    # Master processing
    process_all_scales = "--all-scales" %in% args,
    master_process = "--master" %in% args,
    advanced_analysis = "--advanced" %in% args,
    spatial_viz = "--spatial-viz" %in% args,
    consolidate_viz = "--viz" %in% args
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
  cat("║    FIA Comprehensive Processing Pipeline v2.0            ║\n")
  cat("║    FIA + NEFIN + Climate Analysis                        ║\n")
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
  
  # =========================================================================
  # STAGE 0: CREATE HEX GRIDS
  # =========================================================================
  if (stages$create_hex) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 0: Hex Grid Creation\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/create_hex_grid.R")
    
    hex_settings <- cfg_proc$hex_creation_settings %||% list()
    area_ha_list <- hex_settings$area_ha_list %||% c(100, 500, 1000, 5000, 10000, 50000, 100000)
    clip_to <- hex_settings$clip_to %||% "data/hex/hex_grid.geojson"
    exclude_water <- hex_settings$exclude_water %||% NULL
    
    create_multiple_hex_grids(
      area_ha_list = area_ha_list,
      clip_to = clip_to,
      exclude_water = exclude_water,
      overwrite = overwrite
    )
    cat("\n")
  }
  
  # =========================================================================
  # STAGE 1: INIT
  # =========================================================================
  if (stages$init) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 1: Project Initialization\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/01_init_project.R")
    paths <- init_project(".")
    cat("✓ Project structure initialized\n\n")
  }
  
  # =========================================================================
  # STAGE 2: HEX FILTER
  # =========================================================================
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
  
  # =========================================================================
  # STAGE 3: FIA PULL
  # =========================================================================
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
  
  # =========================================================================
  # STAGE 4: PLOT ASSIGNMENT
  # =========================================================================
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
  
  # =========================================================================
  # STAGE 5: JITTER LIBRARY
  # =========================================================================
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
  
  # =========================================================================
  # STAGE 6: NEFIN PROCESSING
  # =========================================================================
  if (stages$process_nefin) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 6: NEFIN Data Processing\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/process_nefin_data.R")
    
    nefin_cfg <- cfg_proc$nefin %||% list()
    
    process_nefin_data(
      tree_csv = nefin_cfg$tree_csv %||% "data/raw/nefin/TREE_PLOT_DATA.csv",
      plot_csv = nefin_cfg$plot_csv %||% "data/raw/nefin/NEFIN_plots.csv",
      out_dir = "data/processed",
      overwrite = overwrite
    )
    cat("\n")
  }
  
  # =========================================================================
  # STAGE 7: NEFIN ASSIGNMENT
  # =========================================================================
  if (stages$assign_nefin) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 7: Assign NEFIN to Hex Grids (Multi-Scale)\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/assign_nefin_to_hex.R")
    
    assign_nefin_to_hex(
      nefin_csv = "data/processed/nefin_processed.csv",
      hex_grids = cfg_proc$hex_grids,
      out_file = "data/processed/nefin_hex_assignments.csv",
      overwrite = overwrite
    )
    cat("\n")
  }
  
  # =========================================================================
  # STAGE 8: PRISM DOWNLOAD
  # =========================================================================
  if (stages$download_prism) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 8: Download & Process PRISM Climate Data\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/download_process_prism.R")
    
    prism_cfg <- cfg_proc$prism %||% list()
    
    download_process_prism(
      variables = prism_cfg$variables %||% c("tmean", "ppt"),
      years = prism_cfg$years %||% 2015:2020,
      temporal_period = prism_cfg$temporal_period %||% "annual",
      study_region_path = cfg_proc$hex_path %||% "data/hex/hex_grid.geojson",
      out_dir = "data/processed/prism",
      prism_dir = "data/raw/prism",
      overwrite = overwrite
    )
    cat("\n")
  }
  
  # =========================================================================
  # STAGE 9: EXTRACT PRISM TO HEXES
  # =========================================================================
  if (stages$extract_prism) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 9: Extract PRISM to Hex Grids (Multi-Scale)\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/extract_prism_to_hex.R")
    
    extract_prism_to_hex(
      prism_dir = "data/processed/prism",
      hex_grids = cfg_proc$hex_grids,
      out_file = "data/processed/hex_prism_values.csv",
      overwrite = overwrite
    )
    cat("\n")
  }
  
  # =========================================================================
  # STAGE 10: COMPUTE METRICS
  # =========================================================================
  result_info <- NULL
  if (stages$compute_metrics) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 10: Metric Computation\n")
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
  
  # =========================================================================
  # STAGE 11: ERROR ANALYSIS
  # =========================================================================
  if (stages$error_analysis) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 11: Error Analysis & Visualization\n")
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
  
  # =========================================================================
  # STAGE 12: FIA vs NEFIN COMPARISON
  # =========================================================================
  if (stages$compare_nefin) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 12: FIA vs NEFIN Comparison (All Scales)\n")
    cat("══════════════════════════════════════════════════════════\n")
    
    if (!fs::file_exists("data/processed/nefin_hex_assignments.csv")) {
      cat("⚠ NEFIN assignments not found. Run --nefin first.\n")
    } else {
      safe_source("R/compare_fia_nefin.R")
      
      run_all_scales(
        cfg_path = "configs/process.yml",
        nefin_assign_path = "data/processed/nefin_hex_assignments.csv"
      )
    }
    cat("\n")
  }
  
  # =========================================================================
  # STAGE 13: PROCESS ALL SCALES
  # =========================================================================
  if (stages$process_all_scales) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 13: Process All Grid Scales\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/08_process_all_scales.R")
    
    years_vec <- if (!is.null(cfg_proc$years)) {
      if (is.list(cfg_proc$years) && length(cfg_proc$years) == 2) {
        seq(cfg_proc$years[[1]], cfg_proc$years[[2]])
      } else unlist(cfg_proc$years)
    } else 2018:2020
    
    process_all_scales(
      project_dir = cfg_proc$project_dir %||% ".",
      metric = cfg_proc$metric %||% "aglb",
      years = years_vec,
      level_window = cfg_proc$level_window %||% 3,
      skip_existing = !overwrite
    )
    cat("\n")
  }
  
  # =========================================================================
  # STAGE 14: MASTER PROCESS (FIA + NEFIN, all scales)
  # =========================================================================
  if (stages$master_process) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 14: Master Processing (FIA + NEFIN, All Scales)\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/master_process_all.R")
    
    master_process_all(
      project_dir = cfg_proc$project_dir %||% ".",
      include_nefin = TRUE,
      process_nefin_first = TRUE
    )
    cat("\n")
  }
  
  # =========================================================================
  # STAGE 15: CLIMATE-BIOMASS ANALYSIS
  # =========================================================================
  if (stages$climate_analysis) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 15: Climate-Biomass-Uncertainty Analysis\n")
    cat("══════════════════════════════════════════════════════════\n")
    
    if (!fs::file_exists("data/processed/hex_prism_values.csv")) {
      cat("⚠ PRISM data not found. Run --climate first.\n")
    } else {
      safe_source("R/analyze_climate_biomass.R")
      
      analyze_climate_biomass(
        consolidated_dir = NULL,  # Uses most recent
        output_dir = NULL
      )
    }
    cat("\n")
  }
  
  # =========================================================================
  # STAGE 16: ADVANCED ANALYSIS
  # =========================================================================
  if (stages$advanced_analysis) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 16: Advanced Statistical Analysis\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/advanced_analysis.R")
    
    advanced_analysis(
      consolidated_dir = NULL,
      output_dir = NULL
    )
    cat("\n")
  }
  
  # =========================================================================
  # STAGE 17: SPATIAL VISUALIZATIONS
  # =========================================================================
  if (stages$spatial_viz) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 17: Spatial Visualizations with Hex Grids\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/spatial_visualizations.R")
    
    spatial_visualizations(
      consolidated_dir = NULL,
      output_dir = NULL
    )
    cat("\n")
  }
  
  # =========================================================================
  # STAGE 18: CONSOLIDATED VISUALIZATIONS
  # =========================================================================
  if (stages$consolidate_viz) {
    cat("══════════════════════════════════════════════════════════\n")
    cat("STAGE 18: Consolidated Results Visualization\n")
    cat("══════════════════════════════════════════════════════════\n")
    safe_source("R/visualize_consolidated_results.R")
    
    visualize_results(
      consolidated_dir = NULL,
      output_dir = NULL
    )
    cat("\n")
  }
  
  # =========================================================================
  # PIPELINE COMPLETE
  # =========================================================================
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
    cat("\n╔══════════════════════════════════════════════════════════╗\n")
    cat("║  FIA Comprehensive Processing Pipeline v2.0              ║\n")
    cat("╚══════════════════════════════════════════════════════════╝\n")
    cat("\nUsage: Rscript run_pipeline.R [options]\n\n")
    
    cat("═══════════════════════════════════════════════════════════\n")
    cat("CORE STAGES:\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("  --all            Run all core stages (FIA only)\n")
    cat("  --cached         Run cached stages (4-6) + analysis\n")
    cat("  --init           Initialize project structure\n")
    cat("  --create-hex     Create new hex grid(s)\n")
    cat("  --hex            Filter hex grid\n")
    cat("  --fia            Pull FIA data\n")
    cat("  --assign         Assign plots to hexes\n")
    cat("  --jitter         Build jitter library\n")
    cat("  --compute        Compute metrics (default)\n")
    cat("  --analyze        Error analysis & visualization (default)\n")
    cat("\n")
    
    cat("═══════════════════════════════════════════════════════════\n")
    cat("NEFIN COMPARISON:\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("  --nefin          Process & assign NEFIN data\n")
    cat("  --compare-nefin  Compare FIA vs NEFIN (all scales)\n")
    cat("\n")
    
    cat("═══════════════════════════════════════════════════════════\n")
    cat("CLIMATE ANALYSIS:\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("  --climate        Download & extract PRISM data\n")
    cat("  --analyze-climate Analyze climate-biomass relationships\n")
    cat("\n")
    
    cat("═══════════════════════════════════════════════════════════\n")
    cat("ADVANCED PROCESSING:\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("  --all-scales     Process all grid scales sequentially\n")
    cat("  --master         Master process: FIA+NEFIN all scales\n")
    cat("  --advanced       Advanced statistical analysis\n")
    cat("  --spatial-viz    Create spatial maps with hex grids\n")
    cat("  --viz            Consolidated visualizations\n")
    cat("\n")
    
    cat("═══════════════════════════════════════════════════════════\n")
    cat("OPTIONS:\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("  --overwrite      Force regeneration\n")
    cat("  --grid=NAME      Specify grid for compute stage\n")
    cat("\n")
    
    cat("═══════════════════════════════════════════════════════════\n")
    cat("EXAMPLE WORKFLOWS:\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("\n")
    cat("1. First-time full setup (FIA only):\n")
    cat("   Rscript run_pipeline.R --all\n\n")
    
    cat("2. Full setup with NEFIN and climate:\n")
    cat("   Rscript run_pipeline.R --all --nefin --climate\n\n")
    
    cat("3. Create hex grids from config:\n")
    cat("   Rscript run_pipeline.R --create-hex\n\n")
    
    cat("4. Build cache (if FIA data exists):\n")
    cat("   Rscript run_pipeline.R --cached\n\n")
    
    cat("5. Quick compute + analyze (default):\n")
    cat("   Rscript run_pipeline.R\n\n")
    
    cat("6. Compute specific grid scale:\n")
    cat("   Rscript run_pipeline.R --compute --grid=100ha\n\n")
    
    cat("7. Process all scales sequentially:\n")
    cat("   Rscript run_pipeline.R --all-scales\n\n")
    
    cat("8. Complete research workflow:\n")
    cat("   Rscript run_pipeline.R --master --compare-nefin --analyze-climate --advanced\n\n")
    
    cat("9. NEFIN comparison only:\n")
    cat("   Rscript run_pipeline.R --nefin --compare-nefin\n\n")
    
    cat("10. Climate analysis only:\n")
    cat("    Rscript run_pipeline.R --climate --analyze-climate\n\n")
    
    cat("11. Generate all visualizations:\n")
    cat("    Rscript run_pipeline.R --spatial-viz --viz\n\n")
    
    cat("12. Rebuild jitter library:\n")
    cat("    Rscript run_pipeline.R --jitter --overwrite\n\n")
    
    cat("═══════════════════════════════════════════════════════════\n")
  }
}
