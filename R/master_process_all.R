# R/master_process_all.R
# Process FIA and FIA+NEFIN at all scales, consolidate results

suppressPackageStartupMessages({
  library(fs); library(yaml); library(dplyr); library(readr); library(glue)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

master_process_all <- function(project_dir = ".",
                               include_nefin = TRUE,
                               process_nefin_first = TRUE) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Master Processing: FIA + NEFIN, All Grid Scales         ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Load config
  cfg <- if (file.exists("configs/process.yml")) {
    yaml::read_yaml("configs/process.yml")
  } else {
    stop("Missing configs/process.yml")
  }
  
  if (is.null(cfg$hex_grids)) {
    stop("No hex_grids defined in configs/process.yml")
  }
  
  grid_names <- sapply(cfg$hex_grids, function(x) x$name)
  
  cat("Grid scales to process:", length(grid_names), "\n")
  for (gn in grid_names) {
    cat("  -", gn, "\n")
  }
  
  start_time <- Sys.time()
  
  # ========================================================================
  # PHASE 1: Process NEFIN if requested
  # ========================================================================
  
  if (include_nefin && process_nefin_first) {
    cat("\n")
    cat("══════════════════════════════════════════════════════════\n")
    cat("PHASE 1: Processing NEFIN Data\n")
    cat("══════════════════════════════════════════════════════════\n")
    
    # Check if NEFIN files exist
    nefin_tree <- "data/raw/nefin/TREE_PLOT_DATA.csv"
    nefin_plots <- "data/raw/nefin/NEFIN_plots.csv"
    
    if (!fs::file_exists(nefin_tree) || !fs::file_exists(nefin_plots)) {
      cat("⚠ NEFIN files not found, skipping NEFIN processing\n")
      cat("  Missing: ", nefin_tree, " or ", nefin_plots, "\n")
      include_nefin <- FALSE
    } else {
      source("R/process_nefin_data.R")
      source("R/assign_nefin_to_hex.R")
      
      cat("\n→ Step 1: Standardizing NEFIN data...\n")
      process_nefin_data(overwrite = FALSE)
      
      cat("\n→ Step 2: Assigning NEFIN to hex grids...\n")
      assign_nefin_to_hex(overwrite = FALSE)
      
      cat("\n✓ NEFIN preprocessing complete\n")
    }
  }
  
  # ========================================================================
  # PHASE 2: Process FIA at all grid scales
  # ========================================================================
  
  cat("\n")
  cat("══════════════════════════════════════════════════════════\n")
  cat("PHASE 2: Processing FIA at All Grid Scales\n")
  cat("══════════════════════════════════════════════════════════\n")
  
  source("R/06_compute_metrics.R")
  
  fia_results <- list()
  
  for (i in seq_along(grid_names)) {
    grid_name <- grid_names[i]
    
    cat("\n")
    cat("────────────────────────────────────────────────────────\n")
    cat("Processing Grid", i, "of", length(grid_names), ":", grid_name, "\n")
    cat("────────────────────────────────────────────────────────\n")
    
    result <- stage4_compute_metrics(
      project_dir = project_dir,
      hex_grid_name = grid_name,
      hex_path = cfg$hex_path %||% "data/hex/hex_grid.geojson",
      hex_layer = cfg$hex_layer %||% NULL,
      metric = cfg$metric %||% "aglb",
      years = if (!is.null(cfg$years) && is.list(cfg$years) && length(cfg$years) == 2) {
        seq(cfg$years[[1]], cfg$years[[2]])
      } else {
        unlist(cfg$years) %||% 2018:2020
      },
      level_window = cfg$level_window %||% 3,
      run_id = NULL
    )
    
    fia_results[[grid_name]] <- result
  }
  
  # ========================================================================
  # PHASE 3: Consolidate FIA results
  # ========================================================================
  
  cat("\n")
  cat("══════════════════════════════════════════════════════════\n")
  cat("PHASE 3: Consolidating FIA Results\n")
  cat("══════════════════════════════════════════════════════════\n")
  
  fia_combined <- list()
  
  for (grid_name in grid_names) {
    result_file <- fia_results[[grid_name]]$results
    if (fs::file_exists(result_file)) {
      df <- readr::read_csv(result_file, show_col_types = FALSE) |>
        dplyr::mutate(
          grid_scale = grid_name,
          data_source = "FIA"
        )
      fia_combined[[grid_name]] <- df
      cat("  ✓ Loaded", grid_name, ":", nrow(df), "hex-years\n")
    }
  }
  
  fia_all <- dplyr::bind_rows(fia_combined)
  
  # ========================================================================
  # PHASE 4: Process FIA+NEFIN comparisons
  # ========================================================================
  
  nefin_results <- NULL
  
  if (include_nefin) {
    cat("\n")
    cat("══════════════════════════════════════════════════════════\n")
    cat("PHASE 4: Comparing FIA vs NEFIN at All Scales\n")
    cat("══════════════════════════════════════════════════════════\n")
    
    source("R/compare_fia_nefin.R")
    
    comparison_results <- list()
    
    for (grid_name in grid_names) {
      cat("\n→ Comparing at", grid_name, "scale...\n")
      
      comp_result <- tryCatch({
        compare_fia_nefin(
          fia_results = fia_results[[grid_name]]$results,
          nefin_assignments = "data/processed/nefin_hex_assignments.csv",
          hex_grid_name = grid_name,
          years = if (!is.null(cfg$years) && is.list(cfg$years) && length(cfg$years) == 2) {
            seq(cfg$years[[1]], cfg$years[[2]])
          } else {
            unlist(cfg$years) %||% 2018:2020
          },
          level_window = cfg$level_window %||% 3,
          out_dir = fs::path("runs", paste0("fia_nefin_comparison_", grid_name))
        )
        comp_result
      }, error = function(e) {
        cat("  ⚠ Error comparing at", grid_name, ":", e$message, "\n")
        NULL
      })
      
      if (!is.null(comp_result)) {
        comparison_results[[grid_name]] <- comp_result
      }
    }
    
    # Load and combine comparison data
    cat("\n→ Consolidating comparison results...\n")
    
    nefin_combined <- list()
    
    for (grid_name in grid_names) {
      comp_file <- fs::path("runs", paste0("fia_nefin_comparison_", grid_name), 
                            "hex_comparison.csv")
      if (fs::file_exists(comp_file)) {
        df <- readr::read_csv(comp_file, show_col_types = FALSE) |>
          dplyr::mutate(grid_scale = grid_name)
        nefin_combined[[grid_name]] <- df
        cat("  ✓ Loaded", grid_name, "comparison:", nrow(df), "hex-years\n")
      }
    }
    
    if (length(nefin_combined) > 0) {
      nefin_results <- dplyr::bind_rows(nefin_combined)
    }
  }
  
  # ========================================================================
  # PHASE 5: Create consolidated outputs
  # ========================================================================
  
  cat("\n")
  cat("══════════════════════════════════════════════════════════\n")
  cat("PHASE 5: Creating Consolidated Outputs\n")
  cat("══════════════════════════════════════════════════════════\n")
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_dir <- fs::path("runs", paste0("consolidated_", timestamp))
  fs::dir_create(out_dir, recurse = TRUE)
  
  # Write FIA consolidated results
  fia_out <- fs::path(out_dir, "fia_all_scales.csv")
  readr::write_csv(fia_all, fia_out)
  cat("  ✓ Wrote FIA results:", fia_out, "\n")
  cat("    ", nrow(fia_all), "hex-year-scale observations\n")
  
  # Write NEFIN comparison consolidated results
  if (!is.null(nefin_results)) {
    nefin_out <- fs::path(out_dir, "fia_nefin_comparison_all_scales.csv")
    readr::write_csv(nefin_results, nefin_out)
    cat("  ✓ Wrote FIA+NEFIN comparison:", nefin_out, "\n")
    cat("    ", nrow(nefin_results), "hex-year-scale observations\n")
  }
  
  # ========================================================================
  # PHASE 6: Create summary statistics
  # ========================================================================
  
  cat("\n")
  cat("══════════════════════════════════════════════════════════\n")
  cat("PHASE 6: Computing Summary Statistics\n")
  cat("══════════════════════════════════════════════════════════\n")
  
  # FIA summary by grid scale
  fia_summary <- fia_all |>
    dplyr::group_by(grid_scale) |>
    dplyr::summarise(
      n_hex_years = dplyr::n(),
      n_unique_hexes = dplyr::n_distinct(hex_id),
      mean_biomass = mean(mean, na.rm = TRUE),
      mean_sampling_se = mean(se, na.rm = TRUE),
      mean_positional_sd = mean(positional_sd, na.rm = TRUE),
      mean_total_sd = mean(total_sd, na.rm = TRUE),
      mean_n_plots = mean(n_plots, na.rm = TRUE),
      pos_fraction = mean(positional_sd / total_sd, na.rm = TRUE),
      .groups = "drop"
    )
  
  fia_summary_out <- fs::path(out_dir, "fia_summary_by_scale.csv")
  readr::write_csv(fia_summary, fia_summary_out)
  cat("  ✓ Wrote FIA summary:", fia_summary_out, "\n")
  
  # NEFIN comparison summary
  if (!is.null(nefin_results)) {
    nefin_summary <- nefin_results |>
      dplyr::filter(has_both) |>
      dplyr::group_by(grid_scale) |>
      dplyr::summarise(
        n_hex_years_both = dplyr::n(),
        n_unique_hexes_both = dplyr::n_distinct(hex_id),
        mean_fia = mean(mean_fia, na.rm = TRUE),
        mean_nefin = mean(mean_nefin, na.rm = TRUE),
        mean_diff = mean(diff, na.rm = TRUE),
        rmse = sqrt(mean(diff^2, na.rm = TRUE)),
        correlation = cor(mean_fia, mean_nefin, use = "complete.obs"),
        mean_plots_fia = mean(n_plots_fia, na.rm = TRUE),
        mean_plots_nefin = mean(n_plots_nefin, na.rm = TRUE),
        .groups = "drop"
      )
    
    nefin_summary_out <- fs::path(out_dir, "nefin_comparison_summary_by_scale.csv")
    readr::write_csv(nefin_summary, nefin_summary_out)
    cat("  ✓ Wrote NEFIN comparison summary:", nefin_summary_out, "\n")
  }
  
  # ========================================================================
  # PHASE 7: Create text report
  # ========================================================================
  
  cat("\n")
  cat("══════════════════════════════════════════════════════════\n")
  cat("PHASE 7: Creating Summary Report\n")
  cat("══════════════════════════════════════════════════════════\n")
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  
  report_lines <- c(
    "═══════════════════════════════════════════════════════════════",
    "MASTER PROCESSING REPORT: FIA + NEFIN",
    "═══════════════════════════════════════════════════════════════",
    "",
    paste("Generated:", Sys.time()),
    paste("Processing time:", round(elapsed, 1), "minutes"),
    "",
    "───────────────────────────────────────────────────────────────",
    "GRID SCALES PROCESSED",
    "───────────────────────────────────────────────────────────────",
    "",
    paste("Scales:", paste(grid_names, collapse = ", ")),
    "",
    "───────────────────────────────────────────────────────────────",
    "FIA RESULTS BY SCALE",
    "───────────────────────────────────────────────────────────────",
    ""
  )
  
  for (i in seq_len(nrow(fia_summary))) {
    row <- fia_summary[i, ]
    report_lines <- c(report_lines,
                      paste0("Grid: ", row$grid_scale),
                      paste0("  Hex-years: ", row$n_hex_years),
                      paste0("  Unique hexes: ", row$n_unique_hexes),
                      paste0("  Mean biomass: ", round(row$mean_biomass, 2), " Mg/ha"),
                      paste0("  Mean sampling SE: ", round(row$mean_sampling_se, 3), " Mg/ha"),
                      paste0("  Mean positional SD: ", round(row$mean_positional_sd, 3), " Mg/ha"),
                      paste0("  Mean total SD: ", round(row$mean_total_sd, 3), " Mg/ha"),
                      paste0("  Positional fraction: ", round(100 * row$pos_fraction, 1), "%"),
                      paste0("  Mean plots/hex: ", round(row$mean_n_plots, 1)),
                      ""
    )
  }
  
  if (!is.null(nefin_results)) {
    report_lines <- c(report_lines,
                      "───────────────────────────────────────────────────────────────",
                      "FIA vs NEFIN COMPARISON BY SCALE",
                      "───────────────────────────────────────────────────────────────",
                      ""
    )
    
    for (i in seq_len(nrow(nefin_summary))) {
      row <- nefin_summary[i, ]
      report_lines <- c(report_lines,
                        paste0("Grid: ", row$grid_scale),
                        paste0("  Hexes with both: ", row$n_unique_hexes_both),
                        paste0("  Mean FIA: ", round(row$mean_fia, 2), " Mg/ha"),
                        paste0("  Mean NEFIN: ", round(row$mean_nefin, 2), " Mg/ha"),
                        paste0("  Mean difference: ", round(row$mean_diff, 2), " Mg/ha"),
                        paste0("  RMSE: ", round(row$rmse, 2), " Mg/ha"),
                        paste0("  Correlation: ", round(row$correlation, 3)),
                        paste0("  Mean FIA plots: ", round(row$mean_plots_fia, 1)),
                        paste0("  Mean NEFIN plots: ", round(row$mean_plots_nefin, 1)),
                        ""
      )
    }
  }
  
  report_lines <- c(report_lines,
                    "───────────────────────────────────────────────────────────────",
                    "OUTPUT FILES",
                    "───────────────────────────────────────────────────────────────",
                    "",
                    "Consolidated data:",
                    paste0("  ", fia_out),
                    if (!is.null(nefin_results)) paste0("  ", nefin_out) else NULL,
                    "",
                    "Summary statistics:",
                    paste0("  ", fia_summary_out),
                    if (!is.null(nefin_results)) paste0("  ", nefin_summary_out) else NULL,
                    "",
                    "Individual run directories:",
                    paste0("  ", fs::path("runs", "2025-*_aglb_*_W3y")),
                    if (include_nefin) paste0("  ", fs::path("runs", "fia_nefin_comparison_*")) else NULL,
                    "",
                    "═══════════════════════════════════════════════════════════════"
  )
  
  report_file <- fs::path(out_dir, "processing_report.txt")
  writeLines(report_lines, report_file)
  cat("  ✓ Wrote report:", report_file, "\n")
  
  # ========================================================================
  # DONE
  # ========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  PROCESSING COMPLETE                                      ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  cat("Processing time:", round(elapsed, 1), "minutes\n")
  cat("\n")
  cat("Consolidated outputs:\n")
  cat("  ", out_dir, "\n")
  cat("    - fia_all_scales.csv\n")
  if (!is.null(nefin_results)) {
    cat("    - fia_nefin_comparison_all_scales.csv\n")
  }
  cat("    - fia_summary_by_scale.csv\n")
  if (!is.null(nefin_results)) {
    cat("    - nefin_comparison_summary_by_scale.csv\n")
  }
  cat("    - processing_report.txt\n")
  cat("\n")
  cat("Next steps:\n")
  cat("  1. Load consolidated CSVs for visualization\n")
  cat("  2. Compare error across scales\n")
  cat("  3. Analyze FIA vs FIA+NEFIN\n")
  cat("\n")
  
  invisible(list(
    fia_all = fia_all,
    nefin_results = nefin_results,
    fia_summary = fia_summary,
    nefin_summary = if (!is.null(nefin_results)) nefin_summary else NULL,
    out_dir = out_dir
  ))
}

# CLI
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  include_nefin <- !("--fia-only" %in% args)
  
  master_process_all(
    project_dir = ".",
    include_nefin = include_nefin,
    process_nefin_first = TRUE
  )
}