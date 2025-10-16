# R/08_process_all_scales.R
# Process metrics for all hex grid scales and create comparison report

suppressPackageStartupMessages({
  library(fs); library(yaml); library(dplyr); library(readr); library(ggplot2)
})

source("R/06_compute_metrics.R")
source("R/07_error_analysis.R")

`%||%` <- function(a,b) if (!is.null(a)) a else b

process_all_scales <- function(project_dir = ".",
                               metric = "aglb",
                               years = 2018:2020,
                               level_window = 3,
                               skip_existing = TRUE) {
  
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
  
  message("\n")
  message("╔══════════════════════════════════════════════════════════╗")
  message("║       Multi-Scale Processing Pipeline                    ║")
  message("╚══════════════════════════════════════════════════════════╝")
  message("\nProcessing ", length(grid_names), " hex grid scales:")
  for (grid in cfg$hex_grids) {
    message("  - ", grid$name, ": ", grid$path)
  }
  message("\n")
  
  # Process each grid
  results <- list()
  
  for (i in seq_along(grid_names)) {
    grid_name <- grid_names[i]
    grid_info <- cfg$hex_grids[[i]]
    
    message("\n")
    message("══════════════════════════════════════════════════════════")
    message("Processing Grid ", i, "/", length(grid_names), ": ", grid_name)
    message("══════════════════════════════════════════════════════════")
    
    # Check if already processed
    run_id <- paste0(format(Sys.Date(), "%Y-%m-%d"), "_", metric, "_", grid_name, "_W", level_window, "y")
    run_dir <- fs::path(project_dir, "runs", run_id)
    result_file <- fs::path(run_dir, paste0("hex_", metric, "_results.csv"))
    
    if (skip_existing && fs::file_exists(result_file)) {
      message("✓ Already processed: ", result_file)
      message("  Skipping computation (use skip_existing=FALSE to reprocess)")
      
      results[[grid_name]] <- list(
        run_dir = run_dir,
        results_file = result_file,
        grid_name = grid_name
      )
      next
    }
    
    # FIXED: Pass the grid name to compute metrics!
    message("\n→ Computing metrics...")
    result_info <- stage4_compute_metrics(
      project_dir = project_dir,
      hex_grid_name = grid_name,  # CRITICAL: Pass the grid name!
      hex_path = grid_info$path,   # Use grid-specific path
      hex_layer = grid_info$layer,
      metric = metric,
      years = years,
      level_window = level_window,
      run_id = run_id
    )
    
    # Run error analysis with correct hex grid
    message("\n→ Running error analysis...")
    error_analysis(
      run_dir = result_info$run_dir,
      hex_path = grid_info$path,  # Use grid-specific path!
      hex_layer = grid_info$layer,
      create_maps = TRUE,
      create_plots = TRUE
    )
    
    results[[grid_name]] <- list(
      run_dir = result_info$run_dir,
      results_file = result_info$results,
      grid_name = grid_name
    )
  }
  
  # Create comparison report
  message("\n")
  message("══════════════════════════════════════════════════════════")
  message("Creating Multi-Scale Comparison Report")
  message("══════════════════════════════════════════════════════════")
  
  comparison_dir <- fs::path(project_dir, "runs", 
                             paste0(format(Sys.Date(), "%Y-%m-%d"), "_multiscale_comparison"))
  fs::dir_create(comparison_dir, recurse = TRUE)
  
  create_comparison_report(results, comparison_dir, metric, years)
  
  message("\n✓ Multi-scale processing complete!")
  message("  Comparison: ", comparison_dir)
  
  invisible(results)
}

create_comparison_report <- function(results, comparison_dir, metric, years) {
  
  # Load all results
  all_stats <- list()
  all_results <- list()
  
  for (grid_name in names(results)) {
    res <- results[[grid_name]]
    
    # Load overall stats
    stats_file <- fs::path(res$run_dir, "error_analysis", "overall_stats.csv")
    if (fs::file_exists(stats_file)) {
      stats <- readr::read_csv(stats_file, show_col_types = FALSE) |>
        dplyr::mutate(grid = grid_name)
      all_stats[[grid_name]] <- stats
    }
    
    # Load full results
    if (fs::file_exists(res$results_file)) {
      full_res <- readr::read_csv(res$results_file, show_col_types = FALSE) |>
        dplyr::mutate(grid = grid_name)
      all_results[[grid_name]] <- full_res
    }
  }
  
  stats_combined <- dplyr::bind_rows(all_stats)
  results_combined <- dplyr::bind_rows(all_results)
  
  # Write combined tables
  readr::write_csv(stats_combined, fs::path(comparison_dir, "overall_stats_by_grid.csv"))
  readr::write_csv(results_combined, fs::path(comparison_dir, "results_all_grids.csv"))
  
  message("  ✓ Wrote combined statistics")
  
  # Create comparison plots
  message("  → Creating comparison plots...")
  
  # 1. Error components by grid scale
  p1 <- stats_combined |>
    dplyr::select(grid, mean_se, mean_pos_sd, mean_total_sd) |>
    tidyr::pivot_longer(cols = c(mean_se, mean_pos_sd, mean_total_sd),
                        names_to = "error_type", values_to = "value") |>
    dplyr::mutate(
      error_type = factor(error_type,
                          levels = c("mean_se", "mean_pos_sd", "mean_total_sd"),
                          labels = c("Sampling Error", "Positional SD", "Total SD"))
    ) |>
    ggplot(aes(x = grid, y = value, fill = error_type)) +
    geom_col(position = "dodge") +
    scale_fill_viridis_d(option = "D", begin = 0.3, end = 0.9) +
    labs(title = "Error Components by Hex Grid Scale",
         subtitle = paste("Metric:", toupper(metric)),
         x = "Hex Grid Scale", y = "Mean SD (Mg/ha)", fill = "Error Type") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(fs::path(comparison_dir, "error_by_scale.png"), p1, 
         width = 10, height = 6, dpi = 300)
  
  # 2. Positional error fraction by grid scale
  p2 <- stats_combined |>
    ggplot(aes(x = grid, y = mean_pos_fraction * 100)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = paste0(round(mean_pos_fraction * 100, 1), "%")),
              vjust = -0.5, size = 4) +
    labs(title = "Positional Error Contribution by Grid Scale",
         subtitle = "Percentage of total uncertainty due to coordinate fuzzing",
         x = "Hex Grid Scale", y = "Positional Error (% of Total)") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(fs::path(comparison_dir, "positional_fraction_by_scale.png"), p2, 
         width = 10, height = 6, dpi = 300)
  
  # 3. Sample size distribution by grid scale
  p3 <- results_combined |>
    ggplot(aes(x = grid, y = n_plots, fill = grid)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, fill = "white", alpha = 0.8) +
    scale_fill_viridis_d() +
    scale_y_log10() +
    labs(title = "Sample Size Distribution by Grid Scale",
         subtitle = "Distribution of plots per hex across all hex-years",
         x = "Hex Grid Scale", y = "Plots per Hex (log scale)") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "none")
  
  ggsave(fs::path(comparison_dir, "sample_size_by_scale.png"), p3, 
         width = 10, height = 6, dpi = 300)
  
  # 4. Total uncertainty vs sample size by grid
  p4 <- results_combined |>
    ggplot(aes(x = n_plots, y = total_sd, color = grid)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1.5) +
    scale_color_viridis_d() +
    scale_x_log10() +
    labs(title = "Uncertainty vs Sample Size by Grid Scale",
         subtitle = "How total SD relates to plots per hex across grid sizes",
         x = "Plots per Hex (log scale)", y = "Total SD (Mg/ha)", color = "Grid") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(fs::path(comparison_dir, "uncertainty_vs_sample_by_scale.png"), p4, 
         width = 10, height = 6, dpi = 300)
  
  message("  ✓ Created ", length(list.files(comparison_dir, pattern = "\\.png$")), " comparison plots")
  
  # Create summary report
  report_lines <- c(
    "═══════════════════════════════════════════════════════════",
    "MULTI-SCALE COMPARISON REPORT",
    "═══════════════════════════════════════════════════════════",
    "",
    paste("Generated:", Sys.time()),
    paste("Metric:", toupper(metric)),
    paste("Years:", paste(years, collapse = ", ")),
    "",
    "───────────────────────────────────────────────────────────",
    "COMPARISON BY GRID SCALE",
    "───────────────────────────────────────────────────────────",
    ""
  )
  
  # Add table
  for (i in seq_len(nrow(stats_combined))) {
    row <- stats_combined[i, ]
    report_lines <- c(report_lines,
                      paste0("Grid: ", row$grid),
                      paste0("  Total hexes: ", row$n_hexes),
                      paste0("  Mean plots/hex: ", round(row$mean_plots_per_hex, 1)),
                      paste0("  Sampling Error (SE): ", round(row$mean_se, 3), " Mg/ha"),
                      paste0("  Positional SD: ", round(row$mean_pos_sd, 3), " Mg/ha"),
                      paste0("  Total SD: ", round(row$mean_total_sd, 3), " Mg/ha"),
                      paste0("  Positional fraction: ", round(100 * row$mean_pos_fraction, 1), "%"),
                      ""
    )
  }
  
  report_lines <- c(report_lines,
                    "───────────────────────────────────────────────────────────",
                    "KEY FINDINGS",
                    "───────────────────────────────────────────────────────────",
                    "",
                    "Scale Effects:",
                    paste0("  - Finest grid (", stats_combined$grid[1], "): ",
                           round(stats_combined$mean_total_sd[1], 3), " Mg/ha total SD"),
                    paste0("  - Coarsest grid (", stats_combined$grid[nrow(stats_combined)], "): ",
                           round(stats_combined$mean_total_sd[nrow(stats_combined)], 3), " Mg/ha total SD"),
                    "",
                    "Positional Error Trends:",
                    "  (See positional_fraction_by_scale.png)",
                    "",
                    "═══════════════════════════════════════════════════════════"
  )
  
  report_file <- fs::path(comparison_dir, "comparison_report.txt")
  writeLines(report_lines, report_file)
  message("  ✓ Wrote summary: ", report_file)
  
  invisible(comparison_dir)
}

# CLI interface
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
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
  
  skip_existing <- !("--reprocess" %in% args)
  
  process_all_scales(
    project_dir = cfg$project_dir %||% ".",
    metric = cfg$metric %||% "aglb",
    years = years_vec,
    level_window = cfg$level_window %||% 3,
    skip_existing = skip_existing
  )
}