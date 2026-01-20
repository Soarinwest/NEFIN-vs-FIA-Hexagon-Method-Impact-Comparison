#!/usr/bin/env Rscript
# R/process_multi_period.R
# Process multiple time periods and compare results
# This allows analysis of how FIA fuzzing impact changes over time

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(fs)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
})

# Source utilities
if (file.exists("R/utils_scale_names.R")) source("R/utils_scale_names.R")

process_multi_period <- function(cfg_path = "configs/process.yml",
                                  periods = NULL,
                                  skip_jitter = TRUE,
                                  compare_periods = TRUE) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Multi-Period FIA-NEFIN Analysis                         ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Load config
  cfg <- yaml::read_yaml(cfg_path)
  
  # Get periods from config or argument
  if (is.null(periods)) {
    if (!is.null(cfg$analysis_periods)) {
      periods <- cfg$analysis_periods
    } else {
      # Default to single period from years
      periods <- list(
        list(name = paste0(min(cfg$years), "-", max(cfg$years)),
             years = cfg$years,
             label = "default")
      )
    }
  }
  
  cat("Processing", length(periods), "time periods:\n")
  for (p in periods) {
    cat("  -", p$name, ":", paste(range(p$years), collapse = "-"), "\n")
  }
  cat("\n")
  
  # Store results for each period
  all_results <- list()
  
  for (i in seq_along(periods)) {
    period <- periods[[i]]
    period_name <- period$name
    period_years <- period$years
    period_label <- period$label %||% period_name
    
    cat("═══════════════════════════════════════════════════════════\n")
    cat("Processing Period:", period_name, "(", min(period_years), "-", max(period_years), ")\n
")
    cat("═══════════════════════════════════════════════════════════\n")
    
    # Create temporary config for this period
    temp_cfg <- cfg
    temp_cfg$years <- period_years
    
    # Write temporary config
    temp_cfg_path <- paste0("configs/process_", period_label, ".yml")
    yaml::write_yaml(temp_cfg, temp_cfg_path)
    
    # Run pipeline stages for this period
    # Note: We use system() to call the pipeline to ensure clean environment
    
    cat("\n→ Running --all-scales for", period_name, "...\n")
    
    # Modify environment to use this config
    Sys.setenv(FIA_CONFIG = temp_cfg_path)
    
    # Source and run the processing
    # We'll call the key functions directly rather than system()
    
    if (file.exists("R/08_process_all_scales.R")) {
      source("R/08_process_all_scales.R", local = TRUE)
      
      # The function should pick up years from the sourced config
      # We need to temporarily modify the config file
      tryCatch({
        process_all_scales(
          cfg_path = temp_cfg_path,
          skip_jitter = skip_jitter
        )
      }, error = function(e) {
        cat("  ⚠ Error in process_all_scales:", e$message, "\n")
      })
    }
    
    cat("\n→ Running NEFIN comparison for", period_name, "...\n")
    
    if (file.exists("R/compare_fia_nefin.R")) {
      source("R/compare_fia_nefin.R", local = TRUE)
      
      tryCatch({
        run_all_scales(cfg_path = temp_cfg_path)
      }, error = function(e) {
        cat("  ⚠ Error in compare_fia_nefin:", e$message, "\n")
      })
    }
    
    # Store period info
    all_results[[period_label]] <- list(
      name = period_name,
      years = period_years,
      label = period_label,
      cfg_path = temp_cfg_path
    )
    
    cat("\n✓ Completed period:", period_name, "\n\n")
  }
  
  # Compare periods if requested
  if (compare_periods && length(periods) > 1) {
    cat("═══════════════════════════════════════════════════════════\n")
    cat("Comparing Periods\n")
    cat("═══════════════════════════════════════════════════════════\n")
    
    compare_period_results(all_results, cfg)
  }
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Multi-Period Processing Complete                        ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  invisible(all_results)
}

#' Compare results across multiple time periods
compare_period_results <- function(all_results, cfg) {
  
  # Find consolidated directories for each period
  runs_dir <- cfg$runs_dir %||% "runs"
  
  # Collect data from each period
  period_data <- list()
  
  for (label in names(all_results)) {
    result <- all_results[[label]]
    
    # Look for consolidated results
    # Pattern: consolidated_YYYYMMDD or similar
    consol_dirs <- list.dirs(runs_dir, recursive = FALSE)
    consol_dirs <- consol_dirs[grepl("consolidated_", consol_dirs)]
    
    if (length(consol_dirs) > 0) {
      # Get most recent
      consol_dir <- consol_dirs[order(file.mtime(consol_dirs), decreasing = TRUE)][1]
      
      fia_file <- fs::path(consol_dir, "fia_all_scales.csv")
      
      if (fs::file_exists(fia_file)) {
        df <- readr::read_csv(fia_file, show_col_types = FALSE)
        df$period <- result$name
        df$period_label <- label
        period_data[[label]] <- df
      }
    }
  }
  
  if (length(period_data) < 2) {
    cat("  ⚠ Need at least 2 periods with data to compare\n")
    return(invisible(NULL))
  }
  
  # Combine all periods
  combined <- dplyr::bind_rows(period_data)
  
  # Standardize scale names
  if (exists("standardize_scale_name")) {
    combined$grid_scale <- standardize_scale_name(combined$grid_scale)
  }
  
  # Create comparison output directory
  output_dir <- fs::path(runs_dir, "multi_period_comparison")
  fs::dir_create(output_dir, recurse = TRUE)
  
  # Summary by period and scale
  period_summary <- combined %>%
    filter(!is.na(se), is.finite(se), !is.na(positional_sd), is.finite(positional_sd)) %>%
    group_by(period, grid_scale) %>%
    summarise(
      n_hexes = n(),
      mean_biomass = mean(mean, na.rm = TRUE),
      mean_se = mean(se, na.rm = TRUE),
      mean_positional_sd = mean(positional_sd, na.rm = TRUE),
      mean_total_sd = mean(total_sd, na.rm = TRUE),
      pos_fraction = mean(positional_sd / total_sd, na.rm = TRUE),
      .groups = "drop"
    )
  
  readr::write_csv(period_summary, fs::path(output_dir, "period_comparison_summary.csv"))
  
  cat("\n→ Period Comparison Summary:\n")
  print(as.data.frame(period_summary), row.names = FALSE)
  
  # Create comparison plots
  cat("\n→ Creating comparison plots...\n")
  
  # Plot 1: Positional SD by period and scale
  p1 <- ggplot(period_summary, aes(x = grid_scale, y = mean_positional_sd, fill = period)) +
    geom_col(position = "dodge") +
    scale_fill_viridis_d(option = "D") +
    labs(title = "Positional Uncertainty by Period",
         x = "Grid Scale", y = "Mean Positional SD (Mg/ha)",
         fill = "Period") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Plot 2: Sampling error by period
  p2 <- ggplot(period_summary, aes(x = grid_scale, y = mean_se, fill = period)) +
    geom_col(position = "dodge") +
    scale_fill_viridis_d(option = "C") +
    labs(title = "Sampling Error by Period",
         x = "Grid Scale", y = "Mean SE (Mg/ha)",
         fill = "Period") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Plot 3: Positional fraction by period
  p3 <- ggplot(period_summary, aes(x = grid_scale, y = pos_fraction * 100, fill = period)) +
    geom_col(position = "dodge") +
    scale_fill_viridis_d(option = "B") +
    labs(title = "Positional Error Fraction by Period",
         x = "Grid Scale", y = "Positional Fraction (%)",
         fill = "Period") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Combined dashboard
  dashboard <- (p1 | p2) / p3 +
    plot_annotation(
      title = "Multi-Period Comparison: FIA Positional Uncertainty",
      subtitle = paste("Periods:", paste(unique(period_summary$period), collapse = " vs ")),
      theme = theme(plot.title = element_text(face = "bold", size = 16))
    )
  
  ggsave(fs::path(output_dir, "period_comparison_dashboard.png"), dashboard,
         width = 14, height = 10, dpi = 300)
  
  cat("  ✓ Saved: period_comparison_dashboard.png\n")
  cat("  ✓ Saved: period_comparison_summary.csv\n")
  cat("\nOutput directory:", output_dir, "\n")
  
  invisible(period_summary)
}

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Parse arguments
  cfg_path <- "configs/process.yml"
  skip_jitter <- TRUE
  compare <- TRUE
  
  for (arg in args) {
    if (grepl("^--config=", arg)) {
      cfg_path <- sub("^--config=", "", arg)
    }
    if (arg == "--rebuild-jitter") {
      skip_jitter <- FALSE
    }
    if (arg == "--no-compare") {
      compare <- FALSE
    }
  }
  
  if ("--help" %in% args || "-h" %in% args) {
    cat("\nUsage: Rscript R/process_multi_period.R [options]\n\n")
    cat("Options:\n")
    cat("  --config=PATH      Path to config file (default: configs/process.yml)\n")
    cat("  --rebuild-jitter   Rebuild jitter library (WARNING: slow!)\n")
    cat("  --no-compare       Skip period comparison\n")
    cat("  --help, -h         Show this help\n\n")
    quit(status = 0)
  }
  
  process_multi_period(
    cfg_path = cfg_path,
    skip_jitter = skip_jitter,
    compare_periods = compare
  )
}
