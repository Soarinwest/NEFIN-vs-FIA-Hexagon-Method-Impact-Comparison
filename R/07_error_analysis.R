# R/07_error_analysis.R
# Comprehensive error statistics and geospatial visualizations

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(ggplot2); library(fs)
  library(tidyr); library(scales); library(viridis)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

error_analysis <- function(run_dir, hex_path, hex_layer = NULL, 
                          create_maps = TRUE, create_plots = TRUE) {
  
  if (!fs::dir_exists(run_dir)) stop("Run directory not found: ", run_dir)
  
  message("→ Reading results from: ", run_dir)
  
  # Find results file
  res_files <- list.files(run_dir, pattern = "hex_.*_results\\.csv$", full.names = TRUE)
  if (!length(res_files)) stop("No results files found in: ", run_dir)
  
  results <- readr::read_csv(res_files[1], show_col_types = FALSE)
  message("  Loaded ", nrow(results), " hex-year observations")
  
  # Verify required columns
  req_cols <- c("hex_id", "mean", "se", "positional_sd", "total_sd", "n_plots")
  missing <- setdiff(req_cols, names(results))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  
  # Create output directories
  stats_dir <- fs::path(run_dir, "error_analysis")
  maps_dir <- fs::path(stats_dir, "maps")
  plots_dir <- fs::path(stats_dir, "plots")
  
  fs::dir_create(stats_dir, recurse = TRUE)
  if (create_maps) fs::dir_create(maps_dir, recurse = TRUE)
  if (create_plots) fs::dir_create(plots_dir, recurse = TRUE)
  
  # ========================================================================
  # PART 1: ERROR STATISTICS
  # ========================================================================
  
  message("\n→ Computing error statistics...")
  
  # Overall statistics
  overall_stats <- results |>
    dplyr::summarise(
      n_hexes = dplyr::n_distinct(hex_id),
      n_years = dplyr::n_distinct(year_label),
      total_obs = dplyr::n(),
      
      # Mean metric
      mean_value = mean(mean, na.rm = TRUE),
      sd_value = sd(mean, na.rm = TRUE),
      
      # Sampling error
      mean_se = mean(se, na.rm = TRUE),
      median_se = median(se, na.rm = TRUE),
      
      # Positional error
      mean_pos_sd = mean(positional_sd, na.rm = TRUE),
      median_pos_sd = median(positional_sd, na.rm = TRUE),
      
      # Total error
      mean_total_sd = mean(total_sd, na.rm = TRUE),
      median_total_sd = median(total_sd, na.rm = TRUE),
      
      # Error ratios
      mean_pos_fraction = mean(positional_sd / total_sd, na.rm = TRUE),
      median_pos_fraction = median(positional_sd / total_sd, na.rm = TRUE),
      
      # Sample sizes
      mean_plots_per_hex = mean(n_plots, na.rm = TRUE),
      median_plots_per_hex = median(n_plots, na.rm = TRUE)
    )
  
  # Per-hex statistics (across years)
  hex_stats <- results |>
    dplyr::group_by(hex_id) |>
    dplyr::summarise(
      n_years = dplyr::n(),
      mean_value = mean(mean, na.rm = TRUE),
      mean_se = mean(se, na.rm = TRUE),
      mean_pos_sd = mean(positional_sd, na.rm = TRUE),
      mean_total_sd = mean(total_sd, na.rm = TRUE),
      pos_fraction = mean(positional_sd / total_sd, na.rm = TRUE),
      total_plots = sum(n_plots, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Per-year statistics
  year_stats <- results |>
    dplyr::group_by(year_label) |>
    dplyr::summarise(
      n_hexes = dplyr::n(),
      mean_value = mean(mean, na.rm = TRUE),
      mean_se = mean(se, na.rm = TRUE),
      mean_pos_sd = mean(positional_sd, na.rm = TRUE),
      mean_total_sd = mean(total_sd, na.rm = TRUE),
      pos_fraction = mean(positional_sd / total_sd, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Error decomposition quantiles
  error_quantiles <- results |>
    dplyr::reframe(
      source = c("Sampling Error (SE)", "Positional SD", "Total SD"),
      q05 = c(quantile(se, 0.05, na.rm = TRUE),
              quantile(positional_sd, 0.05, na.rm = TRUE),
              quantile(total_sd, 0.05, na.rm = TRUE)),
      q25 = c(quantile(se, 0.25, na.rm = TRUE),
              quantile(positional_sd, 0.25, na.rm = TRUE),
              quantile(total_sd, 0.25, na.rm = TRUE)),
      median = c(median(se, na.rm = TRUE),
                 median(positional_sd, na.rm = TRUE),
                 median(total_sd, na.rm = TRUE)),
      q75 = c(quantile(se, 0.75, na.rm = TRUE),
              quantile(positional_sd, 0.75, na.rm = TRUE),
              quantile(total_sd, 0.75, na.rm = TRUE)),
      q95 = c(quantile(se, 0.95, na.rm = TRUE),
              quantile(positional_sd, 0.95, na.rm = TRUE),
              quantile(total_sd, 0.95, na.rm = TRUE))
    )
  
  # Write statistics
  readr::write_csv(overall_stats, fs::path(stats_dir, "overall_stats.csv"))
  readr::write_csv(hex_stats, fs::path(stats_dir, "hex_stats.csv"))
  readr::write_csv(year_stats, fs::path(stats_dir, "year_stats.csv"))
  readr::write_csv(error_quantiles, fs::path(stats_dir, "error_quantiles.csv"))
  
  message("  ✓ Wrote statistics to: ", stats_dir)
  
  # ========================================================================
  # PART 2: DIAGNOSTIC PLOTS
  # ========================================================================
  
  if (create_plots) {
    message("\n→ Creating diagnostic plots...")
    
    # 1. Error decomposition boxplot
    error_long <- results |>
      dplyr::select(hex_id, se, positional_sd, total_sd) |>
      tidyr::pivot_longer(cols = c(se, positional_sd, total_sd),
                          names_to = "error_type", values_to = "value") |>
      dplyr::mutate(error_type = factor(error_type, 
                                        levels = c("se", "positional_sd", "total_sd"),
                                        labels = c("Sampling Error", "Positional SD", "Total SD")))
    
    p1 <- ggplot(error_long, aes(x = error_type, y = value, fill = error_type)) +
      geom_boxplot(outlier.alpha = 0.3) +
      scale_fill_viridis_d(option = "D", begin = 0.3, end = 0.9) +
      labs(title = "Error Decomposition",
           subtitle = "Distribution of error components across all hex-years",
           x = NULL, y = "Standard Deviation (Mg/ha)") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold"))
    
    ggsave(fs::path(plots_dir, "error_decomposition_boxplot.png"), p1, 
           width = 8, height = 6, dpi = 300)
    
    # 2. Error vs sample size
    p2 <- ggplot(results, aes(x = n_plots, y = total_sd)) +
      geom_hex(bins = 50) +
      geom_smooth(method = "loess", color = "red", se = FALSE) +
      scale_fill_viridis_c(option = "C", trans = "log10") +
      scale_x_log10() +
      labs(title = "Total Error vs Sample Size",
           subtitle = "Relationship between plots per hex and uncertainty",
           x = "Number of Plots (log scale)",
           y = "Total SD (Mg/ha)",
           fill = "Count") +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(fs::path(plots_dir, "error_vs_sample_size.png"), p2, 
           width = 9, height = 6, dpi = 300)
    
    # 3. Positional error fraction
    p3 <- ggplot(results, aes(x = positional_sd / total_sd)) +
      geom_histogram(bins = 50, fill = "steelblue", color = "white") +
      geom_vline(xintercept = median(results$positional_sd / results$total_sd, na.rm = TRUE),
                 color = "red", linetype = "dashed", linewidth = 1) +
      labs(title = "Positional Error Contribution",
           subtitle = paste0("Median: ", 
                            round(100 * median(results$positional_sd / results$total_sd, na.rm = TRUE), 1),
                            "% of total error"),
           x = "Positional SD / Total SD",
           y = "Count") +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(fs::path(plots_dir, "positional_fraction_hist.png"), p3, 
           width = 8, height = 6, dpi = 300)
    
    # 4. Time series (if multiple years)
    if (length(unique(results$year_label)) > 1) {
      p4 <- ggplot(year_stats, aes(x = year_label)) +
        geom_line(aes(y = mean_se, color = "Sampling Error"), linewidth = 1) +
        geom_line(aes(y = mean_pos_sd, color = "Positional SD"), linewidth = 1) +
        geom_line(aes(y = mean_total_sd, color = "Total SD"), linewidth = 1.2) +
        geom_point(aes(y = mean_se, color = "Sampling Error"), size = 3) +
        geom_point(aes(y = mean_pos_sd, color = "Positional SD"), size = 3) +
        geom_point(aes(y = mean_total_sd, color = "Total SD"), size = 3) +
        scale_color_manual(values = c("Sampling Error" = "#440154",
                                      "Positional SD" = "#31688e",
                                      "Total SD" = "#35b779")) +
        labs(title = "Error Components Over Time",
             subtitle = "Mean error across all hexes by year",
             x = "Year", y = "Standard Deviation (Mg/ha)", color = "Error Type") +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(face = "bold"),
              legend.position = "bottom")
      
      ggsave(fs::path(plots_dir, "error_timeseries.png"), p4, 
             width = 10, height = 6, dpi = 300)
    }
    
    # 5. Scatter: SE vs Positional SD
    p5 <- ggplot(results, aes(x = se, y = positional_sd)) +
      geom_hex(bins = 50) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      scale_fill_viridis_c(option = "C", trans = "log10") +
      labs(title = "Sampling Error vs Positional Uncertainty",
           subtitle = "Dashed line: equal contribution",
           x = "Sampling Error (SE)",
           y = "Positional SD",
           fill = "Count") +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(fs::path(plots_dir, "se_vs_positional.png"), p5, 
           width = 8, height = 7, dpi = 300)
    
    message("  ✓ Created ", length(list.files(plots_dir)), " diagnostic plots")
  }
  
  # ========================================================================
  # PART 3: GEOSPATIAL MAPS
  # ========================================================================
  
  if (create_maps) {
    message("\n→ Creating geospatial maps...")
    
    # Load hex grid
    hx <- if (is.null(hex_layer)) {
      sf::st_read(hex_path, quiet = TRUE)
    } else {
      sf::st_read(hex_path, layer = hex_layer, quiet = TRUE)
    }
    
    if (!("hex_id" %in% names(hx))) {
      if ("ID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = ID)
      else if ("OBJECTID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = OBJECTID)
      else hx$hex_id <- seq_len(nrow(hx))
    }
    hx$hex_id <- as.character(hx$hex_id)
    
    # Join with hex stats
    hx_stats <- hx |>
      dplyr::left_join(hex_stats, by = "hex_id")
    
    # Helper function for maps
    make_hex_map <- function(data, var, title, filename, palette = "viridis") {
      
      # Remove infinite/NA values
      data_clean <- data |>
        dplyr::filter(is.finite(.data[[var]]))
      
      if (nrow(data_clean) == 0) {
        message("    Skipping ", filename, " - no valid data")
        return(invisible(NULL))
      }
      
      p <- ggplot(data_clean) +
        geom_sf(aes(fill = .data[[var]]), color = "grey70", linewidth = 0.1) +
        scale_fill_viridis_c(option = palette, na.value = "grey90") +
        labs(title = title, fill = NULL) +
        theme_void(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              legend.position = "right")
      
      ggsave(fs::path(maps_dir, filename), p, width = 10, height = 8, dpi = 300)
    }
    
    # Mean value
    make_hex_map(hx_stats, "mean_value", 
                 "Mean Aboveground Biomass (Mg/ha)", 
                 "mean_biomass.png", palette = "viridis")
    
    # Sampling error
    make_hex_map(hx_stats, "mean_se", 
                 "Mean Sampling Error (SE)", 
                 "sampling_error.png", palette = "plasma")
    
    # Positional SD
    make_hex_map(hx_stats, "mean_pos_sd", 
                 "Mean Positional Uncertainty (SD)", 
                 "positional_sd.png", palette = "inferno")
    
    # Total SD
    make_hex_map(hx_stats, "mean_total_sd", 
                 "Mean Total Uncertainty (SD)", 
                 "total_sd.png", palette = "magma")
    
    # Positional fraction
    make_hex_map(hx_stats, "pos_fraction", 
                 "Positional Error Fraction (Pos SD / Total SD)", 
                 "positional_fraction.png", palette = "cividis")
    
    # Sample size
    make_hex_map(hx_stats, "total_plots", 
                 "Total Plot Count", 
                 "sample_size.png", palette = "mako")
    
    # Create per-year maps if multiple years
    years <- unique(results$year_label)
    if (length(years) > 1) {
      for (yr in years) {
        yr_data <- results |> dplyr::filter(year_label == yr)
        hx_yr <- hx |> dplyr::left_join(yr_data, by = "hex_id")
        
        make_hex_map(hx_yr, "mean", 
                     paste0("Biomass (Mg/ha) - ", yr),
                     paste0("biomass_", yr, ".png"), palette = "viridis")
        
        make_hex_map(hx_yr, "positional_sd", 
                     paste0("Positional SD - ", yr),
                     paste0("positional_sd_", yr, ".png"), palette = "inferno")
      }
    }
    
    message("  ✓ Created ", length(list.files(maps_dir)), " geospatial maps")
  }
  
  # ========================================================================
  # SUMMARY REPORT
  # ========================================================================
  
  message("\n→ Writing summary report...")
  
  report_file <- fs::path(stats_dir, "error_analysis_report.txt")
  
  report_lines <- c(
    "═══════════════════════════════════════════════════════════════",
    "FIA ERROR ANALYSIS REPORT",
    "═══════════════════════════════════════════════════════════════",
    "",
    paste("Generated:", Sys.time()),
    paste("Run directory:", run_dir),
    "",
    "───────────────────────────────────────────────────────────────",
    "OVERALL STATISTICS",
    "───────────────────────────────────────────────────────────────",
    "",
    paste("Total hexes:", overall_stats$n_hexes),
    paste("Total years:", overall_stats$n_years),
    paste("Total observations:", overall_stats$total_obs),
    "",
    "METRIC VALUES:",
    paste("  Mean:", round(overall_stats$mean_value, 2), "Mg/ha"),
    paste("  SD:", round(overall_stats$sd_value, 2), "Mg/ha"),
    "",
    "ERROR COMPONENTS (mean across all hex-years):",
    paste("  Sampling Error (SE):", round(overall_stats$mean_se, 3), "Mg/ha"),
    paste("  Positional SD:", round(overall_stats$mean_pos_sd, 3), "Mg/ha"),
    paste("  Total SD:", round(overall_stats$mean_total_sd, 3), "Mg/ha"),
    "",
    "ERROR COMPONENTS (median):",
    paste("  Sampling Error (SE):", round(overall_stats$median_se, 3), "Mg/ha"),
    paste("  Positional SD:", round(overall_stats$median_pos_sd, 3), "Mg/ha"),
    paste("  Total SD:", round(overall_stats$median_total_sd, 3), "Mg/ha"),
    "",
    "POSITIONAL ERROR CONTRIBUTION:",
    paste("  Mean fraction of total error:", 
          round(100 * overall_stats$mean_pos_fraction, 1), "%"),
    paste("  Median fraction of total error:", 
          round(100 * overall_stats$median_pos_fraction, 1), "%"),
    "",
    "SAMPLE SIZE:",
    paste("  Mean plots per hex:", round(overall_stats$mean_plots_per_hex, 1)),
    paste("  Median plots per hex:", round(overall_stats$median_plots_per_hex, 1)),
    "",
    "───────────────────────────────────────────────────────────────",
    "ERROR QUANTILES",
    "───────────────────────────────────────────────────────────────",
    ""
  )
  
  # Add quantile table
  for (i in seq_len(nrow(error_quantiles))) {
    report_lines <- c(report_lines,
      paste0(error_quantiles$source[i], ":"),
      paste0("  5%: ", round(error_quantiles$q05[i], 3)),
      paste0("  25%: ", round(error_quantiles$q25[i], 3)),
      paste0("  50%: ", round(error_quantiles$median[i], 3)),
      paste0("  75%: ", round(error_quantiles$q75[i], 3)),
      paste0("  95%: ", round(error_quantiles$q95[i], 3)),
      ""
    )
  }
  
  report_lines <- c(report_lines,
    "───────────────────────────────────────────────────────────────",
    "OUTPUTS",
    "───────────────────────────────────────────────────────────────",
    "",
    "Statistics:",
    paste("  ", fs::path(stats_dir, "overall_stats.csv")),
    paste("  ", fs::path(stats_dir, "hex_stats.csv")),
    paste("  ", fs::path(stats_dir, "year_stats.csv")),
    paste("  ", fs::path(stats_dir, "error_quantiles.csv")),
    ""
  )
  
  if (create_plots) {
    report_lines <- c(report_lines,
      "Diagnostic Plots:",
      paste("  ", plots_dir)
    )
  }
  
  if (create_maps) {
    report_lines <- c(report_lines,
      "Geospatial Maps:",
      paste("  ", maps_dir)
    )
  }
  
  report_lines <- c(report_lines,
    "",
    "═══════════════════════════════════════════════════════════════"
  )
  
  writeLines(report_lines, report_file)
  
  message("\n✓ Error analysis complete!")
  message("  Report: ", report_file)
  if (create_plots) message("  Plots: ", plots_dir)
  if (create_maps) message("  Maps: ", maps_dir)
  
  invisible(list(
    stats = overall_stats,
    report = report_file,
    plots_dir = if (create_plots) plots_dir else NULL,
    maps_dir = if (create_maps) maps_dir else NULL
  ))
}

# CLI interface
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  get_arg <- function(flag, default = NULL) {
    hit <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
  }
  
  run_dir <- get_arg("--run", stop("Usage: Rscript R/07_error_analysis.R --run=runs/YOUR_RUN"))
  hex_path <- get_arg("--hex", "data/hex/hex_grid.geojson")
  hex_layer <- get_arg("--layer", NULL)
  
  no_maps <- "--no-maps" %in% args
  no_plots <- "--no-plots" %in% args
  
  error_analysis(
    run_dir = run_dir,
    hex_path = hex_path,
    hex_layer = hex_layer,
    create_maps = !no_maps,
    create_plots = !no_plots
  )
}