# R/visualize_consolidated_results.R
# Comprehensive visualization of FIA error analysis and FIA vs NEFIN comparison

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(tidyr)
  library(scales); library(viridis); library(patchwork); library(fs)
})

visualize_results <- function(consolidated_dir = NULL,
                              output_dir = NULL) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Consolidated Results Visualization                      ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Find most recent consolidated directory if not specified
  if (is.null(consolidated_dir)) {
    all_dirs <- list.dirs("runs", recursive = FALSE)
    consol_dirs <- all_dirs[grepl("^runs/consolidated_", all_dirs)]
    
    if (length(consol_dirs) == 0) {
      stop("No consolidated results found. Run R/master_process_all.R first.")
    }
    
    # Get most recent
    consol_dirs <- consol_dirs[order(file.mtime(consol_dirs), decreasing = TRUE)]
    consolidated_dir <- consol_dirs[1]
  }
  
  cat("Using consolidated results from:\n")
  cat("  ", consolidated_dir, "\n\n")
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- fs::path(consolidated_dir, "visualizations")
  }
  fs::dir_create(output_dir, recurse = TRUE)
  
  # Load data
  cat("Loading data...\n")
  
  fia_file <- fs::path(consolidated_dir, "fia_all_scales.csv")
  nefin_file <- fs::path(consolidated_dir, "fia_nefin_comparison_all_scales.csv")
  fia_summary_file <- fs::path(consolidated_dir, "fia_summary_by_scale.csv")
  nefin_summary_file <- fs::path(consolidated_dir, "nefin_comparison_summary_by_scale.csv")
  
  if (!fs::file_exists(fia_file)) {
    stop("FIA results not found: ", fia_file)
  }
  
  fia <- readr::read_csv(fia_file, show_col_types = FALSE)
  cat("  ✓ FIA results: ", format(nrow(fia), big.mark = ","), " hex-year-scale observations\n")
  
  fia_summary <- if (fs::file_exists(fia_summary_file)) {
    readr::read_csv(fia_summary_file, show_col_types = FALSE)
  } else NULL
  
  has_nefin <- fs::file_exists(nefin_file)
  
  nefin <- if (has_nefin) {
    df <- readr::read_csv(nefin_file, show_col_types = FALSE)
    cat("  ✓ NEFIN comparison: ", format(nrow(df), big.mark = ","), " hex-year-scale observations\n")
    df
  } else {
    cat("  ⚠ NEFIN comparison not found, skipping NEFIN visualizations\n")
    NULL
  }
  
  nefin_summary <- if (has_nefin && fs::file_exists(nefin_summary_file)) {
    readr::read_csv(nefin_summary_file, show_col_types = FALSE)
  } else NULL
  
  cat("\n")
  
  # ========================================================================
  # SECTION 1: FIA ERROR ANALYSIS
  # ========================================================================
  
  cat("Creating FIA error visualizations...\n")
  
  # 1.1 Error decomposition by scale
  cat("  → Error decomposition by scale\n")
  
  error_long <- fia |>
    dplyr::select(grid_scale, se, positional_sd, total_sd) |>
    tidyr::pivot_longer(cols = c(se, positional_sd, total_sd),
                        names_to = "error_type", values_to = "value") |>
    dplyr::mutate(
      error_type = factor(error_type,
                          levels = c("se", "positional_sd", "total_sd"),
                          labels = c("Sampling Error", "Positional SD", "Total SD")),
      grid_scale = factor(grid_scale, levels = unique(fia$grid_scale))
    )
  
  p1 <- ggplot(error_long, aes(x = grid_scale, y = value, fill = error_type)) +
    geom_boxplot(outlier.alpha = 0.3, position = position_dodge(0.8)) +
    scale_fill_viridis_d(option = "D", begin = 0.3, end = 0.9) +
    labs(title = "Error Decomposition by Grid Scale",
         subtitle = "Comparison of sampling vs positional uncertainty",
         x = "Grid Scale", y = "Standard Deviation (Mg/ha)", fill = "Error Type") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom")
  
  ggsave(fs::path(output_dir, "01_error_decomposition_by_scale.png"), p1,
         width = 10, height = 6, dpi = 300)
  
  # 1.2 Positional error fraction by scale
  cat("  → Positional error fraction\n")
  
  p2 <- fia |>
    dplyr::mutate(
      pos_fraction = positional_sd / total_sd,
      grid_scale = factor(grid_scale, levels = unique(fia$grid_scale))
    ) |>
    ggplot(aes(x = grid_scale, y = pos_fraction * 100, fill = grid_scale)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.alpha = 0.3) +
    scale_fill_viridis_d(option = "C") +
    labs(title = "Positional Error Contribution by Grid Scale",
         subtitle = "Percentage of total uncertainty from coordinate fuzzing",
         x = "Grid Scale", y = "Positional Error (% of Total)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "none")
  
  ggsave(fs::path(output_dir, "02_positional_fraction_by_scale.png"), p2,
         width = 10, height = 6, dpi = 300)
  
  # 1.3 Sample size vs total error by scale
  cat("  → Sample size vs error\n")
  
  p3 <- fia |>
    dplyr::mutate(grid_scale = factor(grid_scale, levels = unique(fia$grid_scale))) |>
    ggplot(aes(x = n_plots, y = total_sd, color = grid_scale)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +
    scale_color_viridis_d(option = "D") +
    scale_x_log10(labels = comma) +
    labs(title = "Total Uncertainty vs Sample Size by Grid Scale",
         subtitle = "Relationship between plots per hex and total error",
         x = "Plots per Hex (log scale)", y = "Total SD (Mg/ha)", color = "Grid Scale") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom")
  
  ggsave(fs::path(output_dir, "03_error_vs_sample_size.png"), p3,
         width = 10, height = 6, dpi = 300)
  
  # 1.4 Summary statistics table plot
  if (!is.null(fia_summary)) {
    cat("  → Summary statistics by scale\n")
    
    p4 <- fia_summary |>
      dplyr::select(grid_scale, mean_sampling_se, mean_positional_sd, mean_total_sd, 
                    pos_fraction, mean_n_plots) |>
      dplyr::mutate(
        pos_fraction_pct = round(pos_fraction * 100, 1),
        across(c(mean_sampling_se, mean_positional_sd, mean_total_sd), ~round(., 3)),
        mean_n_plots = round(mean_n_plots, 1)
      ) |>
      dplyr::select(-pos_fraction) |>
      tidyr::pivot_longer(cols = -grid_scale, names_to = "metric", values_to = "value") |>
      dplyr::mutate(
        metric = factor(metric,
                        levels = c("mean_sampling_se", "mean_positional_sd", 
                                   "mean_total_sd", "pos_fraction_pct", "mean_n_plots"),
                        labels = c("Sampling SE", "Positional SD", "Total SD",
                                   "Pos % of Total", "Plots/Hex"))
      ) |>
      ggplot(aes(x = grid_scale, y = metric, fill = value)) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = value), color = "white", fontface = "bold", size = 4) +
      scale_fill_viridis_c(option = "A", direction = -1) +
      labs(title = "Summary Statistics by Grid Scale",
           x = "Grid Scale", y = NULL, fill = "Value") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 14),
            panel.grid = element_blank(),
            legend.position = "none")
    
    ggsave(fs::path(output_dir, "04_summary_heatmap.png"), p4,
           width = 10, height = 5, dpi = 300)
  }
  
  # 1.5 SE vs Positional SD scatter
  cat("  → SE vs Positional SD relationship\n")
  
  p5 <- fia |>
    ggplot(aes(x = se, y = positional_sd, color = grid_scale)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    scale_color_viridis_d(option = "D") +
    labs(title = "Sampling Error vs Positional Uncertainty",
         subtitle = "Dashed line = equal contribution",
         x = "Sampling Error (SE, Mg/ha)", y = "Positional SD (Mg/ha)",
         color = "Grid Scale") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom")
  
  ggsave(fs::path(output_dir, "05_se_vs_positional.png"), p5,
         width = 10, height = 8, dpi = 300)
  
  # ========================================================================
  # SECTION 2: FIA vs NEFIN COMPARISON
  # ========================================================================
  
  if (!is.null(nefin)) {
    cat("\nCreating FIA vs NEFIN visualizations...\n")
    
    # Filter to hexes with both datasets
    nefin_both <- nefin |> dplyr::filter(has_both)
    
    if (nrow(nefin_both) > 0) {
      
      # 2.1 Scatter plot by scale
      cat("  → FIA vs NEFIN scatter by scale\n")
      
      p6 <- nefin_both |>
        dplyr::mutate(grid_scale = factor(grid_scale, levels = unique(nefin$grid_scale))) |>
        ggplot(aes(x = mean_fia, y = mean_nefin)) +
        geom_point(alpha = 0.4, size = 1.5, color = "steelblue") +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
        geom_smooth(method = "lm", color = "darkblue", se = TRUE, linewidth = 1) +
        facet_wrap(~grid_scale, scales = "free", ncol = 2) +
        labs(title = "FIA vs NEFIN Biomass by Grid Scale",
             subtitle = "Red line = perfect agreement, Blue line = actual fit",
             x = "FIA AGLB (Mg/ha)", y = "NEFIN AGLB (Mg/ha)") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(face = "bold", size = 14))
      
      ggsave(fs::path(output_dir, "06_fia_vs_nefin_scatter.png"), p6,
             width = 12, height = 10, dpi = 300)
      
      # 2.2 Difference (bias) by scale
      cat("  → Bias distribution by scale\n")
      
      p7 <- nefin_both |>
        dplyr::mutate(grid_scale = factor(grid_scale, levels = unique(nefin$grid_scale))) |>
        ggplot(aes(x = grid_scale, y = diff, fill = grid_scale)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
        geom_violin(alpha = 0.7) +
        geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.alpha = 0.3) +
        scale_fill_viridis_d(option = "C") +
        labs(title = "Bias (NEFIN - FIA) by Grid Scale",
             subtitle = "Systematic differences between datasets",
             x = "Grid Scale", y = "Difference (Mg/ha)") +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(face = "bold", size = 14),
              legend.position = "none")
      
      ggsave(fs::path(output_dir, "07_bias_by_scale.png"), p7,
             width = 10, height = 6, dpi = 300)
      
      # 2.3 Agreement metrics by scale
      if (!is.null(nefin_summary)) {
        cat("  → Agreement metrics comparison\n")
        
        p8 <- nefin_summary |>
          dplyr::select(grid_scale, rmse, correlation, mean_diff) |>
          tidyr::pivot_longer(cols = -grid_scale, names_to = "metric", values_to = "value") |>
          dplyr::mutate(
            metric = factor(metric,
                            levels = c("correlation", "rmse", "mean_diff"),
                            labels = c("Correlation", "RMSE (Mg/ha)", "Mean Bias (Mg/ha)"))
          ) |>
          ggplot(aes(x = grid_scale, y = value, fill = grid_scale)) +
          geom_col() +
          geom_text(aes(label = round(value, 2)), vjust = -0.5, size = 3.5, fontface = "bold") +
          facet_wrap(~metric, scales = "free_y", ncol = 3) +
          scale_fill_viridis_d(option = "D") +
          labs(title = "FIA vs NEFIN Agreement Metrics by Grid Scale",
               x = "Grid Scale", y = "Value") +
          theme_minimal(base_size = 12) +
          theme(plot.title = element_text(face = "bold", size = 14),
                legend.position = "none",
                strip.text = element_text(face = "bold", size = 11))
        
        ggsave(fs::path(output_dir, "08_agreement_metrics.png"), p8,
               width = 12, height = 4, dpi = 300)
      }
      
      # 2.4 Coverage comparison
      cat("  → Dataset coverage comparison\n")
      
      coverage <- nefin |>
        dplyr::group_by(grid_scale) |>
        dplyr::summarise(
          both = sum(has_both),
          fia_only = sum(has_fia & !has_nefin),
          nefin_only = sum(has_nefin & !has_fia),
          .groups = "drop"
        ) |>
        tidyr::pivot_longer(cols = c(both, fia_only, nefin_only),
                            names_to = "coverage", values_to = "count") |>
        dplyr::mutate(
          coverage = factor(coverage,
                            levels = c("both", "fia_only", "nefin_only"),
                            labels = c("Both FIA & NEFIN", "FIA Only", "NEFIN Only"))
        )
      
      p9 <- ggplot(coverage, aes(x = grid_scale, y = count, fill = coverage)) +
        geom_col(position = "stack") +
        geom_text(aes(label = count), position = position_stack(vjust = 0.5),
                  color = "white", fontface = "bold", size = 3.5) +
        scale_fill_manual(values = c("Both FIA & NEFIN" = "#2E8B57",
                                     "FIA Only" = "#4682B4",
                                     "NEFIN Only" = "#DAA520")) +
        labs(title = "Data Coverage by Grid Scale",
             subtitle = "Number of hex-years with each dataset combination",
             x = "Grid Scale", y = "Number of Hex-Years", fill = "Coverage") +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(face = "bold", size = 14),
              legend.position = "bottom")
      
      ggsave(fs::path(output_dir, "09_coverage_by_scale.png"), p9,
             width = 10, height = 6, dpi = 300)
      
      # 2.5 Percent difference distribution
      cat("  → Percent difference analysis\n")
      
      p10 <- nefin_both |>
        dplyr::filter(is.finite(pct_diff), abs(pct_diff) < 100) |>  # Remove extreme outliers
        dplyr::mutate(grid_scale = factor(grid_scale, levels = unique(nefin$grid_scale))) |>
        ggplot(aes(x = pct_diff, fill = grid_scale)) +
        geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
        facet_wrap(~grid_scale, scales = "free_y", ncol = 2) +
        scale_fill_viridis_d(option = "C") +
        labs(title = "Percent Difference Distribution by Grid Scale",
             subtitle = "Percent difference = 100 × (NEFIN - FIA) / FIA",
             x = "Percent Difference (%)", y = "Count") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(face = "bold", size = 14),
              legend.position = "none")
      
      ggsave(fs::path(output_dir, "10_percent_diff_distribution.png"), p10,
             width = 12, height = 10, dpi = 300)
      
    } else {
      cat("  ⚠ No hexes with both FIA and NEFIN data\n")
    }
  }
  
  # ========================================================================
  # SECTION 3: COMBINED DASHBOARD
  # ========================================================================
  
  cat("\nCreating combined dashboard...\n")
  
  # Create a multi-panel summary figure
  p_error_box <- fia |>
    dplyr::select(grid_scale, se, positional_sd, total_sd) |>
    tidyr::pivot_longer(cols = c(se, positional_sd, total_sd),
                        names_to = "error_type", values_to = "value") |>
    dplyr::mutate(
      error_type = factor(error_type,
                          levels = c("se", "positional_sd", "total_sd"),
                          labels = c("SE", "Pos SD", "Total")),
      grid_scale = factor(grid_scale, levels = unique(fia$grid_scale))
    ) |>
    ggplot(aes(x = error_type, y = value, fill = error_type)) +
    geom_boxplot(outlier.alpha = 0.2) +
    facet_wrap(~grid_scale, ncol = 4) +
    scale_fill_viridis_d(option = "D", begin = 0.3, end = 0.9) +
    labs(title = "A. Error Decomposition", x = NULL, y = "SD (Mg/ha)") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
  
  p_pos_frac <- fia |>
    dplyr::mutate(
      pos_fraction = positional_sd / total_sd * 100,
      grid_scale = factor(grid_scale, levels = unique(fia$grid_scale))
    ) |>
    ggplot(aes(x = grid_scale, y = pos_fraction, fill = grid_scale)) +
    geom_boxplot(outlier.alpha = 0.2) +
    scale_fill_viridis_d(option = "C") +
    labs(title = "B. Positional Error Fraction", x = NULL, y = "% of Total") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
  
  if (!is.null(nefin) && nrow(nefin_both) > 0) {
    p_nefin_scatter <- nefin_both |>
      ggplot(aes(x = mean_fia, y = mean_nefin, color = grid_scale)) +
      geom_point(alpha = 0.3, size = 0.5) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      scale_color_viridis_d(option = "D") +
      labs(title = "C. FIA vs NEFIN", x = "FIA (Mg/ha)", y = "NEFIN (Mg/ha)") +
      theme_minimal(base_size = 9) +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            plot.title = element_text(face = "bold"))
    
    p_nefin_bias <- nefin_both |>
      dplyr::mutate(grid_scale = factor(grid_scale, levels = unique(nefin$grid_scale))) |>
      ggplot(aes(x = grid_scale, y = diff, fill = grid_scale)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_violin(alpha = 0.7) +
      scale_fill_viridis_d(option = "C") +
      labs(title = "D. Bias (NEFIN - FIA)", x = NULL, y = "Difference (Mg/ha)") +
      theme_minimal(base_size = 9) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold"))
    
    dashboard <- (p_error_box / p_pos_frac / (p_nefin_scatter | p_nefin_bias)) +
      plot_annotation(
        title = "FIA Error Analysis & NEFIN Comparison Dashboard",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
      )
  } else {
    dashboard <- (p_error_box / p_pos_frac) +
      plot_annotation(
        title = "FIA Error Analysis Dashboard",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
      )
  }
  
  ggsave(fs::path(output_dir, "00_dashboard.png"), dashboard,
         width = 14, height = if (!is.null(nefin)) 12 else 8, dpi = 300)
  
  # ========================================================================
  # DONE
  # ========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  VISUALIZATION COMPLETE                                   ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  cat("Output directory:\n")
  cat("  ", output_dir, "\n")
  cat("\n")
  cat("Created visualizations:\n")
  
  viz_files <- list.files(output_dir, pattern = "\\.png$", full.names = FALSE)
  for (f in viz_files) {
    cat("  ✓", f, "\n")
  }
  
  cat("\n")
  cat("Key figures:\n")
  cat("  00_dashboard.png - Overview of all results\n")
  cat("  01-05 - FIA error analysis\n")
  if (!is.null(nefin)) {
    cat("  06-10 - FIA vs NEFIN comparison\n")
  }
  cat("\n")
  
  invisible(output_dir)
}

# CLI
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  get_arg <- function(flag, default = NULL) {
    hit <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
  }
  
  visualize_results(
    consolidated_dir = get_arg("--dir", NULL),
    output_dir = get_arg("--out", NULL)
  )
}

