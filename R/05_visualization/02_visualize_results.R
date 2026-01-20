#!/usr/bin/env Rscript
# R/visualize_consolidated_results.R
# Create publication-quality visualizations from consolidated results
# UPDATED: Proper scale ordering with 64kha (formerly "fia")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(scales)
  library(viridis)
  library(patchwork)
  library(fs)
})

# Source utilities
if (file.exists("R/utils_scale_names.R")) {

  source("R/utils_scale_names.R")
} else {
  # Inline fallback
  standardize_scale_name <- function(x) {
    x <- as.character(x)
    x[tolower(x) == "fia"] <- "64kha"
    x
  }
  get_scale_order <- function() {
    c("100ha", "500ha", "1kha", "5kha", "10kha", "50kha", "64kha", "100kha")
  }
  order_scales <- function(df, scale_col = "grid_scale") {
    if (scale_col %in% names(df)) {
      df[[scale_col]] <- standardize_scale_name(df[[scale_col]])
      scale_order <- get_scale_order()
      existing <- unique(df[[scale_col]])
      ordered_levels <- scale_order[scale_order %in% existing]
      extra <- setdiff(existing, scale_order)
      if (length(extra) > 0) ordered_levels <- c(ordered_levels, sort(extra))
      df[[scale_col]] <- factor(df[[scale_col]], levels = ordered_levels, ordered = TRUE)
    }
    df
  }
}

visualize_results <- function(consolidated_dir = NULL,
                               output_dir = NULL) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Consolidated Results Visualization                      ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Find consolidated directory
  if (is.null(consolidated_dir)) {
    all_dirs <- list.dirs("runs", recursive = FALSE)
    consol_dirs <- all_dirs[grepl("consolidated_", all_dirs)]
    
    if (length(consol_dirs) == 0) {
      stop("No consolidated results found. Run master_process_all.R first.")
    }
    
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
  
  # Load FIA results
  cat("Loading data...\n")
  fia_file <- fs::path(consolidated_dir, "fia_all_scales.csv")
  
  if (!fs::file_exists(fia_file)) {
    stop("FIA results not found: ", fia_file)
  }
  
  fia <- readr::read_csv(fia_file, show_col_types = FALSE)
  
  # Standardize scale names and order
  fia <- order_scales(fia, "grid_scale")
  
  cat("  ✓ FIA results: ", format(nrow(fia), big.mark = ","), " hex-year-scale observations\n")
  
  # Load NEFIN comparison if available
  nefin_file <- fs::path(consolidated_dir, "fia_nefin_comparison_all_scales.csv")
  has_nefin <- fs::file_exists(nefin_file)
  
  nefin <- NULL
  if (has_nefin) {
    nefin <- readr::read_csv(nefin_file, show_col_types = FALSE)
    nefin <- order_scales(nefin, "grid_scale")
    cat("  ✓ NEFIN comparison: ", format(nrow(nefin), big.mark = ","), " hex-year-scale observations\n")
  } else {
    cat("  ⚠ NEFIN comparison not found, skipping NEFIN visualizations\n")
  }
  
  # ==========================================================================
  # FIA ERROR VISUALIZATIONS
  # ==========================================================================
  
  cat("\nCreating FIA error visualizations...\n")
  
  # Clean data for plotting
  fia_clean <- fia %>%
    filter(is.finite(se), is.finite(positional_sd), is.finite(total_sd), total_sd > 0)
  
  # 1. Error decomposition by scale
  cat("  → Error decomposition by scale\n")
  
  error_summary <- fia_clean %>%
    group_by(grid_scale) %>%
    summarise(
      mean_se = mean(se, na.rm = TRUE),
      mean_pos_sd = mean(positional_sd, na.rm = TRUE),
      mean_total_sd = mean(total_sd, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(mean_se, mean_pos_sd),
                 names_to = "error_type", values_to = "value") %>%
    mutate(error_type = factor(error_type,
                               levels = c("mean_pos_sd", "mean_se"),
                               labels = c("Positional SD", "Sampling SE")))
  
  p1 <- ggplot(error_summary, aes(x = grid_scale, y = value, fill = error_type)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c("Positional SD" = "#E63946", "Sampling SE" = "#457B9D")) +
    labs(title = "Error Decomposition by Grid Scale",
         subtitle = "Stacked: Positional uncertainty + Sampling error",
         x = "Grid Scale (hectares)", y = "Standard Deviation (Mg/ha)",
         fill = "Error Component") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
  
  ggsave(fs::path(output_dir, "01_error_decomposition_by_scale.png"), p1,
         width = 10, height = 6, dpi = 300)
  
  # 2. Positional fraction by scale
  cat("  → Positional error fraction\n")
  
  pos_frac <- fia_clean %>%
    mutate(pos_fraction = positional_sd / total_sd) %>%
    group_by(grid_scale) %>%
    summarise(
      mean_frac = mean(pos_fraction, na.rm = TRUE),
      median_frac = median(pos_fraction, na.rm = TRUE),
      q25 = quantile(pos_fraction, 0.25, na.rm = TRUE),
      q75 = quantile(pos_fraction, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  p2 <- ggplot(pos_frac, aes(x = grid_scale, y = mean_frac * 100)) +
    geom_col(fill = "#E63946", alpha = 0.8) +
    geom_errorbar(aes(ymin = q25 * 100, ymax = q75 * 100), width = 0.3) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") +
    labs(title = "Positional Error Fraction by Grid Scale",
         subtitle = "Bars = mean, error bars = IQR | Dashed line = 50%",
         x = "Grid Scale", y = "Positional Fraction (%)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(fs::path(output_dir, "02_positional_fraction_by_scale.png"), p2,
         width = 10, height = 6, dpi = 300)
  
  # 3. Sample size vs error
  cat("  → Sample size vs error\n")
  
  p3 <- fia_clean %>%
    filter(n_plots >= 1) %>%
    ggplot(aes(x = n_plots, y = se, color = grid_scale)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_viridis_d(option = "D") +
    labs(title = "Sampling Error vs Plot Count",
         subtitle = "Log-log scale | Lines = LOESS smoothers",
         x = "Number of Plots (log scale)", y = "Sampling SE (Mg/ha, log scale)",
         color = "Grid Scale") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "right")
  
  ggsave(fs::path(output_dir, "03_error_vs_sample_size.png"), p3,
         width = 10, height = 6, dpi = 300)
  
  # 4. Summary heatmap
  cat("  → Summary statistics by scale\n")
  
  scale_summary <- fia_clean %>%
    group_by(grid_scale) %>%
    summarise(
      n_hexes = n_distinct(hex_id),
      mean_biomass = mean(mean, na.rm = TRUE),
      mean_se = mean(se, na.rm = TRUE),
      mean_pos_sd = mean(positional_sd, na.rm = TRUE),
      mean_plots = mean(n_plots, na.rm = TRUE),
      pos_fraction = mean(positional_sd / total_sd, na.rm = TRUE),
      .groups = "drop"
    )
  
  heatmap_data <- scale_summary %>%
    select(grid_scale, mean_se, mean_pos_sd, pos_fraction, mean_plots) %>%
    pivot_longer(-grid_scale, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric,
                           levels = c("mean_plots", "mean_se", "mean_pos_sd", "pos_fraction"),
                           labels = c("Mean Plots/Hex", "Mean SE", "Mean Pos SD", "Pos Fraction")))
  
  # Normalize for heatmap
  heatmap_data <- heatmap_data %>%
    group_by(metric) %>%
    mutate(value_scaled = (value - min(value)) / (max(value) - min(value) + 0.001)) %>%
    ungroup()
  
  p4 <- ggplot(heatmap_data, aes(x = grid_scale, y = metric, fill = value_scaled)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f", value)), color = "white", size = 3.5) +
    scale_fill_viridis_c(option = "D", direction = -1) +
    labs(title = "Summary Statistics by Grid Scale",
         x = "Grid Scale", y = "") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  ggsave(fs::path(output_dir, "04_summary_heatmap.png"), p4,
         width = 10, height = 5, dpi = 300)
  
  # 5. SE vs Positional SD
  cat("  → SE vs Positional SD relationship\n")
  
  p5 <- fia_clean %>%
    ggplot(aes(x = positional_sd, y = se, color = grid_scale)) +
    geom_point(alpha = 0.2, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_color_viridis_d(option = "D") +
    coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100)) +
    labs(title = "Sampling Error vs Positional Uncertainty",
         subtitle = "Dashed line = equality | Points below line: sampling error dominates",
         x = "Positional SD (Mg/ha)", y = "Sampling SE (Mg/ha)",
         color = "Grid Scale") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "right")
  
  ggsave(fs::path(output_dir, "05_se_vs_positional.png"), p5,
         width = 10, height = 8, dpi = 300)
  
  # ==========================================================================
  # NEFIN COMPARISON VISUALIZATIONS
  # ==========================================================================
  
  if (has_nefin) {
    cat("\nCreating FIA vs NEFIN visualizations...\n")
    
    # Filter to hexes with both FIA and NEFIN
    nefin_both <- nefin %>%
      filter(has_both == TRUE)
    
    if (nrow(nefin_both) > 0) {
      
      # 6. FIA vs NEFIN scatter
      cat("  → FIA vs NEFIN scatter by scale\n")
      
      p6 <- nefin_both %>%
        filter(is.finite(mean_fia), is.finite(mean_nefin)) %>%
        ggplot(aes(x = mean_fia, y = mean_nefin, color = grid_scale)) +
        geom_point(alpha = 0.5, size = 2) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
        facet_wrap(~grid_scale, ncol = 4) +
        scale_color_viridis_d(option = "D") +
        labs(title = "FIA vs NEFIN Biomass Estimates",
             subtitle = "Dashed = 1:1 line | Solid = linear fit",
             x = "FIA Mean (Mg/ha)", y = "NEFIN Mean (Mg/ha)") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(face = "bold", size = 14),
              legend.position = "none")
      
      ggsave(fs::path(output_dir, "06_fia_vs_nefin_scatter.png"), p6,
             width = 14, height = 8, dpi = 300)
      
      # 7. Bias distribution
      cat("  → Bias distribution by scale\n")
      
      p7 <- nefin_both %>%
        filter(is.finite(diff)) %>%
        ggplot(aes(x = diff, fill = grid_scale)) +
        geom_histogram(bins = 30, alpha = 0.7, color = "white") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
        facet_wrap(~grid_scale, ncol = 4, scales = "free_y") +
        scale_fill_viridis_d(option = "C") +
        labs(title = "Distribution of FIA-NEFIN Differences",
             subtitle = "NEFIN - FIA (Mg/ha) | Dashed line = zero bias",
             x = "Difference (Mg/ha)", y = "Count") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(face = "bold", size = 14),
              legend.position = "none")
      
      ggsave(fs::path(output_dir, "07_bias_by_scale.png"), p7,
             width = 14, height = 8, dpi = 300)
      
      # 8. Agreement metrics by scale
      cat("  → Dataset coverage comparison\n")
      
      coverage <- nefin %>%
        group_by(grid_scale) %>%
        summarise(
          total = n(),
          fia_only = sum(has_fia & !has_nefin, na.rm = TRUE),
          nefin_only = sum(!has_fia & has_nefin, na.rm = TRUE),
          both = sum(has_both, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        pivot_longer(cols = c(fia_only, nefin_only, both),
                     names_to = "coverage", values_to = "count") %>%
        mutate(coverage = factor(coverage,
                                 levels = c("fia_only", "both", "nefin_only"),
                                 labels = c("FIA Only", "Both", "NEFIN Only")))
      
      p8 <- ggplot(coverage, aes(x = grid_scale, y = count, fill = coverage)) +
        geom_col(position = "stack") +
        scale_fill_manual(values = c("FIA Only" = "#457B9D", "Both" = "#2A9D8F", "NEFIN Only" = "#E9C46A")) +
        labs(title = "Dataset Coverage by Grid Scale",
             x = "Grid Scale", y = "Number of Hexes",
             fill = "Coverage") +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(face = "bold", size = 14),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")
      
      ggsave(fs::path(output_dir, "08_agreement_metrics.png"), p8,
             width = 10, height = 6, dpi = 300)
      
      # 9. Coverage by scale (bar chart)
      cat("  → Percent difference analysis\n")
      
      coverage_pct <- nefin %>%
        group_by(grid_scale) %>%
        summarise(
          pct_both = 100 * sum(has_both, na.rm = TRUE) / n(),
          .groups = "drop"
        )
      
      p9 <- ggplot(coverage_pct, aes(x = grid_scale, y = pct_both, fill = grid_scale)) +
        geom_col() +
        geom_text(aes(label = sprintf("%.1f%%", pct_both)), vjust = -0.5, size = 3.5) +
        scale_fill_viridis_d(option = "D") +
        labs(title = "Percent of Hexes with Both FIA and NEFIN Data",
             x = "Grid Scale", y = "Percent (%)") +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(face = "bold", size = 14),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none") +
        coord_cartesian(ylim = c(0, max(coverage_pct$pct_both, na.rm = TRUE) * 1.15))
      
      ggsave(fs::path(output_dir, "09_coverage_by_scale.png"), p9,
             width = 10, height = 6, dpi = 300)
      
      # 10. Percent difference distribution
      if ("pct_diff" %in% names(nefin_both)) {
        p10 <- nefin_both %>%
          filter(is.finite(pct_diff), abs(pct_diff) < 200) %>%
          ggplot(aes(x = pct_diff, fill = grid_scale)) +
          geom_histogram(bins = 40, alpha = 0.7, color = "white") +
          geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
          facet_wrap(~grid_scale, ncol = 4, scales = "free_y") +
          scale_fill_viridis_d(option = "C") +
          labs(title = "Percent Difference Distribution (NEFIN vs FIA)",
               subtitle = "100 × (NEFIN - FIA) / FIA",
               x = "Percent Difference (%)", y = "Count") +
          theme_minimal(base_size = 11) +
          theme(plot.title = element_text(face = "bold", size = 14),
                legend.position = "none")
        
        ggsave(fs::path(output_dir, "10_percent_diff_distribution.png"), p10,
               width = 14, height = 8, dpi = 300)
      }
    }
  }
  
  # ==========================================================================
  # DASHBOARD
  # ==========================================================================
  
  cat("\nCreating combined dashboard...\n")
  
  # Simple 4-panel dashboard
  dashboard <- (p1 + p2) / (p3 + p5) +
    plot_annotation(
      title = "FIA Positional Uncertainty Analysis Dashboard",
      subtitle = paste("Consolidated results from", basename(consolidated_dir)),
      theme = theme(plot.title = element_text(face = "bold", size = 18),
                    plot.subtitle = element_text(size = 12, color = "gray40"))
    )
  
  ggsave(fs::path(output_dir, "00_dashboard.png"), dashboard,
         width = 16, height = 12, dpi = 300)
  
  # ==========================================================================
  # DONE
  # ==========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  VISUALIZATION COMPLETE                                   ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  cat("Output directory:\n")
  cat("  ", output_dir, "\n\n")
  
  figs <- list.files(output_dir, pattern = "\\.png$")
  cat("Created visualizations:\n")
  for (f in figs) {
    cat("  ✓", f, "\n")
  }
  
  cat("\nKey figures:\n")
  cat("  00_dashboard.png - Overview of all results\n")
  cat("  01-05 - FIA error analysis\n")
  if (has_nefin) {
    cat("  06-10 - FIA vs NEFIN comparison\n")
  }
  
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
