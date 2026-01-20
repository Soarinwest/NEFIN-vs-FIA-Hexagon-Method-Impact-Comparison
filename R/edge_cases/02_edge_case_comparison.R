#!/usr/bin/env Rscript
# =============================================================================
# 02_edge_case_comparison.R
# 
# Compare FIA vs NEFIN distributions for large-tree and mortality metrics
# Focus on tails and extreme conditions
#
# CRITICAL FRAMING:
#   - This is DESCRIPTIVE, not inferential
#   - Differences reflect SAMPLING DESIGN, not population differences
#   - NEFIN is biased toward public/research lands BY DESIGN
#
# Author: Soren Donisvitch
# Date: December 2024
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(fs)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

config <- list(
  # Input
  metrics_path = "data/processed/edge_case_metrics/plot_edge_case_metrics.csv",
  
  # Output
  output_dir = "outputs/edge_cases",
  figures_dir = "outputs/edge_cases/figures",
  tables_dir = "outputs/edge_cases/tables"
)

# Theme for all plots
theme_edge <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.1)),
      plot.subtitle = element_text(color = "gray50", size = rel(0.9)),
      plot.caption = element_text(color = "gray60", hjust = 0, size = rel(0.8)),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}

# Color palette
colors <- c("FIA" = "#e74c3c", "NEFIN" = "#27ae60")

# Standard caption for all figures
bias_caption <- "Note: Differences reflect sampling design and plot placement, not population inference.\nNEFIN over-samples public/research lands by design."

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

# Compute quantile comparison table
compute_quantile_comparison <- function(df, metric, quantiles = c(0.50, 0.75, 0.90, 0.95, 0.99)) {
  df %>%
    group_by(source) %>%
    summarize(
      n = n(),
      mean = mean(.data[[metric]], na.rm = TRUE),
      sd = sd(.data[[metric]], na.rm = TRUE),
      min = min(.data[[metric]], na.rm = TRUE),
      q50 = quantile(.data[[metric]], 0.50, na.rm = TRUE),
      q75 = quantile(.data[[metric]], 0.75, na.rm = TRUE),
      q90 = quantile(.data[[metric]], 0.90, na.rm = TRUE),
      q95 = quantile(.data[[metric]], 0.95, na.rm = TRUE),
      q99 = quantile(.data[[metric]], 0.99, na.rm = TRUE),
      max = max(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(metric = metric) %>%
    select(metric, source, everything())
}

# Compute tail coverage: % of each dataset above FIA's threshold
compute_tail_coverage <- function(df, metric, fia_quantile = 0.95) {
  # Get FIA threshold
  fia_threshold <- df %>%
    filter(source == "FIA") %>%
    pull(!!sym(metric)) %>%
    quantile(fia_quantile, na.rm = TRUE)
  
  # Compute % above threshold for each source
  df %>%
    group_by(source) %>%
    summarize(
      n_total = n(),
      n_above = sum(.data[[metric]] > fia_threshold, na.rm = TRUE),
      pct_above = 100 * n_above / n_total,
      threshold = fia_threshold,
      fia_quantile = fia_quantile,
      .groups = "drop"
    ) %>%
    mutate(metric = metric) %>%
    select(metric, everything())
}

# Create ECDF comparison plot
plot_ecdf <- function(df, metric, title, xlab, log_x = FALSE) {
  p <- ggplot(df, aes(x = .data[[metric]], color = source)) +
    stat_ecdf(linewidth = 1.2, alpha = 0.8) +
    scale_color_manual(values = colors, name = "Dataset") +
    labs(
      title = title,
      subtitle = "Empirical Cumulative Distribution",
      x = xlab,
      y = "Cumulative Proportion",
      caption = bias_caption
    ) +
    theme_edge()
  
  if (log_x) {
    p <- p + scale_x_log10()
  }
  
  # Add reference lines at key quantiles
  p <- p +
    geom_hline(yintercept = c(0.90, 0.95, 0.99), linetype = "dashed", 
               color = "gray70", alpha = 0.5)
  
  p
}

# Create density comparison plot
plot_density <- function(df, metric, title, xlab, log_x = FALSE) {
  p <- ggplot(df, aes(x = .data[[metric]], fill = source)) +
    geom_density(alpha = 0.5, color = NA) +
    geom_density(aes(color = source), fill = NA, linewidth = 0.8) +
    scale_fill_manual(values = colors, name = "Dataset") +
    scale_color_manual(values = colors, guide = "none") +
    labs(
      title = title,
      subtitle = "Density Distribution Comparison",
      x = xlab,
      y = "Density",
      caption = bias_caption
    ) +
    theme_edge()
  
  if (log_x) {
    p <- p + scale_x_log10()
  }
  
  p
}

# Create tail frequency bar chart
plot_tail_frequency <- function(tail_df, title) {
  ggplot(tail_df, aes(x = source, y = pct_above, fill = source)) +
    geom_col(alpha = 0.8, width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", pct_above)), vjust = -0.5, size = 4) +
    scale_fill_manual(values = colors, guide = "none") +
    labs(
      title = title,
      subtitle = sprintf("Plots above FIA 95th percentile (%.1f)", tail_df$threshold[1]),
      x = "",
      y = "% of Plots in Tail",
      caption = bias_caption
    ) +
    theme_edge() +
    coord_cartesian(ylim = c(0, max(tail_df$pct_above) * 1.2))
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

run_edge_case_comparison <- function(cfg = config) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
  cat("║                                                                              ║\n")
  cat("║   Edge Case Comparison: FIA vs NEFIN                                         ║\n")
  cat("║   Large-Tree Structure & Mortality Distributions                             ║\n")
  cat("║                                                                              ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  fs::dir_create(cfg$figures_dir, recurse = TRUE)
  fs::dir_create(cfg$tables_dir, recurse = TRUE)
  
  # Load data
  cat("Loading metrics data...\n")
  
  if (!fs::file_exists(cfg$metrics_path)) {
    stop("Metrics file not found: ", cfg$metrics_path, 
         "\n  Run 01_derive_edge_case_metrics.R first")
  }
  
  df <- read_csv(cfg$metrics_path, show_col_types = FALSE)
  cat("  Loaded", nrow(df), "plots\n")
  cat("  FIA:", sum(df$source == "FIA"), "| NEFIN:", sum(df$source == "NEFIN"), "\n")
  
  # ===========================================================================
  # 1. QUANTILE COMPARISON TABLES
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 1: Quantile Comparisons\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  metrics_to_compare <- c(
    "max_dbh", "max_tree_biomass", "p95_dbh", "p95_tree_biomass",
    "mortality_ratio", "dead_biomass_fraction", "pct_large_trees"
  )
  
  # Filter to metrics that exist
  metrics_to_compare <- intersect(metrics_to_compare, names(df))
  
  quantile_tables <- lapply(metrics_to_compare, function(m) {
    compute_quantile_comparison(df, m)
  })
  
  quantile_summary <- bind_rows(quantile_tables)
  
  cat("Quantile Summary:\n")
  print(quantile_summary %>% select(metric, source, q50, q90, q95, q99, max))
  
  write_csv(quantile_summary, fs::path(cfg$tables_dir, "quantile_comparison.csv"))
  cat("\n  ✓ Saved quantile_comparison.csv\n")
  
  # ===========================================================================
  # 2. TAIL COVERAGE ANALYSIS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 2: Tail Coverage (% Above FIA 95th Percentile)\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  tail_coverage <- lapply(metrics_to_compare, function(m) {
    compute_tail_coverage(df, m, fia_quantile = 0.95)
  })
  
  tail_summary <- bind_rows(tail_coverage)
  
  cat("Tail Coverage (% above FIA 95th percentile):\n")
  tail_wide <- tail_summary %>%
    select(metric, source, pct_above, threshold) %>%
    pivot_wider(names_from = source, values_from = pct_above, names_prefix = "pct_")
  print(tail_wide)
  
  write_csv(tail_summary, fs::path(cfg$tables_dir, "tail_coverage.csv"))
  cat("\n  ✓ Saved tail_coverage.csv\n")
  
  # ===========================================================================
  # 3. DISTRIBUTION FIGURES
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 3: Creating Distribution Figures\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  # --- Max DBH ---
  if ("max_dbh" %in% names(df)) {
    cat("  Creating max_dbh figures...\n")
    
    p1 <- plot_ecdf(df, "max_dbh", 
                    "Maximum Tree DBH per Plot",
                    "Max DBH (cm)")
    
    p2 <- plot_density(df, "max_dbh",
                       "Maximum Tree DBH Distribution",
                       "Max DBH (cm)")
    
    p3 <- plot_tail_frequency(
      tail_summary %>% filter(metric == "max_dbh"),
      "Plots with Large Trees (>FIA 95th)"
    )
    
    fig_dbh <- (p1 | p2) / p3 +
      plot_annotation(
        title = "Large-Tree Structure: Maximum DBH",
        theme = theme(plot.title = element_text(face = "bold", size = 16))
      )
    
    ggsave(fs::path(cfg$figures_dir, "ecdf_max_dbh.png"), 
           fig_dbh, width = 14, height = 10, dpi = 300)
    cat("    ✓ ecdf_max_dbh.png\n")
  }
  
  # --- Max Tree Biomass ---
  if ("max_tree_biomass" %in% names(df)) {
    cat("  Creating max_tree_biomass figures...\n")
    
    p1 <- plot_ecdf(df, "max_tree_biomass",
                    "Maximum Single-Tree Biomass per Plot",
                    "Max Tree Biomass (kg)", log_x = TRUE)
    
    p2 <- plot_density(df, "max_tree_biomass",
                       "Maximum Single-Tree Biomass Distribution",
                       "Max Tree Biomass (kg)", log_x = TRUE)
    
    fig_biomass <- p1 | p2
    
    ggsave(fs::path(cfg$figures_dir, "ecdf_max_tree_biomass.png"),
           fig_biomass, width = 14, height = 6, dpi = 300)
    cat("    ✓ ecdf_max_tree_biomass.png\n")
  }
  
  # --- Mortality ---
  if ("mortality_ratio" %in% names(df)) {
    cat("  Creating mortality figures...\n")
    
    p1 <- plot_ecdf(df, "mortality_ratio",
                    "Mortality Ratio (Dead/Total Trees)",
                    "Mortality Ratio")
    
    p2 <- plot_density(df %>% filter(mortality_ratio < 1), "mortality_ratio",
                       "Mortality Ratio Distribution",
                       "Mortality Ratio")
    
    fig_mort <- p1 | p2
    
    ggsave(fs::path(cfg$figures_dir, "ecdf_mortality_ratio.png"),
           fig_mort, width = 14, height = 6, dpi = 300)
    cat("    ✓ ecdf_mortality_ratio.png\n")
  }
  
  # --- Dead Biomass Fraction ---
  if ("dead_biomass_fraction" %in% names(df)) {
    cat("  Creating dead_biomass figures...\n")
    
    p1 <- plot_ecdf(df, "dead_biomass_fraction",
                    "Dead Biomass Fraction per Plot",
                    "Dead Biomass / Total Biomass")
    
    p2 <- plot_density(df %>% filter(dead_biomass_fraction < 1), 
                       "dead_biomass_fraction",
                       "Dead Biomass Fraction Distribution",
                       "Dead Biomass Fraction")
    
    fig_dead <- p1 | p2
    
    ggsave(fs::path(cfg$figures_dir, "ecdf_dead_biomass_fraction.png"),
           fig_dead, width = 14, height = 6, dpi = 300)
    cat("    ✓ ecdf_dead_biomass_fraction.png\n")
  }
  
  # ===========================================================================
  # 4. LARGE-TREE × MORTALITY INTERACTION
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 4: Large-Tree × Mortality Interaction\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  if (all(c("max_dbh", "dead_biomass_fraction") %in% names(df))) {
    cat("  Creating interaction scatter...\n")
    
    # Sample for visualization if too many points
    df_plot <- if (nrow(df) > 5000) df %>% sample_n(5000) else df
    
    p_scatter <- ggplot(df_plot, aes(x = max_dbh, y = dead_biomass_fraction, color = source)) +
      geom_point(alpha = 0.3, size = 1.5) +
      geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
      scale_color_manual(values = colors, name = "Dataset") +
      labs(
        title = "Large-Tree Structure vs Mortality",
        subtitle = "Does NEFIN capture more high-biomass plots with mortality?",
        x = "Maximum Tree DBH (cm)",
        y = "Dead Biomass Fraction",
        caption = bias_caption
      ) +
      theme_edge() +
      facet_wrap(~source)
    
    ggsave(fs::path(cfg$figures_dir, "scatter_dbh_vs_mortality.png"),
           p_scatter, width = 12, height = 6, dpi = 300)
    cat("    ✓ scatter_dbh_vs_mortality.png\n")
  }
  
  # ===========================================================================
  # 5. FOREST STATE TYPOLOGY
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 5: Forest State Typology\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  # Define thresholds from FIA
  if (all(c("max_dbh", "dead_biomass_fraction") %in% names(df))) {
    
    fia_data <- df %>% filter(source == "FIA")
    
    # Thresholds based on FIA
    large_dbh_thresh <- quantile(fia_data$max_dbh, 0.75, na.rm = TRUE)
    high_mort_thresh <- quantile(fia_data$dead_biomass_fraction, 0.75, na.rm = TRUE)
    
    cat("  Thresholds (FIA 75th percentile):\n")
    cat("    Large DBH:", round(large_dbh_thresh, 1), "cm\n")
    cat("    High Mortality:", round(high_mort_thresh, 3), "\n")
    
    # Classify plots
    df_states <- df %>%
      mutate(
        forest_state = case_when(
          max_dbh >= large_dbh_thresh & dead_biomass_fraction >= high_mort_thresh ~ 
            "Large-tree + High mortality",
          max_dbh >= large_dbh_thresh & dead_biomass_fraction < high_mort_thresh ~ 
            "Large-tree + Low mortality (intact)",
          max_dbh < large_dbh_thresh & dead_biomass_fraction >= high_mort_thresh ~ 
            "Small-tree + High mortality",
          TRUE ~ "Typical (middle)"
        )
      )
    
    state_summary <- df_states %>%
      group_by(source, forest_state) %>%
      summarize(n = n(), .groups = "drop") %>%
      group_by(source) %>%
      mutate(pct = 100 * n / sum(n)) %>%
      ungroup()
    
    cat("\n  Forest State Distribution:\n")
    print(state_summary %>% select(source, forest_state, n, pct) %>% 
            pivot_wider(names_from = source, values_from = c(n, pct)))
    
    write_csv(state_summary, fs::path(cfg$tables_dir, "forest_state_typology.csv"))
    cat("\n  ✓ Saved forest_state_typology.csv\n")
    
    # Plot
    p_states <- ggplot(state_summary, aes(x = forest_state, y = pct, fill = source)) +
      geom_col(position = "dodge", alpha = 0.8) +
      scale_fill_manual(values = colors, name = "Dataset") +
      labs(
        title = "Forest State Composition: FIA vs NEFIN",
        subtitle = "Based on FIA 75th percentile thresholds for max DBH and mortality",
        x = "",
        y = "% of Plots",
        caption = bias_caption
      ) +
      theme_edge() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    
    ggsave(fs::path(cfg$figures_dir, "forest_state_typology.png"),
           p_states, width = 10, height = 6, dpi = 300)
    cat("  ✓ forest_state_typology.png\n")
  }
  
  # ===========================================================================
  # 6. SUMMARY NOTES
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 6: Writing Summary Notes\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
  
  notes <- c(
    "# Edge Case Analysis: FIA vs NEFIN",
    "",
    "## Key Findings",
    "",
    "### 1. Large-Tree Structure",
    sprintf("- FIA 95th percentile max DBH: %.1f cm", 
            quantile(df$max_dbh[df$source == "FIA"], 0.95, na.rm = TRUE)),
    sprintf("- NEFIN plots above this threshold: %.1f%%",
            tail_summary %>% filter(metric == "max_dbh", source == "NEFIN") %>% pull(pct_above)),
    "",
    "### 2. Mortality/Disturbance",
    sprintf("- FIA median mortality ratio: %.3f",
            median(df$mortality_ratio[df$source == "FIA"], na.rm = TRUE)),
    sprintf("- NEFIN median mortality ratio: %.3f",
            median(df$mortality_ratio[df$source == "NEFIN"], na.rm = TRUE)),
    "",
    "## Interpretation Caveats",
    "",
    "**ALL findings must be interpreted with these caveats:**",
    "",
    "1. Differences reflect **sampling design**, not population differences",
    "2. NEFIN is biased toward public lands and research forests **by design**",
    "3. FIA is probability-based and suitable for population inference",
    "4. NEFIN may be better for modeling extremes, but not for unbiased means",
    "",
    "## Implications",
    "",
    "- For models needing extreme conditions: NEFIN may provide better training data",
    "- For population-level inference: FIA remains the gold standard",
    "- Combined use may improve models while maintaining interpretability",
    "",
    paste("Generated:", Sys.time())
  )
  
  writeLines(notes, fs::path(cfg$tables_dir, "edge_case_notes.md"))
  cat("  ✓ Saved edge_case_notes.md\n")
  
  # ===========================================================================
  # SUMMARY
  # ===========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  EDGE CASE ANALYSIS COMPLETE                                                 ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  cat("Outputs:\n")
  cat("  Tables:", cfg$tables_dir, "\n")
  cat("  Figures:", cfg$figures_dir, "\n")
  cat("\n")
  
  invisible(list(
    quantiles = quantile_summary,
    tail_coverage = tail_summary,
    data = df
  ))
}

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  run_edge_case_comparison()
}
