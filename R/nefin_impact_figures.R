#!/usr/bin/env Rscript
# R/nefin_impact_figures.R
# Publication-quality figures answering: Does NEFIN improve hex estimates?
# Author: Soren Donisvitch
# Date: December 2024

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sf)
  library(patchwork)
  library(viridis)
  library(scales)
  library(fs)
})

# Source utilities
if (file.exists("R/utils_scale_names.R")) source("R/utils_scale_names.R")

# =============================================================================
# MAIN FUNCTION
# =============================================================================

create_nefin_impact_figures <- function(
    consolidated_dir = NULL,
    output_dir = NULL
) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  NEFIN Impact Analysis: Does NEFIN Improve Hex Estimates?║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Find consolidated directory
  if (is.null(consolidated_dir)) {
    all_dirs <- list.dirs("runs", recursive = FALSE)
    consol_dirs <- all_dirs[grepl("consolidated_", all_dirs)]
    consol_dirs <- consol_dirs[order(file.mtime(consol_dirs), decreasing = TRUE)]
    consolidated_dir <- consol_dirs[1]
  }
  
  cat("Using:", consolidated_dir, "\n\n")
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- fs::path(consolidated_dir, "nefin_impact_figures")
  }
  fs::dir_create(output_dir, recurse = TRUE)
  
  # =========================================================================
  # Load Data
  # =========================================================================
  
  cat("Loading data...\n")
  
  # FIA results (all scales)
  fia <- read_csv(fs::path(consolidated_dir, "fia_all_scales.csv"), show_col_types = FALSE)
  
  # NEFIN comparison
  nefin <- read_csv(fs::path(consolidated_dir, "fia_nefin_comparison_all_scales.csv"), show_col_types = FALSE)
  
  # Standardize scale names
  if (exists("order_scales")) {
    fia <- order_scales(fia, "grid_scale")
    nefin <- order_scales(nefin, "grid_scale")
  } else {
    # Inline fallback
    scale_order <- c("100ha", "500ha", "1kha", "5kha", "10kha", "50kha", "64kha", "100kha")
    fia$grid_scale <- factor(fia$grid_scale, levels = scale_order, ordered = TRUE)
    nefin$grid_scale <- factor(nefin$grid_scale, levels = scale_order, ordered = TRUE)
  }
  
  cat("  FIA rows:", nrow(fia), "\n")
  cat("  NEFIN comparison rows:", nrow(nefin), "\n")
  
  # Filter to hexes with both FIA and NEFIN data
  nefin_both <- nefin %>% filter(has_both == TRUE)
  cat("  Hexes with both FIA & NEFIN:", nrow(nefin_both), "\n\n")
  
  # =========================================================================
  # FIGURE 1: Coverage - Where does NEFIN add data?
  # =========================================================================
  
  cat("Creating Figure 1: NEFIN Coverage by Scale...\n")
  
  coverage_summary <- nefin %>%
    group_by(grid_scale) %>%
    summarise(
      total_hexes = n(),
      fia_only = sum(has_fia & !has_nefin, na.rm = TRUE),
      nefin_only = sum(!has_fia & has_nefin, na.rm = TRUE),
      both = sum(has_both, na.rm = TRUE),
      pct_with_nefin = 100 * sum(has_nefin, na.rm = TRUE) / n(),
      pct_both = 100 * both / n(),
      .groups = "drop"
    )
  
  # Stacked bar of coverage
  coverage_long <- coverage_summary %>%
    select(grid_scale, fia_only, both, nefin_only) %>%
    pivot_longer(-grid_scale, names_to = "coverage", values_to = "count") %>%
    mutate(coverage = factor(coverage, 
                             levels = c("nefin_only", "both", "fia_only"),
                             labels = c("NEFIN Only", "Both", "FIA Only")))
  
  fig1a <- ggplot(coverage_long, aes(x = grid_scale, y = count, fill = coverage)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c("FIA Only" = "#3498db", "Both" = "#27ae60", "NEFIN Only" = "#f39c12")) +
    scale_y_continuous(labels = comma) +
    labs(
      title = "A) Data Coverage by Grid Scale",
      subtitle = "Number of hexagons with FIA, NEFIN, or both datasets",
      x = "Grid Scale", y = "Number of Hexagons",
      fill = "Data Source"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  # Percent with NEFIN overlap
  fig1b <- ggplot(coverage_summary, aes(x = grid_scale, y = pct_both)) +
    geom_col(fill = "#27ae60", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f%%", pct_both)), vjust = -0.5, size = 3.5) +
    scale_y_continuous(limits = c(0, max(coverage_summary$pct_both) * 1.15)) +
    labs(
      title = "B) NEFIN-FIA Overlap Rate",
      subtitle = "Percent of hexagons with data from both sources",
      x = "Grid Scale", y = "Percent with Both (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  fig1 <- fig1a + fig1b + 
    plot_annotation(
      title = "NEFIN Data Coverage Across Spatial Scales",
      theme = theme(plot.title = element_text(face = "bold", size = 16))
    )
  
  ggsave(fs::path(output_dir, "fig1_nefin_coverage.png"), fig1, 
         width = 14, height = 6, dpi = 300)
  cat("  ✓ Saved fig1_nefin_coverage.png\n")
  
  # =========================================================================
  # FIGURE 2: SE Reduction - Does NEFIN reduce uncertainty?
  # =========================================================================
  
  cat("Creating Figure 2: SE Reduction from NEFIN...\n")
  
  # Calculate SE reduction for hexes with both
  se_analysis <- nefin_both %>%
    filter(!is.na(se_fia), !is.na(se_nefin), se_fia > 0) %>%
    mutate(
      # Combined SE using inverse variance weighting approximation
      se_combined = sqrt(1 / (1/se_fia^2 + 1/se_nefin^2)),
      se_reduction_abs = se_fia - se_combined,
      se_reduction_pct = 100 * (se_fia - se_combined) / se_fia
    ) %>%
    filter(is.finite(se_reduction_pct))
  
  # Summary by scale
  se_summary <- se_analysis %>%
    group_by(grid_scale) %>%
    summarise(
      n_hexes = n(),
      mean_se_fia = mean(se_fia, na.rm = TRUE),
      mean_se_combined = mean(se_combined, na.rm = TRUE),
      mean_reduction_pct = mean(se_reduction_pct, na.rm = TRUE),
      median_reduction_pct = median(se_reduction_pct, na.rm = TRUE),
      sd_reduction_pct = sd(se_reduction_pct, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Bar chart of mean SE reduction
  fig2a <- ggplot(se_summary, aes(x = grid_scale, y = mean_reduction_pct, fill = grid_scale)) +
    geom_col(show.legend = FALSE) +
    geom_errorbar(aes(ymin = mean_reduction_pct - sd_reduction_pct/sqrt(n_hexes),
                      ymax = mean_reduction_pct + sd_reduction_pct/sqrt(n_hexes)),
                  width = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", mean_reduction_pct)), 
              vjust = -0.8, size = 3.5) +
    scale_fill_viridis_d(option = "D") +
    scale_y_continuous(limits = c(0, max(se_summary$mean_reduction_pct, na.rm = TRUE) * 1.2)) +
    labs(
      title = "A) Mean SE Reduction by Scale",
      subtitle = "Percent reduction in sampling error when adding NEFIN",
      x = "Grid Scale", y = "SE Reduction (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Before/After comparison
  se_before_after <- se_summary %>%
    select(grid_scale, mean_se_fia, mean_se_combined) %>%
    pivot_longer(-grid_scale, names_to = "source", values_to = "se") %>%
    mutate(source = factor(source, 
                           levels = c("mean_se_fia", "mean_se_combined"),
                           labels = c("FIA Only", "FIA + NEFIN")))
  
  fig2b <- ggplot(se_before_after, aes(x = grid_scale, y = se, fill = source)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("FIA Only" = "#e74c3c", "FIA + NEFIN" = "#27ae60")) +
    labs(
      title = "B) Standard Error: Before vs After NEFIN",
      subtitle = "Mean SE per hexagon",
      x = "Grid Scale", y = "Standard Error (Mg/ha)",
      fill = "Dataset"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  fig2 <- fig2a + fig2b +
    plot_annotation(
      title = "NEFIN Impact on Sampling Error",
      subtitle = "Does adding NEFIN plots reduce uncertainty in hexagonal estimates?",
      theme = theme(plot.title = element_text(face = "bold", size = 16))
    )
  
  ggsave(fs::path(output_dir, "fig2_se_reduction.png"), fig2,
         width = 14, height = 6, dpi = 300)
  cat("  ✓ Saved fig2_se_reduction.png\n")
  
  # =========================================================================
  # FIGURE 3: Agreement - Do FIA and NEFIN estimates agree?
  # =========================================================================
  
  cat("Creating Figure 3: FIA-NEFIN Agreement...\n")
  
  # Calculate agreement metrics
  agreement <- nefin_both %>%
    filter(!is.na(mean_fia), !is.na(mean_nefin)) %>%
    group_by(grid_scale) %>%
    summarise(
      n = n(),
      rmse = sqrt(mean((mean_fia - mean_nefin)^2, na.rm = TRUE)),
      mae = mean(abs(mean_fia - mean_nefin), na.rm = TRUE),
      bias = mean(mean_nefin - mean_fia, na.rm = TRUE),
      correlation = cor(mean_fia, mean_nefin, use = "complete.obs"),
      mean_biomass = mean((mean_fia + mean_nefin) / 2, na.rm = TRUE),
      rmse_pct = 100 * rmse / mean_biomass,
      .groups = "drop"
    )
  
  # Scatter plots by scale
  fig3a <- nefin_both %>%
    filter(!is.na(mean_fia), !is.na(mean_nefin)) %>%
    ggplot(aes(x = mean_fia, y = mean_nefin)) +
    geom_point(alpha = 0.4, size = 1, color = "#3498db") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "#2c3e50", linewidth = 0.8) +
    facet_wrap(~grid_scale, ncol = 4, scales = "free") +
    labs(
      title = "A) FIA vs NEFIN Biomass Estimates",
      subtitle = "Red dashed = 1:1 line | Dark line = linear fit",
      x = "FIA Estimate (Mg/ha)", y = "NEFIN Estimate (Mg/ha)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )
  
  # RMSE by scale
  fig3b <- ggplot(agreement, aes(x = grid_scale, y = rmse, fill = grid_scale)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.1f", rmse)), vjust = -0.5, size = 3.5) +
    scale_fill_viridis_d(option = "C") +
    labs(
      title = "B) RMSE Between FIA and NEFIN",
      x = "Grid Scale", y = "RMSE (Mg/ha)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Correlation by scale
  fig3c <- ggplot(agreement, aes(x = grid_scale, y = correlation, fill = grid_scale)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.2f", correlation)), vjust = -0.5, size = 3.5) +
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "darkgreen") +
    scale_fill_viridis_d(option = "D") +
    scale_y_continuous(limits = c(0, 1.1)) +
    labs(
      title = "C) Correlation (r)",
      subtitle = "Dashed line = r = 0.9",
      x = "Grid Scale", y = "Pearson r"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  fig3 <- fig3a / (fig3b | fig3c) +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(
      title = "FIA-NEFIN Agreement Across Scales",
      subtitle = "How well do fuzzed FIA and true NEFIN coordinates agree?",
      theme = theme(plot.title = element_text(face = "bold", size = 16))
    )
  
  ggsave(fs::path(output_dir, "fig3_agreement.png"), fig3,
         width = 14, height = 10, dpi = 300)
  cat("  ✓ Saved fig3_agreement.png\n")
  
  # =========================================================================
  # FIGURE 4: Scale Recommendation - Which scale benefits most?
  # =========================================================================
  
  cat("Creating Figure 4: Scale Recommendation...\n")
  
  # Combine all metrics for scale recommendation
  scale_metrics <- coverage_summary %>%
    select(grid_scale, total_hexes, pct_both) %>%
    left_join(se_summary %>% select(grid_scale, n_hexes, mean_reduction_pct), by = "grid_scale") %>%
    left_join(agreement %>% select(grid_scale, rmse, correlation), by = "grid_scale") %>%
    mutate(
      # Normalize metrics for comparison (0-1 scale)
      coverage_score = pct_both / max(pct_both, na.rm = TRUE),
      se_reduction_score = mean_reduction_pct / max(mean_reduction_pct, na.rm = TRUE),
      agreement_score = correlation,
      # Composite score (equal weights)
      composite_score = (coverage_score + se_reduction_score + agreement_score) / 3
    )
  
  # Heatmap of metrics
  metrics_long <- scale_metrics %>%
    select(grid_scale, coverage_score, se_reduction_score, agreement_score, composite_score) %>%
    pivot_longer(-grid_scale, names_to = "metric", values_to = "score") %>%
    mutate(metric = factor(metric,
                           levels = c("coverage_score", "se_reduction_score", "agreement_score", "composite_score"),
                           labels = c("NEFIN Coverage", "SE Reduction", "FIA-NEFIN Agreement", "COMPOSITE")))
  
  fig4a <- ggplot(metrics_long, aes(x = grid_scale, y = metric, fill = score)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.2f", score)), color = "white", fontface = "bold", size = 4) +
    scale_fill_viridis_c(option = "D", limits = c(0, 1)) +
    labs(
      title = "A) NEFIN Benefit Scores by Scale",
      subtitle = "Normalized 0-1 (higher = better)",
      x = "Grid Scale", y = "",
      fill = "Score"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  # Composite score bar chart
  fig4b <- ggplot(scale_metrics, aes(x = grid_scale, y = composite_score, fill = composite_score)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.2f", composite_score)), vjust = -0.5, size = 4, fontface = "bold") +
    scale_fill_viridis_c(option = "D") +
    scale_y_continuous(limits = c(0, 1.1)) +
    labs(
      title = "B) Composite NEFIN Benefit Score",
      subtitle = "Average of coverage, SE reduction, and agreement scores",
      x = "Grid Scale", y = "Composite Score"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Identify best scale
  best_scale <- scale_metrics %>% arrange(desc(composite_score)) %>% slice(1)
  
  fig4 <- fig4a + fig4b +
    plot_annotation(
      title = "Scale Recommendation: Where Does NEFIN Add Most Value?",
      subtitle = sprintf("Recommended scale: %s (composite score = %.2f)", 
                         best_scale$grid_scale, best_scale$composite_score),
      theme = theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12, color = "#27ae60", face = "bold")
      )
    )
  
  ggsave(fs::path(output_dir, "fig4_scale_recommendation.png"), fig4,
         width = 14, height = 7, dpi = 300)
  cat("  ✓ Saved fig4_scale_recommendation.png\n")
  
  # =========================================================================
  # FIGURE 5: Summary Dashboard
  # =========================================================================
  
  cat("Creating Figure 5: Summary Dashboard...\n")
  
  # Key finding text
  key_findings <- sprintf(
    "KEY FINDINGS:\n• NEFIN coverage highest at larger scales (%.0f%% at %s)\n• SE reduction greatest at smaller scales (%.1f%% at %s)\n• FIA-NEFIN agreement improves with scale (r = %.2f at %s)\n• RECOMMENDED SCALE: %s",
    max(coverage_summary$pct_both),
    coverage_summary$grid_scale[which.max(coverage_summary$pct_both)],
    max(se_summary$mean_reduction_pct, na.rm = TRUE),
    se_summary$grid_scale[which.max(se_summary$mean_reduction_pct)],
    max(agreement$correlation, na.rm = TRUE),
    agreement$grid_scale[which.max(agreement$correlation)],
    best_scale$grid_scale
  )
  
  dashboard <- (fig1a | fig2a) / (fig3b | fig3c | fig4b) +
    plot_annotation(
      title = "NEFIN Impact Summary: Does Adding NEFIN Improve Hexagonal Estimates?",
      subtitle = key_findings,
      theme = theme(
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(size = 11, family = "mono", hjust = 0)
      )
    )
  
  ggsave(fs::path(output_dir, "fig5_summary_dashboard.png"), dashboard,
         width = 16, height = 10, dpi = 300)
  cat("  ✓ Saved fig5_summary_dashboard.png\n")
  
  # =========================================================================
  # FIGURE 6: Spatial Maps - Where does NEFIN add value?
  # =========================================================================
  
  cat("Creating Figure 6: Spatial Maps of NEFIN Impact...\n")
  
  # Load hex grids and state boundaries
  hex_grids <- list()
  hex_dir <- "data/hex"
  
  # Load a few representative scales
  scales_to_map <- c("10kha", "50kha", "64kha", "100kha")
  
  for (scale in scales_to_map) {
    hex_file <- if (scale == "64kha") {
      fs::path(hex_dir, "hex_grid.geojson")
    } else {
      fs::path(hex_dir, paste0("hex_grid_", scale, ".geojson"))
    }
    
    if (fs::file_exists(hex_file)) {
      hex_grids[[scale]] <- st_read(hex_file, quiet = TRUE)
      # Ensure consistent ID column
      if (!"hex_id" %in% names(hex_grids[[scale]])) {
        id_col <- names(hex_grids[[scale]])[1]
        hex_grids[[scale]]$hex_id <- hex_grids[[scale]][[id_col]]
      }
    }
  }
  
  # Load state boundaries if available
  states_file <- "data/boundaries/states_5070.geojson"
  states <- if (fs::file_exists(states_file)) {
    st_read(states_file, quiet = TRUE)
  } else {
    NULL
  }
  
  # Create maps for each scale
  map_list <- list()
  
  for (scale in names(hex_grids)) {
    hex <- hex_grids[[scale]]
    
    # Get NEFIN data for this scale
    nefin_scale <- nefin %>% 
      filter(grid_scale == scale) %>%
      mutate(
        data_status = case_when(
          has_both ~ "Both FIA & NEFIN",
          has_fia & !has_nefin ~ "FIA Only",
          !has_fia & has_nefin ~ "NEFIN Only",
          TRUE ~ "No Data"
        )
      )
    
    # Join to hex grid
    hex_joined <- hex %>%
      left_join(nefin_scale, by = "hex_id") %>%
      mutate(data_status = ifelse(is.na(data_status), "No Data", data_status))
    
    # Create map
    p <- ggplot() +
      geom_sf(data = hex_joined, aes(fill = data_status), color = "gray40", linewidth = 0.1) +
      scale_fill_manual(
        values = c(
          "Both FIA & NEFIN" = "#27ae60",
          "FIA Only" = "#3498db", 
          "NEFIN Only" = "#f39c12",
          "No Data" = "gray90"
        ),
        na.value = "gray90"
      ) +
      labs(title = paste0(scale, " Grid"), fill = "Data Source") +
      theme_void(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none"
      )
    
    # Add state boundaries if available
    if (!is.null(states)) {
      p <- p + geom_sf(data = states, fill = NA, color = "black", linewidth = 0.3)
    }
    
    map_list[[scale]] <- p
  }
  
  # Combine maps with shared legend
  if (length(map_list) >= 4) {
    # Create a legend
    legend_data <- data.frame(
      status = factor(c("Both FIA & NEFIN", "FIA Only", "NEFIN Only", "No Data"),
                      levels = c("Both FIA & NEFIN", "FIA Only", "NEFIN Only", "No Data"))
    )
    
    legend_plot <- ggplot(legend_data, aes(x = 1, y = status, fill = status)) +
      geom_tile() +
      scale_fill_manual(
        values = c(
          "Both FIA & NEFIN" = "#27ae60",
          "FIA Only" = "#3498db", 
          "NEFIN Only" = "#f39c12",
          "No Data" = "gray90"
        )
      ) +
      labs(fill = "Data Source") +
      theme_void() +
      theme(legend.position = "bottom")
    
    fig6_coverage <- (map_list[[1]] | map_list[[2]]) / (map_list[[3]] | map_list[[4]]) +
      plot_annotation(
        title = "Spatial Distribution of NEFIN Coverage by Scale",
        subtitle = "Green = hexagons with both FIA and NEFIN data",
        theme = theme(
          plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 12)
        )
      )
    
    ggsave(fs::path(output_dir, "fig6_spatial_coverage.png"), fig6_coverage,
           width = 12, height = 10, dpi = 300)
    cat("  ✓ Saved fig6_spatial_coverage.png\n")
  }
  
  # =========================================================================
  # FIGURE 7: SE Reduction Map - Where does NEFIN help most?
  # =========================================================================
  
  cat("Creating Figure 7: SE Reduction Maps...\n")
  
  # Use a mid-range scale for SE reduction visualization
  map_scale <- "50kha"  
  if (!map_scale %in% names(hex_grids)) {
    map_scale <- names(hex_grids)[1]
  }
  
  hex <- hex_grids[[map_scale]]
  
  # Get SE data
  se_map_data <- se_analysis %>%
    filter(grid_scale == map_scale) %>%
    select(hex_id, se_fia, se_combined, se_reduction_pct)
  
  if (nrow(se_map_data) > 0) {
    hex_se <- hex %>%
      left_join(se_map_data, by = "hex_id")
    
    # Map of SE reduction
    fig7a <- ggplot() +
      geom_sf(data = hex_se, aes(fill = se_reduction_pct), color = "gray50", linewidth = 0.1) +
      scale_fill_viridis_c(
        option = "D", 
        na.value = "gray90",
        limits = c(0, 50),
        oob = scales::squish,
        labels = function(x) paste0(x, "%")
      ) +
      labs(
        title = paste0("A) SE Reduction (", map_scale, ")"),
        subtitle = "Percent reduction in standard error when adding NEFIN",
        fill = "SE Reduction"
      ) +
      theme_void(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm")
      )
    
    if (!is.null(states)) {
      fig7a <- fig7a + geom_sf(data = states, fill = NA, color = "black", linewidth = 0.4)
    }
    
    # Map of original SE (FIA only)
    fig7b <- ggplot() +
      geom_sf(data = hex_se, aes(fill = se_fia), color = "gray50", linewidth = 0.1) +
      scale_fill_viridis_c(
        option = "C", 
        na.value = "gray90",
        limits = c(0, 100),
        oob = scales::squish
      ) +
      labs(
        title = "B) Original SE (FIA Only)",
        subtitle = "Standard error before adding NEFIN",
        fill = "SE (Mg/ha)"
      ) +
      theme_void(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm")
      )
    
    if (!is.null(states)) {
      fig7b <- fig7b + geom_sf(data = states, fill = NA, color = "black", linewidth = 0.4)
    }
    
    fig7 <- fig7a + fig7b +
      plot_annotation(
        title = "Where Does NEFIN Reduce Uncertainty?",
        theme = theme(plot.title = element_text(face = "bold", size = 16))
      )
    
    ggsave(fs::path(output_dir, "fig7_se_reduction_map.png"), fig7,
           width = 14, height = 7, dpi = 300)
    cat("  ✓ Saved fig7_se_reduction_map.png\n")
  } else {
    cat("  ⚠ No SE data available for mapping\n")
  }
  
  # =========================================================================
  # FIGURE 8: FIA-NEFIN Difference Map
  # =========================================================================
  
  cat("Creating Figure 8: FIA-NEFIN Difference Map...\n")
  
  # Get biomass difference data
  diff_map_data <- nefin_both %>%
    filter(grid_scale == map_scale, !is.na(mean_fia), !is.na(mean_nefin)) %>%
    mutate(
      biomass_diff = mean_nefin - mean_fia,
      pct_diff = 100 * biomass_diff / ((mean_fia + mean_nefin) / 2)
    ) %>%
    select(hex_id, mean_fia, mean_nefin, biomass_diff, pct_diff)
  
  if (nrow(diff_map_data) > 0) {
    hex_diff <- hex %>%
      left_join(diff_map_data, by = "hex_id")
    
    # Percent difference map
    fig8a <- ggplot() +
      geom_sf(data = hex_diff, aes(fill = pct_diff), color = "gray50", linewidth = 0.1) +
      scale_fill_gradient2(
        low = "#e74c3c", mid = "white", high = "#3498db",
        midpoint = 0,
        na.value = "gray90",
        limits = c(-50, 50),
        oob = scales::squish,
        labels = function(x) paste0(x, "%")
      ) +
      labs(
        title = paste0("A) NEFIN vs FIA Difference (", map_scale, ")"),
        subtitle = "Blue = NEFIN higher, Red = FIA higher",
        fill = "% Difference"
      ) +
      theme_void(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm")
      )
    
    if (!is.null(states)) {
      fig8a <- fig8a + geom_sf(data = states, fill = NA, color = "black", linewidth = 0.4)
    }
    
    # Mean biomass map (average of FIA and NEFIN)
    hex_diff <- hex_diff %>%
      mutate(mean_biomass = (mean_fia + mean_nefin) / 2)
    
    fig8b <- ggplot() +
      geom_sf(data = hex_diff, aes(fill = mean_biomass), color = "gray50", linewidth = 0.1) +
      scale_fill_viridis_c(
        option = "G", 
        na.value = "gray90",
        direction = -1
      ) +
      labs(
        title = "B) Mean Biomass (FIA+NEFIN Average)",
        subtitle = "Hexagons with both data sources",
        fill = "Biomass\n(Mg/ha)"
      ) +
      theme_void(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm")
      )
    
    if (!is.null(states)) {
      fig8b <- fig8b + geom_sf(data = states, fill = NA, color = "black", linewidth = 0.4)
    }
    
    fig8 <- fig8a + fig8b +
      plot_annotation(
        title = "FIA vs NEFIN Biomass Comparison",
        subtitle = "Spatial patterns in estimate differences",
        theme = theme(
          plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 12)
        )
      )
    
    ggsave(fs::path(output_dir, "fig8_difference_map.png"), fig8,
           width = 14, height = 7, dpi = 300)
    cat("  ✓ Saved fig8_difference_map.png\n")
  } else {
    cat("  ⚠ No difference data available for mapping\n")
  }
  
  # =========================================================================
  # Save Summary Table
  # =========================================================================
  
  cat("\nSaving summary tables...\n")
  
  # Build summary table more carefully
  summary_table <- scale_metrics %>%
    select(grid_scale, composite_score, coverage_score, se_reduction_score, agreement_score)
  
  # Add coverage info if available
  if ("pct_both" %in% names(coverage_summary)) {
    summary_table <- summary_table %>%
      left_join(coverage_summary %>% select(grid_scale, any_of(c("total_hexes", "both", "pct_both"))), by = "grid_scale")
  }
  
  # Add SE info if available
  if (nrow(se_summary) > 0) {
    summary_table <- summary_table %>%
      left_join(se_summary %>% select(grid_scale, any_of(c("n_hexes", "mean_se_fia", "mean_se_combined", "mean_reduction_pct"))), by = "grid_scale")
  }
  
  # Add agreement info if available
  if (nrow(agreement) > 0) {
    summary_table <- summary_table %>%
      left_join(agreement %>% select(grid_scale, any_of(c("rmse", "bias", "correlation"))), by = "grid_scale")
  }
  
  summary_table <- summary_table %>% arrange(desc(composite_score))
  
  write_csv(summary_table, fs::path(output_dir, "nefin_impact_summary.csv"))
  cat("  ✓ Saved nefin_impact_summary.csv\n")
  
  # =========================================================================
  # Print Summary
  # =========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  ANALYSIS COMPLETE                                        ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  cat(key_findings, "\n\n")
  
  cat("Output files in:", output_dir, "\n")
  cat("  • fig1_nefin_coverage.png      - Coverage by scale (bar charts)\n")
  cat("  • fig2_se_reduction.png        - SE reduction by scale\n")
  cat("  • fig3_agreement.png           - FIA-NEFIN agreement (scatter + metrics)\n")
  cat("  • fig4_scale_recommendation.png - Composite scores by scale\n")
  cat("  • fig5_summary_dashboard.png   - All key findings combined\n")
  cat("  • fig6_spatial_coverage.png    - Maps of NEFIN coverage\n")
  cat("  • fig7_se_reduction_map.png    - Maps of where NEFIN reduces SE\n")
  cat("  • fig8_difference_map.png      - Maps of FIA-NEFIN differences\n")
  cat("  • nefin_impact_summary.csv     - Summary statistics table\n")
  
  invisible(summary_table)
}

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  create_nefin_impact_figures()
}
