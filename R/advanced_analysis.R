# R/advanced_analysis.R
# Advanced statistical analysis of FIA error and NEFIN comparison

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(tidyr)
  library(broom); library(scales); library(viridis); library(patchwork); library(fs)
})

advanced_analysis <- function(consolidated_dir = NULL,
                              output_dir = NULL) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Advanced Statistical Analysis                           ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Find most recent consolidated directory if not specified
  if (is.null(consolidated_dir)) {
    all_dirs <- list.dirs("runs", recursive = FALSE)
    consol_dirs <- all_dirs[grepl("^runs/consolidated_", all_dirs)]
    
    if (length(consol_dirs) == 0) {
      stop("No consolidated results found. Run R/master_process_all.R first.")
    }
    
    consol_dirs <- consol_dirs[order(file.mtime(consol_dirs), decreasing = TRUE)]
    consolidated_dir <- consol_dirs[1]
  }
  
  cat("Using consolidated results from:\n")
  cat("  ", consolidated_dir, "\n\n")
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- fs::path(consolidated_dir, "advanced_analysis")
  }
  fs::dir_create(output_dir, recurse = TRUE)
  
  # Load data
  cat("Loading data...\n")
  fia_file <- fs::path(consolidated_dir, "fia_all_scales.csv")
  nefin_file <- fs::path(consolidated_dir, "fia_nefin_comparison_all_scales.csv")
  
  fia <- readr::read_csv(fia_file, show_col_types = FALSE)
  cat("  ✓ FIA results loaded\n")
  
  has_nefin <- fs::file_exists(nefin_file)
  nefin <- if (has_nefin) {
    readr::read_csv(nefin_file, show_col_types = FALSE)
  } else NULL
  
  if (has_nefin) cat("  ✓ NEFIN comparison loaded\n")
  
  # ========================================================================
  # ANALYSIS 1: Does NEFIN reduce sampling error?
  # ========================================================================
  
  if (!is.null(nefin)) {
    cat("\n")
    cat("══════════════════════════════════════════════════════════\n")
    cat("ANALYSIS 1: NEFIN Impact on Sampling Error\n")
    cat("══════════════════════════════════════════════════════════\n")
    
    nefin_impact <- nefin |>
      dplyr::filter(has_both, !is.na(se_fia), !is.na(n_plots_fia), !is.na(n_plots_nefin)) |>
      dplyr::mutate(
        n_combined = n_plots_fia + n_plots_nefin,
        
        # Expected SE reduction from larger sample size
        # SE is proportional to 1/sqrt(n)
        se_combined_expected = se_fia * sqrt(n_plots_fia / n_combined),
        
        # Percent reduction in SE
        se_reduction_pct = 100 * (1 - se_combined_expected / se_fia),
        
        # Absolute reduction
        se_reduction_abs = se_fia - se_combined_expected
      )
    
    # Summary by scale
    impact_summary <- nefin_impact |>
      dplyr::group_by(grid_scale) |>
      dplyr::summarise(
        n_hexes = dplyr::n(),
        mean_fia_plots = mean(n_plots_fia),
        mean_nefin_plots = mean(n_plots_nefin),
        mean_se_fia = mean(se_fia),
        mean_se_combined = mean(se_combined_expected),
        mean_reduction_pct = mean(se_reduction_pct),
        median_reduction_pct = median(se_reduction_pct),
        mean_reduction_abs = mean(se_reduction_abs),
        .groups = "drop"
      )
    
    readr::write_csv(impact_summary, fs::path(output_dir, "nefin_se_reduction_summary.csv"))
    
    cat("\n→ SE Reduction Summary by Scale:\n")
    print(impact_summary, n = 100)
    
    # Plot 1: SE reduction distribution
    p1 <- ggplot(nefin_impact, aes(x = grid_scale, y = se_reduction_pct, fill = grid_scale)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.2, fill = "white", outlier.alpha = 0.3) +
      scale_fill_viridis_d(option = "C") +
      labs(title = "Expected Sampling Error Reduction from Adding NEFIN",
           subtitle = "Based on combined sample size (assumes independent samples)",
           x = "Grid Scale", y = "SE Reduction (%)") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 14),
            legend.position = "none")
    
    ggsave(fs::path(output_dir, "01_nefin_se_reduction.png"), p1,
           width = 10, height = 6, dpi = 300)
    
    # Plot 2: Relationship between NEFIN plots and SE reduction
    p2 <- ggplot(nefin_impact, aes(x = n_plots_nefin, y = se_reduction_pct, color = grid_scale)) +
      geom_point(alpha = 0.4) +
      geom_smooth(method = "loess", se = TRUE) +
      scale_color_viridis_d(option = "D") +
      scale_x_log10() +
      labs(title = "SE Reduction vs Number of NEFIN Plots",
           x = "NEFIN Plots (log scale)", y = "SE Reduction (%)", color = "Grid Scale") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 14),
            legend.position = "bottom")
    
    ggsave(fs::path(output_dir, "02_se_reduction_vs_nefin_plots.png"), p2,
           width = 10, height = 6, dpi = 300)
  }
  
  # ========================================================================
  # ANALYSIS 2: Scale effects on positional uncertainty
  # ========================================================================
  
  cat("\n")
  cat("══════════════════════════════════════════════════════════\n")
  cat("ANALYSIS 2: Grid Scale Effects on Positional Error\n")
  cat("══════════════════════════════════════════════════════════\n")
  
  # Map grid names to approximate areas
  scale_effects <- fia |>
    dplyr::mutate(
      grid_area_acres = dplyr::case_when(
        grepl("fia", grid_scale, ignore.case = TRUE) ~ 6000,
        grepl("1.5", grid_scale) ~ 1500,
        grepl("3", grid_scale) ~ 3000,
        grepl("6", grid_scale) ~ 6000,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::filter(!is.na(grid_area_acres))
  
  # Regression: positional SD ~ log(grid area)
  scale_model <- scale_effects |>
    dplyr::filter(is.finite(positional_sd), positional_sd > 0, grid_area_acres > 0) |>
    dplyr::mutate(log_area = log10(grid_area_acres))
  
  if (nrow(scale_model) > 0) {
    lm_fit <- lm(positional_sd ~ log_area, data = scale_model)
    lm_summary <- broom::tidy(lm_fit)
    
    readr::write_csv(lm_summary, fs::path(output_dir, "scale_effect_model.csv"))
    
    cat("\n→ Linear Model: Positional SD ~ log10(Grid Area)\n")
    print(lm_summary)
    cat("\n  R²: ", round(summary(lm_fit)$r.squared, 4), "\n")
    
    # Summary by scale
    scale_summary <- scale_effects |>
      dplyr::group_by(grid_scale, grid_area_acres) |>
      dplyr::summarise(
        n = dplyr::n(),
        mean_pos_sd = mean(positional_sd, na.rm = TRUE),
        median_pos_sd = median(positional_sd, na.rm = TRUE),
        sd_pos_sd = sd(positional_sd, na.rm = TRUE),
        mean_se = mean(se, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(grid_area_acres)
    
    readr::write_csv(scale_summary, fs::path(output_dir, "positional_sd_by_scale.csv"))
    
    cat("\n→ Positional SD by Grid Scale:\n")
    print(scale_summary, n = 100)
    
    # Plot 3: Scale effect
    p3 <- ggplot(scale_summary, aes(x = grid_area_acres, y = mean_pos_sd)) +
      geom_line(linewidth = 1.2, color = "steelblue") +
      geom_point(size = 4, color = "steelblue") +
      geom_errorbar(aes(ymin = mean_pos_sd - sd_pos_sd/sqrt(n),
                        ymax = mean_pos_sd + sd_pos_sd/sqrt(n)),
                    width = 0.1, linewidth = 1) +
      geom_text(aes(label = grid_scale), vjust = -1, size = 4) +
      scale_x_log10(breaks = c(1500, 3000, 6000), labels = comma) +
      labs(title = "Positional Uncertainty vs Grid Size",
           subtitle = "Error bars show ±1 SE",
           x = "Grid Area (acres, log scale)", y = "Mean Positional SD (Mg/ha)") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 14))
    
    ggsave(fs::path(output_dir, "03_scale_effect.png"), p3,
           width = 10, height = 6, dpi = 300)
  }
  
  # ========================================================================
  # ANALYSIS 3: Where does NEFIN provide most value?
  # ========================================================================
  
  if (!is.null(nefin)) {
    cat("\n")
    cat("══════════════════════════════════════════════════════════\n")
    cat("ANALYSIS 3: Where NEFIN Adds Most Value\n")
    cat("══════════════════════════════════════════════════════════\n")
    
    # Identify hexes where NEFIN supplements sparse FIA data
    nefin_value <- nefin |>
      dplyr::filter(has_both) |>
      dplyr::mutate(
        fia_sparse = n_plots_fia < 5,
        nefin_substantial = n_plots_nefin >= 3,
        adds_value = fia_sparse & nefin_substantial,
        
        # Potential SE improvement
        se_improvement = ifelse(
          !is.na(se_fia) & !is.na(n_plots_fia) & !is.na(n_plots_nefin),
          se_fia * (1 - sqrt(n_plots_fia / (n_plots_fia + n_plots_nefin))),
          NA_real_
        )
      )
    
    value_hexes <- nefin_value |>
      dplyr::filter(adds_value) |>
      dplyr::arrange(desc(se_improvement)) |>
      dplyr::select(hex_id, grid_scale, year_label, n_plots_fia, n_plots_nefin,
                    mean_fia, mean_nefin, diff, se_fia, se_improvement) |>
      dplyr::slice_head(n = 50)
    
    readr::write_csv(value_hexes, fs::path(output_dir, "high_value_nefin_hexes.csv"))
    
    cat("\n→ High-Value NEFIN Hexes (FIA sparse, NEFIN supplements):\n")
    cat("  Total hexes where NEFIN adds value: ", sum(nefin_value$adds_value), "\n")
    cat("  Top 20 hexes by potential SE improvement:\n")
    print(as.data.frame(head(value_hexes, 20)), row.names = FALSE)
    
    # Plot 4: NEFIN value map
    value_summary <- nefin_value |>
      dplyr::group_by(grid_scale) |>
      dplyr::summarise(
        total_hexes = dplyr::n(),
        sparse_fia = sum(fia_sparse),
        nefin_supplements = sum(adds_value),
        pct_supplemented = 100 * nefin_supplements / total_hexes,
        .groups = "drop"
      )
    
    p4 <- ggplot(value_summary, aes(x = grid_scale, y = pct_supplemented, fill = grid_scale)) +
      geom_col() +
      geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", pct_supplemented, 
                                    nefin_supplements, total_hexes)),
                vjust = -0.5, size = 4) +
      scale_fill_viridis_d(option = "C") +
      labs(title = "Where NEFIN Supplements Sparse FIA Data",
           subtitle = "Hexes with <5 FIA plots but ≥3 NEFIN plots",
           x = "Grid Scale", y = "% of Hexes Supplemented") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 14),
            legend.position = "none")
    
    ggsave(fs::path(output_dir, "04_nefin_value_coverage.png"), p4,
           width = 10, height = 6, dpi = 300)
  }
  
  # ========================================================================
  # ANALYSIS 4: Temporal patterns (if multiple years)
  # ========================================================================
  
  years <- unique(fia$year_label)
  
  if (length(years) > 1) {
    cat("\n")
    cat("══════════════════════════════════════════════════════════\n")
    cat("ANALYSIS 4: Temporal Trends\n")
    cat("══════════════════════════════════════════════════════════\n")
    
    # FIA temporal trends
    fia_temporal <- fia |>
      dplyr::group_by(year_label, grid_scale) |>
      dplyr::summarise(
        mean_biomass = mean(mean, na.rm = TRUE),
        mean_se = mean(se, na.rm = TRUE),
        mean_pos_sd = mean(positional_sd, na.rm = TRUE),
        mean_total_sd = mean(total_sd, na.rm = TRUE),
        n_hexes = dplyr::n(),
        .groups = "drop"
      )
    
    readr::write_csv(fia_temporal, fs::path(output_dir, "fia_temporal_trends.csv"))
    
    # Plot 5: Temporal error trends
    p5 <- fia_temporal |>
      tidyr::pivot_longer(cols = c(mean_se, mean_pos_sd, mean_total_sd),
                          names_to = "error_type", values_to = "value") |>
      dplyr::mutate(
        error_type = factor(error_type,
                            levels = c("mean_se", "mean_pos_sd", "mean_total_sd"),
                            labels = c("Sampling Error", "Positional SD", "Total SD"))
      ) |>
      ggplot(aes(x = year_label, y = value, color = error_type, group = error_type)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      facet_wrap(~grid_scale, ncol = 2) +
      scale_color_viridis_d(option = "D", begin = 0.3, end = 0.9) +
      labs(title = "Error Components Over Time by Grid Scale",
           x = "Year", y = "Standard Deviation (Mg/ha)", color = "Error Type") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 14),
            legend.position = "bottom")
    
    ggsave(fs::path(output_dir, "05_temporal_error_trends.png"), p5,
           width = 12, height = 8, dpi = 300)
    
    # NEFIN temporal agreement
    if (!is.null(nefin)) {
      nefin_temporal <- nefin |>
        dplyr::filter(has_both) |>
        dplyr::group_by(year_label, grid_scale) |>
        dplyr::summarise(
          n_hexes = dplyr::n(),
          rmse = sqrt(mean(diff^2, na.rm = TRUE)),
          correlation = cor(mean_fia, mean_nefin, use = "complete.obs"),
          mean_bias = mean(diff, na.rm = TRUE),
          .groups = "drop"
        )
      
      readr::write_csv(nefin_temporal, fs::path(output_dir, "nefin_temporal_agreement.csv"))
      
      # Plot 6: Temporal agreement
      p6 <- ggplot(nefin_temporal, aes(x = year_label, y = correlation, 
                                       color = grid_scale, group = grid_scale)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 3) +
        scale_color_viridis_d(option = "D") +
        labs(title = "FIA-NEFIN Agreement Over Time",
             x = "Year", y = "Correlation", color = "Grid Scale") +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(face = "bold", size = 14),
              legend.position = "bottom")
      
      ggsave(fs::path(output_dir, "06_temporal_agreement.png"), p6,
             width = 10, height = 6, dpi = 300)
    }
  }
  
  # ========================================================================
  # ANALYSIS 5: Top disagreement hexes
  # ========================================================================
  
  if (!is.null(nefin)) {
    cat("\n")
    cat("══════════════════════════════════════════════════════════\n")
    cat("ANALYSIS 5: Largest FIA-NEFIN Disagreements\n")
    cat("══════════════════════════════════════════════════════════\n")
    
    top_disagreements <- nefin |>
      dplyr::filter(has_both) |>
      dplyr::mutate(abs_diff = abs(diff)) |>
      dplyr::group_by(grid_scale) |>
      dplyr::arrange(desc(abs_diff)) |>
      dplyr::slice_head(n = 20) |>
      dplyr::ungroup() |>
      dplyr::select(hex_id, grid_scale, year_label, mean_fia, mean_nefin, diff, abs_diff,
                    n_plots_fia, n_plots_nefin, se_fia, se_nefin)
    
    readr::write_csv(top_disagreements, fs::path(output_dir, "top_disagreement_hexes.csv"))
    
    cat("\n→ Top 10 disagreements by grid scale:\n")
    print(as.data.frame(head(top_disagreements, 10)), row.names = FALSE)
    
    # Plot 7: Top disagreements
    p7 <- top_disagreements |>
      dplyr::group_by(grid_scale) |>
      dplyr::slice_head(n = 10) |>
      dplyr::mutate(hex_id_short = substr(hex_id, 1, 8)) |>
      ggplot(aes(x = reorder(hex_id_short, abs_diff), y = abs_diff, fill = grid_scale)) +
      geom_col() +
      coord_flip() +
      facet_wrap(~grid_scale, scales = "free_y") +
      scale_fill_viridis_d(option = "C") +
      labs(title = "Top 10 Hexes with Largest FIA-NEFIN Disagreement",
           subtitle = "By grid scale",
           x = "Hex ID", y = "Absolute Difference (Mg/ha)") +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold", size = 14),
            legend.position = "none")
    
    ggsave(fs::path(output_dir, "07_top_disagreements.png"), p7,
           width = 12, height = 10, dpi = 300)
  }
  
  # ========================================================================
  # ANALYSIS 6: Error decomposition statistics
  # ========================================================================
  
  cat("\n")
  cat("══════════════════════════════════════════════════════════\n")
  cat("ANALYSIS 6: Detailed Error Decomposition\n")
  cat("══════════════════════════════════════════════════════════\n")
  
  error_decomp <- fia |>
    dplyr::mutate(
      pos_fraction = positional_sd / total_sd,
      sampling_fraction = se / total_sd,
      
      # Classify dominant error source
      dominant_error = dplyr::case_when(
        pos_fraction > 0.6 ~ "Positional",
        sampling_fraction > 0.6 ~ "Sampling",
        TRUE ~ "Mixed"
      )
    ) |>
    dplyr::group_by(grid_scale, dominant_error) |>
    dplyr::summarise(
      n_hexes = dplyr::n(),
      mean_pos_sd = mean(positional_sd, na.rm = TRUE),
      mean_se = mean(se, na.rm = TRUE),
      mean_total_sd = mean(total_sd, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::group_by(grid_scale) |>
    dplyr::mutate(pct_hexes = 100 * n_hexes / sum(n_hexes))
  
  readr::write_csv(error_decomp, fs::path(output_dir, "error_decomposition_detailed.csv"))
  
  cat("\n→ Error Source Classification:\n")
  print(error_decomp, n = 100)
  
  # Plot 8: Dominant error source
  p8 <- ggplot(error_decomp, aes(x = grid_scale, y = pct_hexes, fill = dominant_error)) +
    geom_col(position = "stack") +
    geom_text(aes(label = sprintf("%.1f%%", pct_hexes)),
              position = position_stack(vjust = 0.5),
              color = "white", fontface = "bold", size = 3.5) +
    scale_fill_manual(values = c("Positional" = "#E63946",
                                 "Sampling" = "#457B9D",
                                 "Mixed" = "#F1FAEE")) +
    labs(title = "Dominant Error Source by Grid Scale",
         subtitle = "Positional >60% of total | Sampling >60% | Mixed otherwise",
         x = "Grid Scale", y = "% of Hexes", fill = "Dominant Error") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom")
  
  ggsave(fs::path(output_dir, "08_dominant_error_source.png"), p8,
         width = 10, height = 6, dpi = 300)
  
  # ========================================================================
  # SUMMARY REPORT
  # ========================================================================
  
  cat("\n")
  cat("══════════════════════════════════════════════════════════\n")
  cat("Creating Summary Report\n")
  cat("══════════════════════════════════════════════════════════\n")
  
  report_lines <- c(
    "═══════════════════════════════════════════════════════════════",
    "ADVANCED STATISTICAL ANALYSIS REPORT",
    "═══════════════════════════════════════════════════════════════",
    "",
    paste("Generated:", Sys.time()),
    paste("Data source:", consolidated_dir),
    ""
  )
  
  if (!is.null(nefin) && exists("impact_summary")) {
    report_lines <- c(report_lines,
                      "───────────────────────────────────────────────────────────────",
                      "1. NEFIN IMPACT ON SAMPLING ERROR",
                      "───────────────────────────────────────────────────────────────",
                      "",
                      "Expected SE reduction from combining FIA + NEFIN:",
                      ""
    )
    
    for (i in seq_len(nrow(impact_summary))) {
      row <- impact_summary[i,]
      report_lines <- c(report_lines,
                        paste0("Grid: ", row$grid_scale),
                        paste0("  Mean SE reduction: ", round(row$mean_reduction_pct, 1), "%"),
                        paste0("  Absolute reduction: ", round(row$mean_reduction_abs, 3), " Mg/ha"),
                        paste0("  FIA plots/hex: ", round(row$mean_fia_plots, 1)),
                        paste0("  NEFIN plots/hex: ", round(row$mean_nefin_plots, 1)),
                        ""
      )
    }
  }
  
  if (exists("scale_summary")) {
    report_lines <- c(report_lines,
                      "───────────────────────────────────────────────────────────────",
                      "2. GRID SCALE EFFECTS ON POSITIONAL ERROR",
                      "───────────────────────────────────────────────────────────────",
                      ""
    )
    
    for (i in seq_len(nrow(scale_summary))) {
      row <- scale_summary[i,]
      report_lines <- c(report_lines,
                        paste0("Grid: ", row$grid_scale, " (", format(row$grid_area_acres, big.mark=","), " acres)"),
                        paste0("  Mean positional SD: ", round(row$mean_pos_sd, 3), " Mg/ha"),
                        paste0("  Mean sampling SE: ", round(row$mean_se, 3), " Mg/ha"),
                        ""
      )
    }
    
    if (exists("lm_summary")) {
      slope <- lm_summary$estimate[lm_summary$term == "log_area"]
      pval <- lm_summary$p.value[lm_summary$term == "log_area"]
      report_lines <- c(report_lines,
                        "Linear model: Positional SD ~ log10(Grid Area)",
                        paste0("  Slope: ", round(slope, 4), " (p = ", format.pval(pval, digits=3), ")"),
                        paste0("  R²: ", round(summary(lm_fit)$r.squared, 4)),
                        ""
      )
    }
  }
  
  report_lines <- c(report_lines,
                    "───────────────────────────────────────────────────────────────",
                    "OUTPUT FILES",
                    "───────────────────────────────────────────────────────────────",
                    "",
                    "CSV Data:",
                    paste0("  ", list.files(output_dir, pattern = "\\.csv$")),
                    "",
                    "Figures:",
                    paste0("  ", list.files(output_dir, pattern = "\\.png$")),
                    "",
                    "═══════════════════════════════════════════════════════════════"
  )
  
  report_file <- fs::path(output_dir, "analysis_report.txt")
  writeLines(report_lines, report_file)
  
  # ========================================================================
  # DONE
  # ========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  ADVANCED ANALYSIS COMPLETE                               ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  cat("Output directory:\n")
  cat("  ", output_dir, "\n")
  cat("\n")
  cat("Key findings:\n")
  
  if (!is.null(nefin) && exists("impact_summary")) {
    overall_reduction <- mean(impact_summary$mean_reduction_pct)
    cat("  • NEFIN reduces sampling error by ", round(overall_reduction, 1), "% on average\n", sep="")
  }
  
  if (exists("scale_summary")) {
    finest <- scale_summary$grid_scale[which.min(scale_summary$grid_area_acres)]
    coarsest <- scale_summary$grid_scale[which.max(scale_summary$grid_area_acres)]
    cat("  • Finest grid (", finest, ") has lowest positional SD\n", sep="")
    cat("  • Coarsest grid (", coarsest, ") has highest positional SD\n", sep="")
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
  
  advanced_analysis(
    consolidated_dir = get_arg("--dir", NULL),
    output_dir = get_arg("--out", NULL)
  )
}