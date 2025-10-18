# R/analyze_scale_threshold.R
# Analyze error components across hex scales to find thresholds

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(tidyr)
  library(scales); library(fs); library(mgcv); library(broom)
})

analyze_scale_threshold <- function(consolidated_dir = NULL,
                                    output_dir = NULL) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Scale Threshold Analysis                                 ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Find consolidated directory
  if (is.null(consolidated_dir)) {
    all_dirs <- list.dirs("runs", recursive = FALSE)
    consol_dirs <- all_dirs[grepl("^runs/consolidated_", all_dirs)]
    if (length(consol_dirs) == 0) {
      stop("No consolidated results found.")
    }
    consolidated_dir <- consol_dirs[order(file.mtime(consol_dirs), decreasing = TRUE)][1]
  }
  
  if (is.null(output_dir)) {
    output_dir <- fs::path(consolidated_dir, "scale_threshold_analysis")
  }
  fs::dir_create(output_dir, recurse = TRUE)
  
  # Load data
  cat("Loading data from:", consolidated_dir, "\n")
  fia_file <- fs::path(consolidated_dir, "fia_all_scales.csv")
  fia <- readr::read_csv(fia_file, show_col_types = FALSE)
  
  # Add numeric scale value (extract from grid_scale name)
  fia <- fia |>
    dplyr::mutate(
      scale_ha = dplyr::case_when(      
        grid_scale == "200ha" ~ 200,
        grid_scale == "400ha" ~ 400,
        grid_scale == "800ha" ~ 800,
        grid_scale == "1.5kha" ~ 1500,
        grid_scale == "2.5kha" ~ 2500,
        grid_scale == "3.5kha" ~ 3500,
        grid_scale == "6kha" ~ 6000,
        grid_scale == "10kha" ~ 10000,
        grid_scale == "20kha" ~ 20000,
        grid_scale == "40kha" ~ 40000,
        grid_scale == "80kha" ~ 80000,
        grid_scale == "150kha" ~ 150000,
        grid_scale == "250kha" ~ 250000,
        grid_scale == "600kha" ~ 600000,
        grid_scale == "fia" ~ 640000,
        grid_scale == "1000kha" ~ 1000000,
        TRUE ~ NA_real_
      ),
      scale_km2 = scale_ha / 100  # Convert ha to km²
    )
  
  
  # ========================================================================
  # ANALYSIS 1: Error components by scale
  # ========================================================================
  
  cat("\nAnalysis 1: Error components across scales\n")
  
  error_by_scale <- fia |>
    dplyr::group_by(grid_scale, scale_acres, scale_km) |>
    dplyr::summarise(
      n_hexes = dplyr::n_distinct(hex_id),
      n_obs = dplyr::n(),
      
      # Mean error components
      mean_se = mean(se, na.rm = TRUE),
      median_se = median(se, na.rm = TRUE),
      
      mean_pos_sd = mean(positional_sd, na.rm = TRUE),
      median_pos_sd = median(positional_sd, na.rm = TRUE),
      
      mean_total_sd = mean(total_sd, na.rm = TRUE),
      median_total_sd = median(total_sd, na.rm = TRUE),
      
      # Positional fraction
      mean_pos_fraction = mean(positional_sd / total_sd, na.rm = TRUE),
      median_pos_fraction = median(positional_sd / total_sd, na.rm = TRUE),
      
      # CV of errors
      cv_se = 100 * sd(se, na.rm = TRUE) / mean(se, na.rm = TRUE),
      cv_pos = 100 * sd(positional_sd, na.rm = TRUE) / mean(positional_sd, na.rm = TRUE),
      
      # Sample size
      mean_n_plots = mean(n_plots, na.rm = TRUE),
      
      .groups = "drop"
    ) |>
    dplyr::arrange(scale_acres)
  
  readr::write_csv(error_by_scale, fs::path(output_dir, "error_by_scale.csv"))
  
  # ========================================================================
  # ANALYSIS 2: Find threshold using GAM
  # ========================================================================
  
  cat("\nAnalysis 2: Identifying scale thresholds\n")
  
  # Fit GAM to identify where positional error contribution changes
  gam_pos_fraction <- mgcv::gam(mean_pos_fraction ~ s(log(scale_km)), 
                                data = error_by_scale)
  
  # Predict at fine scale
  scale_pred <- data.frame(
    scale_km = exp(seq(log(min(error_by_scale$scale_km)), 
                       log(max(error_by_scale$scale_km)), 
                       length.out = 200))
  )
  scale_pred$pos_fraction_pred <- predict(gam_pos_fraction, scale_pred)
  
  # Find derivative to identify inflection points
  scale_pred$derivative <- c(NA, diff(scale_pred$pos_fraction_pred) / diff(log(scale_pred$scale_km)))
  
  # Find where derivative is maximum (steepest increase)
  threshold_idx <- which.max(abs(scale_pred$derivative[-1])) + 1
  threshold_km <- scale_pred$scale_km[threshold_idx]
  
  cat("\n  Identified threshold at approximately", round(threshold_km, 1), "km\n")
  
  # ========================================================================
  # PLOT 1: Main threshold plot
  # ========================================================================
  
  cat("\nCreating visualizations...\n")
  
  p1 <- ggplot() +
    # Points for actual data
    geom_point(data = error_by_scale, 
               aes(x = scale_km, y = mean_pos_fraction * 100),
               size = 4, color = "darkblue", alpha = 0.8) +
    
    # GAM smooth
    geom_line(data = scale_pred,
              aes(x = scale_km, y = pos_fraction_pred * 100),
              color = "red", linewidth = 1.2) +
    
    # Threshold line
    geom_vline(xintercept = threshold_km, 
               linetype = "dashed", color = "darkgreen", linewidth = 1) +
    
    # Threshold annotation
    annotate("text", x = threshold_km * 1.2, y = 80, 
             label = paste("Threshold:", round(threshold_km, 1), "km"),
             angle = 90, vjust = -0.5, size = 4, color = "darkgreen") +
    
    scale_x_log10(breaks = c(2, 5, 10, 20, 50, 100, 200),
                  labels = comma) +
    labs(title = "Positional Error Contribution vs Hex Scale",
         subtitle = "Identifying the scale threshold where coordinate fuzzing dominates",
         x = "Hex Scale (km across, log scale)",
         y = "Positional Error as % of Total Error") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(fs::path(output_dir, "01_scale_threshold_main.png"), p1,
         width = 10, height = 7, dpi = 300)
  
  # ========================================================================
  # PLOT 2: All error components
  # ========================================================================
  
  error_long <- error_by_scale |>
    dplyr::select(scale_km, mean_se, mean_pos_sd, mean_total_sd) |>
    tidyr::pivot_longer(cols = c(mean_se, mean_pos_sd, mean_total_sd),
                        names_to = "error_type", values_to = "value") |>
    dplyr::mutate(
      error_type = factor(error_type,
                          levels = c("mean_se", "mean_pos_sd", "mean_total_sd"),
                          labels = c("Sampling Error", "Positional SD", "Total SD"))
    )
  
  p2 <- ggplot(error_long, aes(x = scale_km, y = value, color = error_type)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(linewidth = 1.2, alpha = 0.6) +
    geom_vline(xintercept = threshold_km, 
               linetype = "dashed", color = "darkgray", linewidth = 1) +
    scale_x_log10(breaks = c(2, 5, 10, 20, 50, 100, 200),
                  labels = comma) +
    scale_color_manual(values = c("Sampling Error" = "#E69F00",
                                  "Positional SD" = "#56B4E9",
                                  "Total SD" = "#009E73")) +
    labs(title = "Error Components Across Hex Scales",
         subtitle = "How sampling, positional, and total errors change with scale",
         x = "Hex Scale (km across, log scale)",
         y = "Standard Deviation (Mg/ha)",
         color = "Error Type") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  
  ggsave(fs::path(output_dir, "02_error_components_by_scale.png"), p2,
         width = 10, height = 7, dpi = 300)
  
  # ========================================================================
  # PLOT 3: Error reduction efficiency
  # ========================================================================
  
  # Calculate error reduction relative to smallest scale
  baseline_total_sd <- error_by_scale$mean_total_sd[1]
  
  error_by_scale <- error_by_scale |>
    dplyr::mutate(
      error_reduction = (baseline_total_sd - mean_total_sd) / baseline_total_sd * 100,
      error_reduction_per_doubling = error_reduction / log2(scale_km / min(scale_km))
    )
  
  p3 <- ggplot(error_by_scale, aes(x = scale_km)) +
    geom_line(aes(y = error_reduction), color = "darkblue", linewidth = 1.2) +
    geom_point(aes(y = error_reduction), color = "darkblue", size = 3) +
    geom_vline(xintercept = threshold_km, 
               linetype = "dashed", color = "darkgreen", linewidth = 1) +
    
    # Add shaded regions
    annotate("rect", xmin = 0, xmax = threshold_km, 
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "blue") +
    annotate("rect", xmin = threshold_km, xmax = Inf, 
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "red") +
    
    # Labels for regions
    annotate("text", x = threshold_km / 3, y = max(error_by_scale$error_reduction) * 0.9,
             label = "Sampling Error\nDominated", size = 4, color = "darkblue") +
    annotate("text", x = threshold_km * 3, y = max(error_by_scale$error_reduction) * 0.9,
             label = "Positional Error\nDominated", size = 4, color = "darkred") +
    
    scale_x_log10(breaks = c(2, 5, 10, 20, 50, 100, 200),
                  labels = comma) +
    labs(title = "Error Reduction Efficiency Across Scales",
         subtitle = "Percentage reduction in total error compared to finest scale",
         x = "Hex Scale (km across, log scale)",
         y = "Error Reduction (%)") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(fs::path(output_dir, "03_error_reduction_efficiency.png"), p3,
         width = 10, height = 7, dpi = 300)
  
  # ========================================================================
  # PLOT 4: Sample size effects
  # ========================================================================
  
  p4 <- ggplot(error_by_scale, aes(x = scale_km)) +
    geom_col(aes(y = mean_n_plots), fill = "steelblue", alpha = 0.7) +
    geom_line(aes(y = mean_pos_fraction * max(mean_n_plots) * 2), 
              color = "red", linewidth = 1.2) +
    geom_point(aes(y = mean_pos_fraction * max(mean_n_plots) * 2), 
               color = "red", size = 3) +
    geom_vline(xintercept = threshold_km, 
               linetype = "dashed", color = "darkgreen", linewidth = 1) +
    scale_x_log10(breaks = c(2, 5, 10, 20, 50, 100, 200),
                  labels = comma) +
    scale_y_continuous(
      name = "Mean Plots per Hex",
      sec.axis = sec_axis(~ . / (max(error_by_scale$mean_n_plots) * 2),
                          name = "Positional Error Fraction",
                          labels = percent)
    ) +
    labs(title = "Sample Size and Positional Error Across Scales",
         subtitle = "Trade-off between plots per hex and positional uncertainty",
         x = "Hex Scale (km across, log scale)") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"),
          axis.title.y.right = element_text(color = "red"),
          axis.text.y.right = element_text(color = "red"))
  
  ggsave(fs::path(output_dir, "04_sample_size_vs_positional.png"), p4,
         width = 10, height = 7, dpi = 300)
  
  # ========================================================================
  # SUMMARY REPORT
  # ========================================================================
  
  report_lines <- c(
    "═══════════════════════════════════════════════════════════════",
    "SCALE THRESHOLD ANALYSIS REPORT",
    "═══════════════════════════════════════════════════════════════",
    "",
    paste("Generated:", Sys.time()),
    paste("Data source:", consolidated_dir),
    "",
    "───────────────────────────────────────────────────────────────",
    "KEY FINDINGS",
    "───────────────────────────────────────────────────────────────",
    "",
    paste("1. SCALE THRESHOLD:", round(threshold_km, 1), "km"),
    paste("   - Below this scale: Sampling error dominates"),
    paste("   - Above this scale: Positional error dominates"),
    "",
    "2. ERROR COMPONENTS BY SCALE:",
    ""
  )
  
  # Add table of key scales
  key_scales <- error_by_scale |>
    dplyr::filter(scale_km %in% c(2, 5, 10, 20, 50, 100)) |>
    dplyr::select(scale_km, mean_se, mean_pos_sd, mean_total_sd, mean_pos_fraction)
  
  for (i in seq_len(nrow(key_scales))) {
    row <- key_scales[i,]
    report_lines <- c(report_lines,
                      paste0("Scale: ", row$scale_km, " km"),
                      paste0("  Sampling Error: ", round(row$mean_se, 3), " Mg/ha"),
                      paste0("  Positional SD: ", round(row$mean_pos_sd, 3), " Mg/ha"),
                      paste0("  Total SD: ", round(row$mean_total_sd, 3), " Mg/ha"),
                      paste0("  Positional %: ", round(row$mean_pos_fraction * 100, 1), "%"),
                      ""
    )
  }
  
  report_lines <- c(report_lines,
                    "───────────────────────────────────────────────────────────────",
                    "RECOMMENDATIONS",
                    "───────────────────────────────────────────────────────────────",
                    "",
                    paste("1. For maximum precision, use hex scales <", round(threshold_km, 0), "km"),
                    "   where sampling error dominates and coordinate fuzzing has",
                    "   minimal impact.",
                    "",
                    paste("2. At scales >", round(threshold_km, 0), "km, coordinate fuzzing"),
                    "   becomes the dominant source of uncertainty. Consider:",
                    "   - Using exact coordinates if available",
                    "   - Implementing spatial smoothing techniques",
                    "   - Reporting larger confidence intervals",
                    "",
                    "3. The 'sweet spot' appears to be around", 
                    round(threshold_km * 0.8, 0), "-", round(threshold_km * 1.2, 0), "km",
                    "   balancing sample size, total error, and computational",
                    "   efficiency.",
                    "",
                    "═══════════════════════════════════════════════════════════════"
  )
  
  report_file <- fs::path(output_dir, "scale_threshold_report.txt")
  writeLines(report_lines, report_file)
  
  cat("\n✓ Scale threshold analysis complete!\n")
  cat("  Threshold identified at:", round(threshold_km, 1), "km\n")
  cat("  Output directory:", output_dir, "\n")
  
  invisible(list(
    threshold_km = threshold_km,
    error_by_scale = error_by_scale,
    output_dir = output_dir
  ))
}

# CLI
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  analyze_scale_threshold()
}