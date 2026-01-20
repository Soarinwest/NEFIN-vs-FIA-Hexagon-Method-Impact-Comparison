# R/advanced_analysis.R
# Advanced statistical analysis of FIA error and NEFIN comparison
# FIXED: Proper NA handling for division, correlation, and classification
# FIXED: Scale naming (fia → 64kha) for proper ordering in figures

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(tidyr)
  library(broom); library(scales); library(viridis); library(patchwork); library(fs)
})

# Source scale name utilities
if (file.exists("R/utils_scale_names.R")) {
  source("R/utils_scale_names.R")
} else {
  # Inline fallback if utils not available
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
      df[[scale_col]] <- factor(df[[scale_col]], levels = ordered_levels, ordered = TRUE)
    }
    df
  }
}

# =============================================================================
# HELPER FUNCTIONS FOR SAFE CALCULATIONS
# =============================================================================

# Safe division - returns NA when denominator is 0, NA, or would produce Inf/NaN
safe_fraction <- function(num, denom, default = NA_real_) {
  result <- num / denom
  result[!is.finite(result) | is.na(denom) | denom <= 0] <- default
  result
}

# Safe percentage reduction calculation
# Returns 100 * (1 - reduced/original), NA when original <= 0
safe_pct_reduction <- function(original, reduced, default = NA_real_) {
  result <- 100 * (1 - reduced / original)
  result[!is.finite(result) | is.na(original) | original <= 0] <- default
  result
}

# Safe correlation - requires minimum number of pairs
safe_cor <- function(x, y, min_pairs = 5) {
  valid <- !is.na(x) & !is.na(y)
  if (sum(valid) >= min_pairs) {
    cor(x[valid], y[valid])
  } else {
    NA_real_
  }
}

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
  
  # Standardize scale names (fia → 64kha) and order properly
  fia <- order_scales(fia, "grid_scale")
  
  cat("  ✓ FIA results loaded: ", nrow(fia), " rows\n")
  
  has_nefin <- fs::file_exists(nefin_file)
  nefin <- if (has_nefin) {
    df <- readr::read_csv(nefin_file, show_col_types = FALSE)
    order_scales(df, "grid_scale")
  } else NULL
  
  if (has_nefin) cat("  ✓ NEFIN comparison loaded: ", nrow(nefin), " rows\n")
  
  # ==========================================================================
  # DATA VALIDATION AND CLEANING
  # ==========================================================================
  cat("\n→ Validating data quality...\n")
  
  # Check FIA data quality by scale
  fia_quality <- fia |>
    dplyr::group_by(grid_scale) |>
    dplyr::summarise(
      n_total = dplyr::n(),
      n_na_se = sum(is.na(se)),
      n_zero_se = sum(se == 0, na.rm = TRUE),
      n_invalid_total_sd = sum(is.na(total_sd) | total_sd <= 0, na.rm = TRUE),
      n_single_plot = sum(n_plots == 1, na.rm = TRUE),
      pct_valid = 100 * sum(
        !is.na(se) & is.finite(se) & 
        !is.na(total_sd) & is.finite(total_sd) & total_sd > 0, 
        na.rm = TRUE
      ) / dplyr::n(),
      .groups = "drop"
    )
  
  cat("\n  FIA Data Quality by Scale:\n")
  print(as.data.frame(fia_quality), row.names = FALSE)
  
  # Warn about problematic scales
  problem_scales <- fia_quality |> dplyr::filter(pct_valid < 80)
  if (nrow(problem_scales) > 0) {
    cat("\n  ⚠ WARNING: Scales with >20% invalid rows (often single-plot hexes):\n")
    for (i in seq_len(nrow(problem_scales))) {
      cat("    -", problem_scales$grid_scale[i], ":", 
          round(100 - problem_scales$pct_valid[i], 1), "% invalid,",
          problem_scales$n_single_plot[i], "single-plot hexes\n")
    }
  }
  
  # Create cleaned FIA dataset for analyses requiring valid SE/total_sd
  fia_clean <- fia |>
    dplyr::filter(
      is.finite(se),
      is.finite(positional_sd),
      is.finite(total_sd),
      total_sd > 0
    )
  
  cat("\n  FIA rows after cleaning: ", nrow(fia_clean), " of ", nrow(fia), 
      " (", round(100 * nrow(fia_clean) / nrow(fia), 1), "%)\n")
  
  # ========================================================================
  # ANALYSIS 1: Does NEFIN reduce sampling error?
  # ========================================================================
  
  if (!is.null(nefin)) {
    cat("\n")
    cat("══════════════════════════════════════════════════════════\n")
    cat("ANALYSIS 1: NEFIN Impact on Sampling Error\n")
    cat("══════════════════════════════════════════════════════════\n")
    
    # FIXED: Filter to rows with valid SE (must be > 0 to calculate % reduction)
    nefin_impact <- nefin |>
      dplyr::filter(
        has_both, 
        !is.na(se_fia), 
        is.finite(se_fia),
        se_fia > 0,              # CRITICAL: Can't calculate % reduction from 0
        !is.na(n_plots_fia), 
        n_plots_fia >= 1,
        !is.na(n_plots_nefin),
        n_plots_nefin >= 1
      ) |>
      dplyr::mutate(
        n_combined = n_plots_fia + n_plots_nefin,
        
        # Expected SE reduction from larger sample size
        # SE is proportional to 1/sqrt(n)
        se_combined_expected = se_fia * sqrt(n_plots_fia / n_combined),
        
        # FIXED: Use safe percentage calculation
        se_reduction_pct = safe_pct_reduction(se_fia, se_combined_expected),
        
        # Absolute reduction (always calculable)
        se_reduction_abs = se_fia - se_combined_expected
      )
    
    cat("\n  Rows with valid SE for reduction analysis: ", nrow(nefin_impact), "\n")
    
    # Summary by scale
    impact_summary <- nefin_impact |>
      dplyr::group_by(grid_scale) |>
      dplyr::summarise(
        n_hexes = dplyr::n(),
        mean_fia_plots = mean(n_plots_fia, na.rm = TRUE),
        mean_nefin_plots = mean(n_plots_nefin, na.rm = TRUE),
        mean_se_fia = mean(se_fia, na.rm = TRUE),
        mean_se_combined = mean(se_combined_expected, na.rm = TRUE),
        # Only calculate mean if we have valid values
        mean_reduction_pct = if(sum(!is.na(se_reduction_pct)) > 0) {
          mean(se_reduction_pct, na.rm = TRUE)
        } else NA_real_,
        median_reduction_pct = if(sum(!is.na(se_reduction_pct)) > 0) {
          median(se_reduction_pct, na.rm = TRUE)
        } else NA_real_,
        mean_reduction_abs = mean(se_reduction_abs, na.rm = TRUE),
        .groups = "drop"
      )
    
    readr::write_csv(impact_summary, fs::path(output_dir, "nefin_se_reduction_summary.csv"))
    
    cat("\n→ SE Reduction Summary by Scale:\n")
    print(as.data.frame(impact_summary), row.names = FALSE)
    
    # Plot 1: SE reduction distribution (only scales with data)
    plot_data <- nefin_impact |>
      dplyr::filter(!is.na(se_reduction_pct), is.finite(se_reduction_pct))
    
    if (nrow(plot_data) > 0) {
      p1 <- ggplot(plot_data, aes(x = grid_scale, y = se_reduction_pct, fill = grid_scale)) +
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
    } else {
      cat("  ⚠ No valid data for SE reduction plot\n")
    }
    
    # Plot 2: Relationship between NEFIN plots and SE reduction
    if (nrow(plot_data) > 0) {
      p2 <- ggplot(plot_data, aes(x = n_plots_nefin, y = se_reduction_pct, color = grid_scale)) +
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
  }
  
  # ========================================================================
  # ANALYSIS 2: Scale effects on positional uncertainty
  # ========================================================================
  
  cat("\n")
  cat("══════════════════════════════════════════════════════════\n")
  cat("ANALYSIS 2: Grid Scale Effects on Positional Error\n")
  cat("══════════════════════════════════════════════════════════\n")
  
  # Map grid names to approximate areas
  scale_effects <- fia_clean |>
    dplyr::mutate(
      grid_area_ha = dplyr::case_when(
        grepl("100ha", grid_scale, ignore.case = TRUE) ~ 100,
        grepl("500ha", grid_scale, ignore.case = TRUE) ~ 500,
        grepl("1kha", grid_scale, ignore.case = TRUE) ~ 1000,
        grepl("5kha", grid_scale, ignore.case = TRUE) ~ 5000,
        grepl("10kha", grid_scale, ignore.case = TRUE) ~ 10000,
        grepl("50kha", grid_scale, ignore.case = TRUE) ~ 50000,
        grepl("100kha", grid_scale, ignore.case = TRUE) ~ 100000,
        grepl("fia", grid_scale, ignore.case = TRUE) ~ 6000,
        grepl("1.5", grid_scale) ~ 1500,
        grepl("3k", grid_scale) ~ 3000,
        grepl("6k", grid_scale) ~ 6000,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::filter(!is.na(grid_area_ha))
  
  if (nrow(scale_effects) > 0) {
    # Summary by scale
    scale_summary <- scale_effects |>
      dplyr::group_by(grid_scale, grid_area_ha) |>
      dplyr::summarise(
        n = dplyr::n(),
        mean_pos_sd = mean(positional_sd, na.rm = TRUE),
        median_pos_sd = median(positional_sd, na.rm = TRUE),
        sd_pos_sd = sd(positional_sd, na.rm = TRUE),
        mean_se = mean(se, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(grid_area_ha)
    
    readr::write_csv(scale_summary, fs::path(output_dir, "positional_sd_by_scale.csv"))
    
    cat("\n→ Positional SD by Grid Scale:\n")
    print(as.data.frame(scale_summary), row.names = FALSE)
    
    # Regression: positional SD ~ log(grid area)
    scale_model <- scale_effects |>
      dplyr::filter(is.finite(positional_sd), positional_sd > 0, grid_area_ha > 0) |>
      dplyr::mutate(log_area = log10(grid_area_ha))
    
    if (nrow(scale_model) > 10) {
      lm_fit <- lm(positional_sd ~ log_area, data = scale_model)
      lm_summary <- broom::tidy(lm_fit)
      
      readr::write_csv(lm_summary, fs::path(output_dir, "scale_effect_model.csv"))
      
      cat("\n→ Linear Model: Positional SD ~ log10(Grid Area)\n")
      print(as.data.frame(lm_summary), row.names = FALSE)
      cat("\n  R²: ", round(summary(lm_fit)$r.squared, 4), "\n")
    }
    
    # Plot 3: Scale effect
    if (nrow(scale_summary) >= 2) {
      p3 <- ggplot(scale_summary, aes(x = grid_area_ha, y = mean_pos_sd)) +
        geom_line(linewidth = 1.2, color = "steelblue") +
        geom_point(size = 4, color = "steelblue") +
        geom_errorbar(aes(ymin = pmax(0, mean_pos_sd - sd_pos_sd/sqrt(n)),
                          ymax = mean_pos_sd + sd_pos_sd/sqrt(n)),
                      width = 0.1, linewidth = 1) +
        geom_text(aes(label = grid_scale), vjust = -1.5, size = 3.5) +
        scale_x_log10(labels = scales::comma) +
        labs(title = "Positional Uncertainty vs Grid Size",
             subtitle = "Error bars show ±1 SE",
             x = "Grid Area (hectares, log scale)", y = "Mean Positional SD (Mg/ha)") +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(face = "bold", size = 14))
      
      ggsave(fs::path(output_dir, "03_scale_effect.png"), p3,
             width = 10, height = 6, dpi = 300)
    }
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
      dplyr::mutate(
        # FIA is sparse if < 5 plots (or NA/0)
        fia_sparse = is.na(n_plots_fia) | n_plots_fia < 5,
        # NEFIN is substantial if >= 3 plots
        nefin_substantial = !is.na(n_plots_nefin) & n_plots_nefin >= 3,
        # NEFIN adds value where FIA is sparse AND NEFIN is substantial
        adds_value = fia_sparse & nefin_substantial,
        
        # Potential SE improvement (only for valid SE values)
        se_improvement = dplyr::case_when(
          is.na(se_fia) | se_fia <= 0 ~ NA_real_,
          is.na(n_plots_fia) | is.na(n_plots_nefin) ~ NA_real_,
          TRUE ~ se_fia * (1 - sqrt(n_plots_fia / (n_plots_fia + n_plots_nefin)))
        )
      )
    
    # Summary by scale
    value_summary <- nefin_value |>
      dplyr::group_by(grid_scale) |>
      dplyr::summarise(
        total_hexes = dplyr::n(),
        n_sparse_fia = sum(fia_sparse, na.rm = TRUE),
        n_nefin_supplements = sum(adds_value, na.rm = TRUE),
        pct_supplemented = 100 * n_nefin_supplements / total_hexes,
        .groups = "drop"
      )
    
    cat("\n→ NEFIN Value Summary:\n")
    print(as.data.frame(value_summary), row.names = FALSE)
    
    # High-value hexes
    value_hexes <- nefin_value |>
      dplyr::filter(adds_value, !is.na(se_improvement)) |>
      dplyr::arrange(desc(se_improvement)) |>
      dplyr::select(hex_id, grid_scale, year_label, n_plots_fia, n_plots_nefin,
                    mean_fia, mean_nefin, diff, se_fia, se_improvement) |>
      dplyr::slice_head(n = 50)
    
    if (nrow(value_hexes) > 0) {
      readr::write_csv(value_hexes, fs::path(output_dir, "high_value_nefin_hexes.csv"))
      cat("\n  Top high-value NEFIN hexes saved\n")
    }
    
    # Plot 4: NEFIN value coverage
    p4 <- ggplot(value_summary, aes(x = grid_scale, y = pct_supplemented, fill = grid_scale)) +
      geom_col() +
      geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", pct_supplemented, 
                                    n_nefin_supplements, total_hexes)),
                vjust = -0.2, size = 3.5) +
      scale_fill_viridis_d(option = "C") +
      coord_cartesian(ylim = c(0, max(value_summary$pct_supplemented, na.rm = TRUE) * 1.3)) +
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
  
  years <- unique(fia_clean$year_label)
  
  if (length(years) > 1) {
    cat("\n")
    cat("══════════════════════════════════════════════════════════\n")
    cat("ANALYSIS 4: Temporal Trends\n")
    cat("══════════════════════════════════════════════════════════\n")
    
    # FIA temporal trends
    fia_temporal <- fia_clean |>
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
      facet_wrap(~grid_scale, ncol = 2, scales = "free_y") +
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
      # FIXED: Use safe correlation with minimum pairs requirement
      nefin_temporal <- nefin |>
        dplyr::filter(has_both) |>
        dplyr::group_by(year_label, grid_scale) |>
        dplyr::summarise(
          n_hexes = dplyr::n(),
          n_valid_pairs = sum(!is.na(mean_fia) & !is.na(mean_nefin)),
          # FIXED: Only calculate correlation with enough data points
          correlation = if (sum(!is.na(mean_fia) & !is.na(mean_nefin)) >= 5) {
            cor(mean_fia, mean_nefin, use = "complete.obs")
          } else {
            NA_real_
          },
          rmse = sqrt(mean(diff^2, na.rm = TRUE)),
          mean_bias = mean(diff, na.rm = TRUE),
          .groups = "drop"
        )
      
      readr::write_csv(nefin_temporal, fs::path(output_dir, "nefin_temporal_agreement.csv"))
      
      # Plot 6: Temporal agreement (only scales/years with enough data)
      plot_data <- nefin_temporal |>
        dplyr::filter(!is.na(correlation), n_valid_pairs >= 5)
      
      if (nrow(plot_data) > 0) {
        p6 <- ggplot(plot_data, aes(x = year_label, y = correlation, 
                                    color = grid_scale, group = grid_scale)) +
          geom_line(linewidth = 1.2) +
          geom_point(size = 3) +
          scale_color_viridis_d(option = "D") +
          labs(title = "FIA-NEFIN Agreement Over Time",
               subtitle = "Only scale-years with ≥5 paired observations shown",
               x = "Year", y = "Correlation", color = "Grid Scale") +
          theme_minimal(base_size = 12) +
          theme(plot.title = element_text(face = "bold", size = 14),
                legend.position = "bottom")
        
        ggsave(fs::path(output_dir, "06_temporal_agreement.png"), p6,
               width = 10, height = 6, dpi = 300)
      } else {
        cat("  ⚠ Insufficient data for temporal agreement plot\n")
      }
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
      dplyr::filter(has_both, !is.na(diff)) |>
      dplyr::mutate(abs_diff = abs(diff)) |>
      dplyr::group_by(grid_scale) |>
      dplyr::arrange(desc(abs_diff)) |>
      dplyr::slice_head(n = 20) |>
      dplyr::ungroup() |>
      dplyr::select(hex_id, grid_scale, year_label, mean_fia, mean_nefin, diff, abs_diff,
                    n_plots_fia, n_plots_nefin, se_fia, se_nefin)
    
    readr::write_csv(top_disagreements, fs::path(output_dir, "top_disagreement_hexes.csv"))
    
    cat("\n→ Top disagreements saved to CSV\n")
    
    # Plot 7: Top disagreements
    if (nrow(top_disagreements) > 0) {
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
  }
  
  # ========================================================================
  # ANALYSIS 6: Error decomposition statistics
  # ========================================================================
  
  cat("\n")
  cat("══════════════════════════════════════════════════════════\n")
  cat("ANALYSIS 6: Detailed Error Decomposition\n")
  cat("══════════════════════════════════════════════════════════\n")
  
  # FIXED: Use cleaned data and safe fraction calculation
  error_decomp <- fia_clean |>
    dplyr::mutate(
      # FIXED: Use safe division
      pos_fraction = safe_fraction(positional_sd, total_sd),
      sampling_fraction = safe_fraction(se, total_sd),
      
      # FIXED: Handle NA explicitly in classification
      dominant_error = dplyr::case_when(
        is.na(pos_fraction) | is.na(sampling_fraction) ~ NA_character_,
        pos_fraction > 0.6 ~ "Positional",
        sampling_fraction > 0.6 ~ "Sampling",
        TRUE ~ "Mixed"
      )
    ) |>
    # Remove unclassifiable rows
    dplyr::filter(!is.na(dominant_error)) |>
    dplyr::group_by(grid_scale, dominant_error) |>
    dplyr::summarise(
      n_hexes = dplyr::n(),
      mean_pos_sd = mean(positional_sd, na.rm = TRUE),
      mean_se = mean(se, na.rm = TRUE),
      mean_total_sd = mean(total_sd, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::group_by(grid_scale) |>
    dplyr::mutate(pct_hexes = 100 * n_hexes / sum(n_hexes)) |>
    dplyr::ungroup()
  
  readr::write_csv(error_decomp, fs::path(output_dir, "error_decomposition_detailed.csv"))
  
  cat("\n→ Error Source Classification:\n")
  print(as.data.frame(error_decomp), row.names = FALSE)
  
  # Plot 8: Dominant error source
  if (nrow(error_decomp) > 0) {
    p8 <- ggplot(error_decomp, aes(x = grid_scale, y = pct_hexes, fill = dominant_error)) +
      geom_col(position = "stack") +
      geom_text(aes(label = sprintf("%.1f%%", pct_hexes)),
                position = position_stack(vjust = 0.5),
                color = "white", fontface = "bold", size = 3.5) +
      scale_fill_manual(values = c("Positional" = "#E63946",
                                   "Sampling" = "#457B9D",
                                   "Mixed" = "#A8DADC")) +
      labs(title = "Dominant Error Source by Grid Scale",
           subtitle = "Positional >60% of total | Sampling >60% | Mixed otherwise",
           x = "Grid Scale", y = "% of Hexes", fill = "Dominant Error") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 14),
            legend.position = "bottom")
    
    ggsave(fs::path(output_dir, "08_dominant_error_source.png"), p8,
           width = 10, height = 6, dpi = 300)
  }
  
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
    "",
    "───────────────────────────────────────────────────────────────",
    "DATA QUALITY NOTES",
    "───────────────────────────────────────────────────────────────",
    "",
    paste("Total FIA rows:", nrow(fia)),
    paste("Valid rows for analysis:", nrow(fia_clean)),
    paste("Rows excluded (invalid SE/total_sd):", nrow(fia) - nrow(fia_clean)),
    "",
    "Note: Single-plot hexes cannot have SE calculated (SE = NA).",
    "These are excluded from analyses requiring valid SE values.",
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
      reduction_str <- if (!is.na(row$mean_reduction_pct)) {
        sprintf("%.1f%%", row$mean_reduction_pct)
      } else {
        "NA (insufficient valid SE values)"
      }
      report_lines <- c(report_lines,
                        paste0("Grid: ", row$grid_scale),
                        paste0("  Mean SE reduction: ", reduction_str),
                        paste0("  Absolute reduction: ", round(row$mean_reduction_abs, 3), " Mg/ha"),
                        paste0("  FIA plots/hex: ", round(row$mean_fia_plots, 1)),
                        paste0("  NEFIN plots/hex: ", round(row$mean_nefin_plots, 1)),
                        ""
      )
    }
  }
  
  if (exists("scale_summary") && nrow(scale_summary) > 0) {
    report_lines <- c(report_lines,
                      "───────────────────────────────────────────────────────────────",
                      "2. GRID SCALE EFFECTS ON POSITIONAL ERROR",
                      "───────────────────────────────────────────────────────────────",
                      ""
    )
    
    for (i in seq_len(nrow(scale_summary))) {
      row <- scale_summary[i,]
      report_lines <- c(report_lines,
                        paste0("Grid: ", row$grid_scale, " (", format(row$grid_area_ha, big.mark=","), " ha)"),
                        paste0("  Mean positional SD: ", round(row$mean_pos_sd, 3), " Mg/ha"),
                        paste0("  Mean sampling SE: ", round(row$mean_se, 3), " Mg/ha"),
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
                    paste0("  ", paste(list.files(output_dir, pattern = "\\.csv$"), collapse = "\n  ")),
                    "",
                    "Figures:",
                    paste0("  ", paste(list.files(output_dir, pattern = "\\.png$"), collapse = "\n  ")),
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
