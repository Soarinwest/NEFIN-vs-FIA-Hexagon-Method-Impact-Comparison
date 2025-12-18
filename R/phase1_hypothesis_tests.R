#!/usr/bin/env Rscript
# =============================================================================
# phase1_hypothesis_tests.R
# Statistical tests and summary statistics for Phase I NEFIN evaluation
# 
# Run: Rscript R/phase1_hypothesis_tests.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  Phase I: Statistical Analysis of NEFIN Impact on Summaries         ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# =============================================================================
# Configuration
# =============================================================================
CONSOLIDATED_DIR <- "runs/consolidated_20251211_120654"
OUTPUT_DIR <- file.path(CONSOLIDATED_DIR, "phase1_statistics")
N_BOOTSTRAP <- 1000
SIGNIFICANCE_LEVEL <- 0.05

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# =============================================================================
# Load Data
# =============================================================================
cat("Loading data...\n")

comparison_path <- file.path(CONSOLIDATED_DIR, "fia_nefin_comparison_all_scales.csv")
if (!file.exists(comparison_path)) {
  stop("Comparison file not found: ", comparison_path)
}

comparison <- read_csv(comparison_path, show_col_types = FALSE)
cat("  Loaded", nrow(comparison), "hex-scale observations\n")
cat("  Scales:", paste(unique(comparison$grid_scale), collapse = ", "), "\n\n")

# =============================================================================
# Compute delta_se from augmented_se and se_fia
# =============================================================================
cat("Computing delta_se...\n")

# Your data has augmented_se and se_fia (or fia_only_se)
se_fia_col <- if ("se_fia" %in% names(comparison)) "se_fia" else "fia_only_se"
se_aug_col <- if ("augmented_se" %in% names(comparison)) "augmented_se" else NULL

if (is.null(se_aug_col)) {
  stop("Cannot find augmented_se column")
}

comparison <- comparison %>%
  mutate(
    delta_se = .data[[se_aug_col]] - .data[[se_fia_col]]
  )

cat("  ✓ delta_se computed from", se_aug_col, "-", se_fia_col, "\n")
cat("  delta_se range:", round(min(comparison$delta_se, na.rm = TRUE), 3), "to",
    round(max(comparison$delta_se, na.rm = TRUE), 3), "\n")
cat("  Non-NA delta_se:", sum(!is.na(comparison$delta_se)), "\n\n")

# =============================================================================
# Helper Functions
# =============================================================================
boot_ci <- function(x, stat_fn = median, n_boot = N_BOOTSTRAP, conf = 0.95) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(c(lower = NA, upper = NA, estimate = NA))
  boot_stats <- replicate(n_boot, stat_fn(sample(x, replace = TRUE)))
  alpha <- 1 - conf
  c(lower = quantile(boot_stats, alpha / 2),
    upper = quantile(boot_stats, 1 - alpha / 2),
    estimate = stat_fn(x))
}

interpret_d <- function(d) {
  if (is.na(d)) return("N/A")
  d <- abs(d)
  if (d < 0.2) return("negligible")
  if (d < 0.5) return("small")
  if (d < 0.8) return("medium")
  return("large")
}

scale_order <- c("100ha", "500ha", "1kha", "5kha", "10kha", "50kha", "64kha", "100kha")

order_scales <- function(df, col = "grid_scale") {
  if (col %in% names(df)) {
    df[[col]] <- factor(df[[col]], levels = intersect(scale_order, unique(df[[col]])), ordered = TRUE)
  }
  df
}

# =============================================================================
# 1. Scale-by-Scale Summary Statistics
# =============================================================================
cat("Computing scale-by-scale summary statistics...\n")

scale_summary <- comparison %>%
  filter(!is.na(delta_se)) %>%
  group_by(grid_scale) %>%
  summarize(
    n_hexes = n(),
    n_with_both = sum(has_both, na.rm = TRUE),
    
    # SE statistics
    mean_se_fia = mean(.data[[se_fia_col]], na.rm = TRUE),
    median_se_fia = median(.data[[se_fia_col]], na.rm = TRUE),
    mean_se_augmented = mean(.data[[se_aug_col]], na.rm = TRUE),
    median_se_augmented = median(.data[[se_aug_col]], na.rm = TRUE),
    
    # Delta SE
    mean_delta_se = mean(delta_se, na.rm = TRUE),
    median_delta_se = median(delta_se, na.rm = TRUE),
    sd_delta_se = sd(delta_se, na.rm = TRUE),
    
    # Directional counts
    n_se_reduced = sum(delta_se < 0, na.rm = TRUE),
    n_se_increased = sum(delta_se > 0, na.rm = TRUE),
    pct_se_reduced = mean(delta_se < 0, na.rm = TRUE) * 100,
    
    # Relative change
    mean_pct_delta_se = mean(delta_se / .data[[se_fia_col]] * 100, na.rm = TRUE),
    median_pct_delta_se = median(delta_se / .data[[se_fia_col]] * 100, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  order_scales()

cat("  ✓ Scale summary computed\n")

# =============================================================================
# 2. Hypothesis Tests
# =============================================================================
cat("Running hypothesis tests...\n")

hypothesis_tests <- comparison %>%
  filter(!is.na(delta_se)) %>%
  group_by(grid_scale) %>%
  summarize(
    n = n(),
    
    # Wilcoxon signed-rank test: H0: median delta_se = 0
    wilcox_p = if(n() >= 5) {
      tryCatch(
        wilcox.test(delta_se, mu = 0, alternative = "two.sided")$p.value,
        error = function(e) NA
      )
    } else NA,
    
    # One-sample t-test: H0: mean delta_se = 0
    ttest_p = if(n() >= 5) {
      tryCatch(
        t.test(delta_se, mu = 0, alternative = "two.sided")$p.value,
        error = function(e) NA
      )
    } else NA,
    
    # Sign test
    sign_test_p = if(n() >= 5) {
      n_neg <- sum(delta_se < 0, na.rm = TRUE)
      n_pos <- sum(delta_se > 0, na.rm = TRUE)
      n_total <- n_neg + n_pos
      if (n_total > 0) binom.test(n_neg, n_total, p = 0.5)$p.value else NA
    } else NA,
    
    .groups = "drop"
  ) %>%
  mutate(
    wilcox_sig = wilcox_p < SIGNIFICANCE_LEVEL,
    ttest_sig = ttest_p < SIGNIFICANCE_LEVEL,
    sign_test_sig = sign_test_p < SIGNIFICANCE_LEVEL
  ) %>%
  order_scales()

cat("  ✓ Hypothesis tests completed\n")

# =============================================================================
# 3. Effect Sizes
# =============================================================================
cat("Computing effect sizes...\n")

effect_sizes <- comparison %>%
  filter(!is.na(delta_se)) %>%
  group_by(grid_scale) %>%
  summarize(
    cohens_d = mean(delta_se, na.rm = TRUE) / sd(delta_se, na.rm = TRUE),
    relative_effect_pct = mean(delta_se / .data[[se_fia_col]] * 100, na.rm = TRUE),
    pct_meaningful_reduction = mean(delta_se / .data[[se_fia_col]] < -0.10, na.rm = TRUE) * 100,
    pct_meaningful_increase = mean(delta_se / .data[[se_fia_col]] > 0.10, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(effect_interpretation = sapply(cohens_d, interpret_d)) %>%
  order_scales()

cat("  ✓ Effect sizes computed\n")

# =============================================================================
# 4. Bootstrap Confidence Intervals
# =============================================================================
cat("Computing bootstrap confidence intervals...\n")

bootstrap_cis <- comparison %>%
  filter(!is.na(delta_se)) %>%
  group_by(grid_scale) %>%
  summarize(
    boot_result = list(boot_ci(delta_se, stat_fn = median, n_boot = N_BOOTSTRAP)),
    .groups = "drop"
  ) %>%
  mutate(
    ci_lower = sapply(boot_result, function(x) x["lower"]),
    ci_upper = sapply(boot_result, function(x) x["upper"]),
    median_estimate = sapply(boot_result, function(x) x["estimate"]),
    ci_excludes_zero = (ci_lower > 0) | (ci_upper < 0)
  ) %>%
  select(-boot_result) %>%
  order_scales()

cat("  ✓ Bootstrap CIs computed\n")

# =============================================================================
# 5. Recommendation Distribution
# =============================================================================
cat("Analyzing recommendation distribution...\n")

if ("recommendation" %in% names(comparison)) {
  recommendation_summary <- comparison %>%
    group_by(grid_scale, recommendation) %>%
    summarize(n = n(), .groups = "drop") %>%
    group_by(grid_scale) %>%
    mutate(pct = n / sum(n) * 100, total = sum(n)) %>%
    ungroup() %>%
    order_scales()
  
  recommendation_wide <- recommendation_summary %>%
    select(grid_scale, recommendation, pct) %>%
    pivot_wider(names_from = recommendation, values_from = pct, values_fill = 0)
  
  cat("  ✓ Recommendation analysis completed\n")
} else {
  recommendation_summary <- NULL
  recommendation_wide <- NULL
  cat("  ⚠ No recommendation column found\n")
}

# =============================================================================
# 6. Estimate Stability Analysis
# =============================================================================
cat("Analyzing estimate stability...\n")

# Use est_change_pct or pct_diff
pct_col <- if ("est_change_pct" %in% names(comparison)) "est_change_pct" else "pct_diff"
diff_col <- if ("est_change" %in% names(comparison)) "est_change" else "diff"

stability_summary <- comparison %>%
  filter(!is.na(.data[[pct_col]])) %>%
  group_by(grid_scale) %>%
  summarize(
    n_hexes = n(),
    mean_abs_diff = mean(abs(.data[[diff_col]]), na.rm = TRUE),
    median_abs_diff = median(abs(.data[[diff_col]]), na.rm = TRUE),
    max_abs_diff = max(abs(.data[[diff_col]]), na.rm = TRUE),
    mean_abs_pct_diff = mean(abs(.data[[pct_col]]), na.rm = TRUE),
    median_abs_pct_diff = median(abs(.data[[pct_col]]), na.rm = TRUE),
    pct_diff_gt_5 = mean(abs(.data[[pct_col]]) > 5, na.rm = TRUE) * 100,
    pct_diff_gt_10 = mean(abs(.data[[pct_col]]) > 10, na.rm = TRUE) * 100,
    pct_diff_gt_20 = mean(abs(.data[[pct_col]]) > 20, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  order_scales()

cat("  ✓ Stability analysis completed\n")

# =============================================================================
# 7. Save Results
# =============================================================================
cat("\nSaving results...\n")

write_csv(scale_summary, file.path(OUTPUT_DIR, "phase1_scale_summary.csv"))
cat("  ✓ phase1_scale_summary.csv\n")

write_csv(hypothesis_tests, file.path(OUTPUT_DIR, "phase1_hypothesis_tests.csv"))
cat("  ✓ phase1_hypothesis_tests.csv\n")

write_csv(effect_sizes, file.path(OUTPUT_DIR, "phase1_effect_sizes.csv"))
cat("  ✓ phase1_effect_sizes.csv\n")

write_csv(bootstrap_cis, file.path(OUTPUT_DIR, "phase1_bootstrap_cis.csv"))
cat("  ✓ phase1_bootstrap_cis.csv\n")

if (!is.null(recommendation_wide)) {
  write_csv(recommendation_wide, file.path(OUTPUT_DIR, "phase1_recommendation_distribution.csv"))
  cat("  ✓ phase1_recommendation_distribution.csv\n")
}

write_csv(stability_summary, file.path(OUTPUT_DIR, "phase1_stability_summary.csv"))
cat("  ✓ phase1_stability_summary.csv\n")

# =============================================================================
# 8. Generate Figures
# =============================================================================
cat("\nGenerating figures...\n")

# Figure 1: Delta SE distribution by scale (violin plot)
p1 <- comparison %>%
  filter(!is.na(delta_se)) %>%
  order_scales() %>%
  ggplot(aes(x = grid_scale, y = delta_se, fill = grid_scale)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.8) +
  scale_fill_viridis_d(option = "D", guide = "none") +
  labs(
    title = "Distribution of SE Change (FIA+NEFIN vs FIA-only)",
    subtitle = "Negative values = NEFIN reduced uncertainty | Dashed line = no change",
    x = "Hex Grid Scale",
    y = expression(Delta*"SE (Mg/ha)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(OUTPUT_DIR, "fig1_delta_se_distribution.png"), p1, 
       width = 10, height = 6, dpi = 300)
cat("  ✓ fig1_delta_se_distribution.png\n")

# Figure 2: Recommendation breakdown by scale
if (!is.null(recommendation_summary)) {
  rec_colors <- c(
    "USE_FIA_ONLY" = "#2ecc71",
    "USE_AUGMENTED" = "#3498db",
    "NEEDS_REVIEW_FOR_BIAS" = "#e74c3c",
    "MIXED_SIGNAL" = "#95a5a6"
  )
  
  p2 <- recommendation_summary %>%
    ggplot(aes(x = grid_scale, y = pct, fill = recommendation)) +
    geom_col(position = "stack", alpha = 0.85) +
    scale_fill_manual(values = rec_colors, name = "Recommendation") +
    labs(
      title = "NEFIN Recommendation Distribution by Scale",
      subtitle = "Percentage of hexes in each recommendation category",
      x = "Hex Grid Scale",
      y = "Percent of Hexes"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    guides(fill = guide_legend(nrow = 2))
  
  ggsave(file.path(OUTPUT_DIR, "fig2_recommendation_by_scale.png"), p2,
         width = 10, height = 7, dpi = 300)
  cat("  ✓ fig2_recommendation_by_scale.png\n")
}

# Figure 3: Effect size and significance summary
combined_stats <- hypothesis_tests %>%
  left_join(effect_sizes, by = "grid_scale") %>%
  left_join(bootstrap_cis, by = "grid_scale")

p3 <- combined_stats %>%
  ggplot(aes(x = grid_scale)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "gray40") +
  geom_point(aes(y = median_estimate, color = ci_excludes_zero), size = 4) +
  scale_color_manual(
    values = c("TRUE" = "#e74c3c", "FALSE" = "#3498db"),
    labels = c("TRUE" = "Significant (CI excludes 0)", "FALSE" = "Not Significant"),
    name = "95% CI"
  ) +
  labs(
    title = "Median ΔSE with 95% Bootstrap Confidence Intervals",
    subtitle = "Negative = NEFIN reduced uncertainty",
    x = "Hex Grid Scale",
    y = expression("Median "*Delta*"SE (Mg/ha)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(OUTPUT_DIR, "fig3_effect_sizes_with_ci.png"), p3,
       width = 10, height = 6, dpi = 300)
cat("  ✓ fig3_effect_sizes_with_ci.png\n")

# Figure 4: Estimate shift histogram
p4 <- comparison %>%
  filter(!is.na(.data[[pct_col]]), abs(.data[[pct_col]]) < 100) %>%
  order_scales() %>%
  ggplot(aes(x = .data[[pct_col]])) +
  geom_histogram(bins = 50, fill = "#3498db", alpha = 0.7, color = "white") +
  geom_vline(xintercept = c(-10, 10), linetype = "dashed", color = "#e74c3c") +
  facet_wrap(~ grid_scale, scales = "free_y", ncol = 4) +
  labs(
    title = "Distribution of Estimate Shifts (Augmented vs FIA-only Mean)",
    subtitle = "Red dashed lines at ±10% threshold",
    x = "Percent Difference in Mean Biomass",
    y = "Count"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(OUTPUT_DIR, "fig4_estimate_shift_distribution.png"), p4,
       width = 12, height = 8, dpi = 300)
cat("  ✓ fig4_estimate_shift_distribution.png\n")

# =============================================================================
# 9. Summary Report
# =============================================================================
cat("\nGenerating summary report...\n")

report_lines <- c(
  "# Phase I Statistical Summary: NEFIN Impact on Regional Biomass Summaries",
  "",
  paste("Generated:", Sys.time()),
  paste("Significance level:", SIGNIFICANCE_LEVEL),
  paste("Bootstrap replicates:", N_BOOTSTRAP),
  "",
  "## Overview",
  "",
  paste("Total hex-scale observations:", sum(scale_summary$n_hexes)),
  paste("Scales analyzed:", nrow(scale_summary)),
  "",
  "## Key Statistical Findings",
  ""
)

sig_scales <- hypothesis_tests %>% filter(wilcox_sig == TRUE) %>% pull(grid_scale) %>% as.character()
nonsig_scales <- hypothesis_tests %>% filter(wilcox_sig == FALSE | is.na(wilcox_sig)) %>% pull(grid_scale) %>% as.character()

report_lines <- c(report_lines,
  paste("Scales with significant delta-SE (p <", SIGNIFICANCE_LEVEL, "):", 
        if(length(sig_scales) > 0) paste(sig_scales, collapse = ", ") else "None"),
  paste("Scales with non-significant delta-SE:", 
        if(length(nonsig_scales) > 0) paste(nonsig_scales, collapse = ", ") else "None"),
  "",
  "## Scale-by-Scale Summary",
  ""
)

for (i in 1:nrow(scale_summary)) {
  row <- scale_summary[i, ]
  report_lines <- c(report_lines,
    paste0("### ", row$grid_scale),
    paste0("- n hexes: ", row$n_hexes),
    paste0("- Median delta-SE: ", round(row$median_delta_se, 3), " Mg/ha"),
    paste0("- % with reduced SE: ", round(row$pct_se_reduced, 1), "%"),
    paste0("- Mean % change in SE: ", round(row$mean_pct_delta_se, 1), "%"),
    ""
  )
}

report_lines <- c(report_lines,
  "## Effect Sizes (Cohen's d)",
  ""
)

for (i in 1:nrow(effect_sizes)) {
  row <- effect_sizes[i, ]
  report_lines <- c(report_lines,
    paste0("- ", row$grid_scale, ": d = ", round(row$cohens_d, 3), 
           " (", row$effect_interpretation, ")")
  )
}

report_lines <- c(report_lines,
  "",
  "## Interpretation",
  "",
  "- Negative delta-SE = NEFIN reduced uncertainty (improvement)",
  "- Positive delta-SE = NEFIN increased uncertainty (degradation)", 
  "- Cohen's d < 0.2 = negligible, 0.2-0.5 = small, 0.5-0.8 = medium, > 0.8 = large",
  ""
)

writeLines(report_lines, file.path(OUTPUT_DIR, "phase1_summary_report.txt"))
cat("  ✓ phase1_summary_report.txt\n")

# =============================================================================
# Final Summary
# =============================================================================
cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("PHASE I ANALYSIS COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("\nOutput directory:", OUTPUT_DIR, "\n")
cat("Files generated:", length(list.files(OUTPUT_DIR)), "\n\n")

cat("Quick Summary - % of hexes where NEFIN reduced SE:\n")
for (i in 1:nrow(scale_summary)) {
  cat(sprintf("  %8s: %5.1f%% (%d hexes)\n", 
              as.character(scale_summary$grid_scale[i]),
              scale_summary$pct_se_reduced[i],
              scale_summary$n_hexes[i]))
}

cat("\nStatistical Significance (Wilcoxon test, p <", SIGNIFICANCE_LEVEL, "):\n")
for (i in 1:nrow(hypothesis_tests)) {
  sig_marker <- if (!is.na(hypothesis_tests$wilcox_sig[i]) && hypothesis_tests$wilcox_sig[i]) "***" else ""
  cat(sprintf("  %8s: p = %.4f %s\n",
              as.character(hypothesis_tests$grid_scale[i]),
              hypothesis_tests$wilcox_p[i],
              sig_marker))
}
cat("\n")
