#!/usr/bin/env Rscript
# =============================================================================
# nefin_dominance_analysis.R
# Analyze potential bias risk from NEFIN source dominance in hexes
# 
# Run: Rscript R/nefin_dominance_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  NEFIN Dominance & Bias Risk Analysis                               ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# =============================================================================
# Configuration
# =============================================================================
COMPARISON_DATA <- "runs/consolidated_20251211_120654/fia_nefin_comparison_all_scales.csv"
OUTPUT_DIR <- "runs/consolidated_20251211_120654/dominance_analysis"

DOMINANCE_THRESHOLD <- 0.50   # Flag if single source contributes >50% of plots
ESTIMATE_SHIFT_THRESHOLD <- 10  # Flag if estimate shifts >10% (in percent units)

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

scale_order <- c("100ha", "500ha", "1kha", "5kha", "10kha", "50kha", "64kha", "100kha")

order_scales <- function(df, col = "grid_scale") {
  if (col %in% names(df)) {
    df[[col]] <- factor(df[[col]], levels = intersect(scale_order, unique(df[[col]])), ordered = TRUE)
  }
  df
}

# =============================================================================
# Load Data
# =============================================================================
cat("Loading data...\n")

if (!file.exists(COMPARISON_DATA)) {
  stop("Comparison file not found: ", COMPARISON_DATA)
}

comparison <- read_csv(COMPARISON_DATA, show_col_types = FALSE)
cat("  ✓ Loaded", nrow(comparison), "hex-scale observations\n")

# Identify key columns
cat("\nIdentifying columns...\n")

# Plot count columns
n_fia_col <- if ("n_plots_fia" %in% names(comparison)) "n_plots_fia" else "fia_only_n"
n_nefin_col <- if ("n_plots_nefin" %in% names(comparison)) "n_plots_nefin" else "nefin_n_5y"

# Percent difference column
pct_col <- if ("est_change_pct" %in% names(comparison)) "est_change_pct" else "pct_diff"

cat("  FIA plot count column:", n_fia_col, "\n")
cat("  NEFIN plot count column:", n_nefin_col, "\n")
cat("  Percent difference column:", pct_col, "\n\n")

# =============================================================================
# Analysis 1: Plot Count Dominance by Scale
# =============================================================================
cat("Analyzing plot count dominance...\n")

plot_dominance <- comparison %>%
  mutate(
    n_fia = .data[[n_fia_col]],
    n_nefin = .data[[n_nefin_col]],
    pct_diff_val = .data[[pct_col]]
  ) %>%
  filter(!is.na(n_fia) | !is.na(n_nefin)) %>%
  mutate(
    n_fia = replace_na(n_fia, 0),
    n_nefin = replace_na(n_nefin, 0),
    n_total = n_fia + n_nefin,
    pct_fia = ifelse(n_total > 0, n_fia / n_total * 100, NA),
    pct_nefin = ifelse(n_total > 0, n_nefin / n_total * 100, NA),
    
    # Dominance flags
    fia_dominant = pct_fia > DOMINANCE_THRESHOLD * 100,
    nefin_dominant = pct_nefin > DOMINANCE_THRESHOLD * 100,
    balanced = !fia_dominant & !nefin_dominant,
    
    # Single source flags
    fia_only = n_nefin == 0 & n_fia > 0,
    nefin_only = n_fia == 0 & n_nefin > 0
  )

# Summarize by scale
dominance_summary <- plot_dominance %>%
  group_by(grid_scale) %>%
  summarize(
    n_hexes = n(),
    n_with_both = sum(n_fia > 0 & n_nefin > 0, na.rm = TRUE),
    n_fia_only = sum(fia_only, na.rm = TRUE),
    n_nefin_only = sum(nefin_only, na.rm = TRUE),
    n_fia_dominant = sum(fia_dominant & !fia_only, na.rm = TRUE),
    n_nefin_dominant = sum(nefin_dominant & !nefin_only, na.rm = TRUE),
    n_balanced = sum(balanced, na.rm = TRUE),
    pct_fia_dominant = mean(fia_dominant, na.rm = TRUE) * 100,
    pct_nefin_dominant = mean(nefin_dominant, na.rm = TRUE) * 100,
    pct_balanced = mean(balanced, na.rm = TRUE) * 100,
    mean_pct_nefin = mean(pct_nefin, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  order_scales()

cat("  ✓ Plot dominance analysis complete\n")

# =============================================================================
# Analysis 2: Estimate Shift vs Dominance Relationship
# =============================================================================
cat("Analyzing estimate shift vs dominance...\n")

shift_dominance <- plot_dominance %>%
  filter(!is.na(pct_diff_val)) %>%
  mutate(
    large_shift = abs(pct_diff_val) > ESTIMATE_SHIFT_THRESHOLD,
    
    # Risk category
    bias_risk = case_when(
      nefin_dominant & large_shift ~ "HIGH_RISK",
      nefin_dominant & !large_shift ~ "MODERATE_RISK",
      !nefin_dominant & large_shift ~ "INVESTIGATE",
      TRUE ~ "LOW_RISK"
    )
  )

# Summarize shift-dominance relationship
shift_summary <- shift_dominance %>%
  group_by(grid_scale) %>%
  summarize(
    n_hexes = n(),
    n_large_shift = sum(large_shift, na.rm = TRUE),
    pct_large_shift = mean(large_shift, na.rm = TRUE) * 100,
    n_high_risk = sum(bias_risk == "HIGH_RISK", na.rm = TRUE),
    n_moderate_risk = sum(bias_risk == "MODERATE_RISK", na.rm = TRUE),
    n_investigate = sum(bias_risk == "INVESTIGATE", na.rm = TRUE),
    n_low_risk = sum(bias_risk == "LOW_RISK", na.rm = TRUE),
    cor_nefin_pct_shift = cor(pct_nefin, abs(pct_diff_val), use = "complete.obs"),
    .groups = "drop"
  ) %>%
  order_scales()

cat("  ✓ Shift-dominance analysis complete\n")

# =============================================================================
# Analysis 3: High-Risk Hex Identification
# =============================================================================
cat("Identifying high-risk hexes...\n")

high_risk_hexes <- shift_dominance %>%
  filter(bias_risk %in% c("HIGH_RISK", "MODERATE_RISK")) %>%
  select(
    hex_id, grid_scale, 
    n_fia, n_nefin, pct_nefin,
    pct_diff_val, bias_risk,
    any_of(c("recommendation", "max_source_frac"))
  ) %>%
  arrange(bias_risk, desc(abs(pct_diff_val)))

cat("  ✓ Identified", nrow(high_risk_hexes), "high/moderate risk hexes\n")

# =============================================================================
# Analysis 4: NEFIN Source Concentration
# =============================================================================
cat("Analyzing NEFIN source concentration...\n")

if ("max_source_frac" %in% names(comparison)) {
  source_concentration <- comparison %>%
    filter(!is.na(max_source_frac), .data[[n_nefin_col]] > 0) %>%
    group_by(grid_scale) %>%
    summarize(
      n_hexes_with_nefin = n(),
      mean_max_source_frac = mean(max_source_frac, na.rm = TRUE),
      median_max_source_frac = median(max_source_frac, na.rm = TRUE),
      pct_single_source_gt50 = mean(max_source_frac > 0.5, na.rm = TRUE) * 100,
      pct_single_source_gt75 = mean(max_source_frac > 0.75, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    order_scales()
  
  cat("  ✓ Source concentration analysis complete\n")
} else {
  source_concentration <- NULL
  cat("  ⚠ max_source_frac not available\n")
}

# =============================================================================
# Save Results
# =============================================================================
cat("\nSaving results...\n")

write_csv(dominance_summary, file.path(OUTPUT_DIR, "plot_dominance_by_scale.csv"))
cat("  ✓ plot_dominance_by_scale.csv\n")

write_csv(shift_summary, file.path(OUTPUT_DIR, "shift_dominance_summary.csv"))
cat("  ✓ shift_dominance_summary.csv\n")

if (nrow(high_risk_hexes) > 0) {
  write_csv(high_risk_hexes, file.path(OUTPUT_DIR, "high_risk_hexes.csv"))
  cat("  ✓ high_risk_hexes.csv (", nrow(high_risk_hexes), " hexes)\n")
}

if (!is.null(source_concentration)) {
  write_csv(source_concentration, file.path(OUTPUT_DIR, "nefin_source_concentration.csv"))
  cat("  ✓ nefin_source_concentration.csv\n")
}

write_csv(shift_dominance, file.path(OUTPUT_DIR, "hex_bias_risk_detail.csv"))
cat("  ✓ hex_bias_risk_detail.csv\n")

# =============================================================================
# Generate Figures
# =============================================================================
cat("\nGenerating figures...\n")

# Figure 1: Dominance breakdown by scale
dom_long <- dominance_summary %>%
  select(grid_scale, n_fia_only, n_nefin_only, n_fia_dominant, n_nefin_dominant, n_balanced) %>%
  pivot_longer(-grid_scale, names_to = "category", values_to = "n") %>%
  mutate(
    category = factor(category, 
      levels = c("n_fia_only", "n_fia_dominant", "n_balanced", "n_nefin_dominant", "n_nefin_only"),
      labels = c("FIA Only", "FIA Dominant", "Balanced", "NEFIN Dominant", "NEFIN Only"))
  )

p1 <- ggplot(dom_long, aes(x = grid_scale, y = n, fill = category)) +
  geom_col(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = c(
    "FIA Only" = "#e74c3c",
    "FIA Dominant" = "#f39c12", 
    "Balanced" = "#27ae60",
    "NEFIN Dominant" = "#3498db",
    "NEFIN Only" = "#9b59b6"
  ), name = "Plot Source\nBalance") +
  labs(
    title = "Hex Plot Source Balance by Scale",
    subtitle = paste0("Dominant = >", DOMINANCE_THRESHOLD * 100, "% of plots from single source"),
    x = "Hex Grid Scale",
    y = "Number of Hexes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(OUTPUT_DIR, "fig1_dominance_by_scale.png"), p1,
       width = 10, height = 6, dpi = 300)
cat("  ✓ fig1_dominance_by_scale.png\n")

# Figure 2: NEFIN fraction vs estimate shift
p2 <- shift_dominance %>%
  filter(n_nefin > 0) %>%
  ggplot(aes(x = pct_nefin, y = abs(pct_diff_val), color = bias_risk)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = TRUE, color = "gray40", linetype = "dashed") +
  geom_hline(yintercept = ESTIMATE_SHIFT_THRESHOLD, linetype = "dashed", color = "#e74c3c") +
  geom_vline(xintercept = DOMINANCE_THRESHOLD * 100, linetype = "dashed", color = "#e74c3c") +
  scale_color_manual(values = c(
    "HIGH_RISK" = "#e74c3c",
    "MODERATE_RISK" = "#f39c12",
    "INVESTIGATE" = "#3498db",
    "LOW_RISK" = "#27ae60"
  ), name = "Bias Risk") +
  facet_wrap(~ grid_scale, ncol = 4, scales = "free") +
  labs(
    title = "NEFIN Plot Fraction vs Estimate Shift",
    subtitle = paste0("Dashed lines: ", ESTIMATE_SHIFT_THRESHOLD, "% shift & ", 
                      DOMINANCE_THRESHOLD * 100, "% dominance thresholds"),
    x = "NEFIN Plot Fraction (%)",
    y = "|% Difference in Mean|"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(OUTPUT_DIR, "fig2_nefin_fraction_vs_shift.png"), p2,
       width = 12, height = 8, dpi = 300)
cat("  ✓ fig2_nefin_fraction_vs_shift.png\n")

# Figure 3: Risk category distribution by scale
risk_long <- shift_summary %>%
  select(grid_scale, n_high_risk, n_moderate_risk, n_investigate, n_low_risk) %>%
  pivot_longer(-grid_scale, names_to = "risk", values_to = "n") %>%
  mutate(
    risk = factor(risk,
      levels = c("n_high_risk", "n_moderate_risk", "n_investigate", "n_low_risk"),
      labels = c("High Risk", "Moderate Risk", "Investigate", "Low Risk"))
  )

p3 <- ggplot(risk_long, aes(x = grid_scale, y = n, fill = risk)) +
  geom_col(position = "fill", alpha = 0.85) +
  scale_fill_manual(values = c(
    "High Risk" = "#e74c3c",
    "Moderate Risk" = "#f39c12",
    "Investigate" = "#3498db",
    "Low Risk" = "#27ae60"
  ), name = "Bias Risk\nCategory") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Bias Risk Category Distribution by Scale",
    subtitle = paste0("High Risk = NEFIN dominant + large estimate shift (>", ESTIMATE_SHIFT_THRESHOLD, "%)"),
    x = "Hex Grid Scale",
    y = "Proportion of Hexes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(OUTPUT_DIR, "fig3_bias_risk_by_scale.png"), p3,
       width = 10, height = 6, dpi = 300)
cat("  ✓ fig3_bias_risk_by_scale.png\n")

# =============================================================================
# Generate Report
# =============================================================================
cat("\nGenerating dominance analysis report...\n")

report_lines <- c(
  "# NEFIN Dominance & Bias Risk Analysis Report",
  "",
  paste("Generated:", Sys.time()),
  paste("Dominance threshold:", DOMINANCE_THRESHOLD * 100, "%"),
  paste("Estimate shift threshold:", ESTIMATE_SHIFT_THRESHOLD, "%"),
  "",
  "## Summary",
  "",
  "### Plot Balance by Scale",
  "",
  "| Scale | n_hexes | FIA Only | NEFIN Only | FIA Dom | NEFIN Dom | Balanced |",
  "|-------|---------|----------|------------|---------|-----------|----------|"
)

for (i in 1:nrow(dominance_summary)) {
  row <- dominance_summary[i, ]
  report_lines <- c(report_lines,
    sprintf("| %s | %d | %d | %d | %d | %d | %d |",
            as.character(row$grid_scale), row$n_hexes, row$n_fia_only, row$n_nefin_only,
            row$n_fia_dominant, row$n_nefin_dominant, row$n_balanced)
  )
}

report_lines <- c(report_lines, "",
  "### Bias Risk Summary by Scale",
  "",
  "| Scale | High Risk | Moderate | Investigate | Low Risk | Corr |",
  "|-------|-----------|----------|-------------|----------|------|"
)

for (i in 1:nrow(shift_summary)) {
  row <- shift_summary[i, ]
  report_lines <- c(report_lines,
    sprintf("| %s | %d | %d | %d | %d | %.3f |",
            as.character(row$grid_scale), row$n_high_risk, row$n_moderate_risk, 
            row$n_investigate, row$n_low_risk, row$cor_nefin_pct_shift)
  )
}

if (nrow(high_risk_hexes) > 0) {
  report_lines <- c(report_lines, "",
    "### High-Risk Hexes (Top 10 by Estimate Shift)",
    "",
    paste("Total high/moderate risk hexes:", nrow(high_risk_hexes)),
    ""
  )
  
  top10 <- high_risk_hexes %>% head(10)
  for (i in 1:nrow(top10)) {
    row <- top10[i, ]
    report_lines <- c(report_lines,
      sprintf("- %s (%s): %.1f%% NEFIN plots, %.1f%% shift, %s",
              row$hex_id, as.character(row$grid_scale), row$pct_nefin, row$pct_diff_val, row$bias_risk)
    )
  }
}

report_lines <- c(report_lines, "",
  "## Interpretation",
  "",
  "- **High Risk**: NEFIN dominant (>50% of plots) AND large estimate shift (>10%)",
  "- **Moderate Risk**: NEFIN dominant but estimate shift within bounds",
  "- **Investigate**: Large shift but FIA dominant (other factors causing shift)",
  "- **Low Risk**: Balanced plots or small shifts",
  "",
  "Hexes flagged as High Risk warrant manual review before including NEFIN",
  "data in regional summaries.",
  ""
)

writeLines(report_lines, file.path(OUTPUT_DIR, "dominance_analysis_report.txt"))
cat("  ✓ dominance_analysis_report.txt\n")

# =============================================================================
# Final Summary
# =============================================================================
cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("DOMINANCE ANALYSIS COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("\nOutput directory:", OUTPUT_DIR, "\n")
cat("Files generated:", length(list.files(OUTPUT_DIR)), "\n\n")

cat("High-Risk Hexes by Scale:\n")
for (i in 1:nrow(shift_summary)) {
  cat(sprintf("  %8s: %3d high-risk, %3d moderate-risk, %3d investigate\n",
              as.character(shift_summary$grid_scale[i]),
              shift_summary$n_high_risk[i],
              shift_summary$n_moderate_risk[i],
              shift_summary$n_investigate[i]))
}

cat("\nCorrelation (NEFIN % vs |Estimate Shift|) by Scale:\n")
for (i in 1:nrow(shift_summary)) {
  cat(sprintf("  %8s: r = %+.3f\n",
              as.character(shift_summary$grid_scale[i]),
              shift_summary$cor_nefin_pct_shift[i]))
}
cat("\n")
