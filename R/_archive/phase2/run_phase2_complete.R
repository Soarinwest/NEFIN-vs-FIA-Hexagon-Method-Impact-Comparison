#!/usr/bin/env Rscript
# R/run_phase2_complete.R
# Master script: Run all Phase 2 analyses and generate publication figures
# Author: Soren Donisvitch
# Date: December 2024

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(fs)
})

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
cat("║                                                                              ║\n")
cat("║   NEFIN vs FIA: Phase 2 Complete Analysis                                    ║\n")
cat("║   Does Coordinate Fuzzing Degrade Biomass Prediction?                        ║\n")
cat("║                                                                              ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Create output directory
output_dir <- "runs/phase2_complete"
fs::dir_create(output_dir, recurse = TRUE)
fs::dir_create(fs::path(output_dir, "figures"), recurse = TRUE)

# =============================================================================
# RUN INDIVIDUAL ANALYSES
# =============================================================================

cat("Running individual analyses...\n")
cat("═══════════════════════════════════════════════════════════════════════════════\n\n")

# Track which analyses succeeded
results <- list()

# 1. Spatial Model Comparison (holdout test)
cat("1. Spatial Model Comparison...\n")
script_path <- "R/_archive/phase2/spatial_model_comparison_v2.R"
if (!fs::file_exists(script_path)) {
  script_path <- "R/_archive/phase2/spatial_model_comparison.R"
}
if (fs::file_exists(script_path)) {
  tryCatch({
    source(script_path)
    # Call the main function after sourcing
    if (exists("run_model_comparison")) {
      run_model_comparison()
    } else if (exists("run_spatial_model_comparison")) {
      run_spatial_model_comparison()
    }
    results$model_comparison <- TRUE
    cat("   ✓ Complete\n\n")
  }, error = function(e) {
    cat("   ✗ Failed:", e$message, "\n\n")
    results$model_comparison <- FALSE
  })
} else {
  cat("   ⚠ Script not found\n\n")
}

# 2. Fuzzing Effect Analysis
cat("2. Fuzzing Effect Analysis...\n")
script_path <- "R/_archive/phase2/fuzzing_effect_analysis.R"
if (fs::file_exists(script_path)) {
  tryCatch({
    source(script_path)
    # Call the main function after sourcing
    if (exists("analyze_fuzzing_effects")) {
      analyze_fuzzing_effects()
    }
    results$fuzzing <- TRUE
    cat("   ✓ Complete\n\n")
  }, error = function(e) {
    cat("   ✗ Failed:", e$message, "\n\n")
    results$fuzzing <- FALSE
  })
} else {
  cat("   ⚠ Script not found\n\n")
}

# 3. Sensor Resolution Comparison
cat("3. Sensor Resolution Comparison (MODIS vs S2)...\n")
script_path <- "R/_archive/phase2/sensor_resolution_comparison.R"
if (fs::file_exists(script_path)) {
  tryCatch({
    source(script_path)
    # Call the main function after sourcing
    if (exists("compare_sensor_resolution")) {
      compare_sensor_resolution()
    }
    results$sensor <- TRUE
    cat("   ✓ Complete\n\n")
  }, error = function(e) {
    cat("   ✗ Failed:", e$message, "\n\n")
    results$sensor <- FALSE
  })
} else {
  cat("   ⚠ Script not found\n\n")
}

# 4. Landscape Heterogeneity Analysis
cat("4. Landscape Heterogeneity Analysis...\n")
script_path <- "R/_archive/phase2/landscape_heterogeneity_analysis.R"
if (fs::file_exists(script_path)) {
  tryCatch({
    source(script_path)
    # Call the main function after sourcing
    if (exists("analyze_heterogeneity_effect")) {
      analyze_heterogeneity_effect()
    }
    results$heterogeneity <- TRUE
    cat("   ✓ Complete\n\n")
  }, error = function(e) {
    cat("   ✗ Failed:", e$message, "\n\n")
    results$heterogeneity <- FALSE
  })
} else {
  cat("   ⚠ Script not found\n\n")
}

# 5. Prediction Maps
cat("5. MODIS Prediction Maps...\n")
script_path <- "R/_archive/phase2/create_prediction_maps.R"
if (fs::file_exists(script_path)) {
  tryCatch({
    source(script_path)
    # Call the main function after sourcing
    if (exists("create_all_prediction_maps")) {
      create_all_prediction_maps()
    }
    results$maps <- TRUE
    cat("   ✓ Complete\n\n")
  }, error = function(e) {
    cat("   ✗ Failed:", e$message, "\n\n")
    results$maps <- FALSE
  })
} else {
  cat("   ⚠ Script not found\n\n")
}

# =============================================================================
# COMPILE RESULTS
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════════\n")
cat("COMPILING RESULTS\n")
cat("═══════════════════════════════════════════════════════════════════════════════\n\n")

# Load results from each analysis
compiled <- list()

# Model comparison results
if (fs::file_exists("runs/spatial_model_comparison/holdout_prediction_results.csv")) {
  compiled$holdout <- read_csv("runs/spatial_model_comparison/holdout_prediction_results.csv", 
                                show_col_types = FALSE)
  cat("  ✓ Loaded holdout results\n")
}

# Fuzzing results
if (fs::file_exists("runs/fuzzing_effect_analysis/covariate_uncertainty_summary.csv")) {
  compiled$uncertainty <- read_csv("runs/fuzzing_effect_analysis/covariate_uncertainty_summary.csv",
                                    show_col_types = FALSE)
  cat("  ✓ Loaded fuzzing uncertainty\n")
}

if (fs::file_exists("runs/fuzzing_effect_analysis/relationship_attenuation.csv")) {
  compiled$attenuation <- read_csv("runs/fuzzing_effect_analysis/relationship_attenuation.csv",
                                    show_col_types = FALSE)
  cat("  ✓ Loaded attenuation results\n")
}

# Sensor comparison
if (fs::file_exists("runs/sensor_resolution_comparison/sensor_comparison_summary.csv")) {
  compiled$sensor <- read_csv("runs/sensor_resolution_comparison/sensor_comparison_summary.csv",
                               show_col_types = FALSE)
  cat("  ✓ Loaded sensor comparison\n")
}

# Heterogeneity results
if (fs::file_exists("runs/heterogeneity_analysis/rmse_by_heterogeneity_class.csv")) {
  compiled$heterogeneity <- read_csv("runs/heterogeneity_analysis/rmse_by_heterogeneity_class.csv",
                                      show_col_types = FALSE)
  cat("  ✓ Loaded heterogeneity results\n")
}

# Hex ΔRMSE
if (fs::file_exists("runs/prediction_maps/hex_delta_rmse.csv")) {
  compiled$hex_rmse <- read_csv("runs/prediction_maps/hex_delta_rmse.csv",
                                 show_col_types = FALSE)
  cat("  ✓ Loaded hex ΔRMSE\n")
}

# =============================================================================
# CREATE SUMMARY TABLE
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════════\n")
cat("KEY FINDINGS SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════════════════\n\n")

summary_lines <- c()

# Holdout RMSE improvement
if (!is.null(compiled$holdout)) {
  xgb_results <- compiled$holdout %>% filter(model == "XGBoost")
  if (nrow(xgb_results) == 2) {
    fia_rmse <- xgb_results$holdout_rmse[xgb_results$trained_on == "FIA (fuzzed)"]
    nefin_rmse <- xgb_results$holdout_rmse[xgb_results$trained_on == "NEFIN (true)"]
    delta_rmse <- fia_rmse - nefin_rmse
    pct_improve <- 100 * delta_rmse / fia_rmse
    
    cat("  1. HOLDOUT PREDICTION (XGBoost):\n")
    cat("     FIA-trained RMSE:   ", round(fia_rmse, 1), "Mg/ha\n")
    cat("     NEFIN-trained RMSE: ", round(nefin_rmse, 1), "Mg/ha\n")
    cat("     Improvement:        ", round(delta_rmse, 1), "Mg/ha (", round(pct_improve, 1), "%)\n\n")
    
    summary_lines <- c(summary_lines, 
                       paste0("Holdout RMSE improvement: ", round(pct_improve, 1), "%"))
  }
}

# Attenuation
if (!is.null(compiled$attenuation)) {
  atten <- compiled$attenuation %>% filter(metric == "Attenuation %")
  if (nrow(atten) > 0) {
    cat("  2. RELATIONSHIP ATTENUATION:\n")
    cat("     NDVI-biomass slope attenuated by:", round(atten$value, 1), "%\n\n")
    
    summary_lines <- c(summary_lines,
                       paste0("NDVI slope attenuation: ", round(atten$value, 1), "%"))
  }
}

# Sensor comparison
if (!is.null(compiled$sensor)) {
  modis_sd <- compiled$sensor %>% filter(sensor == "MODIS", metric == "NDVI_SD")
  s2_sd <- compiled$sensor %>% filter(sensor == "Sentinel-2", metric == "NDVI_SD")
  
  if (nrow(modis_sd) > 0 && nrow(s2_sd) > 0) {
    ratio <- s2_sd$value / modis_sd$value
    cat("  3. SENSOR RESOLUTION EFFECT:\n")
    cat("     MODIS NDVI uncertainty:  ", round(modis_sd$value, 4), "\n")
    cat("     S2 NDVI uncertainty:     ", round(s2_sd$value, 4), "\n")
    cat("     S2/MODIS ratio:          ", round(ratio, 1), "x\n\n")
    
    summary_lines <- c(summary_lines,
                       paste0("S2 ", round(ratio, 1), "x more sensitive than MODIS"))
  }
}

# Heterogeneity effect
if (!is.null(compiled$heterogeneity)) {
  low_het <- compiled$heterogeneity %>% filter(het_class == "Low")
  high_het <- compiled$heterogeneity %>% filter(het_class == "High")
  
  if (nrow(low_het) > 0 && nrow(high_het) > 0) {
    cat("  4. LANDSCAPE HETEROGENEITY EFFECT:\n")
    cat("     Low heterogeneity ΔRMSE:  +", round(low_het$delta_rmse, 1), "Mg/ha\n")
    cat("     High heterogeneity ΔRMSE: +", round(high_het$delta_rmse, 1), "Mg/ha\n\n")
    
    summary_lines <- c(summary_lines,
                       paste0("High heterogeneity: +", round(high_het$delta_rmse, 1), " vs +",
                              round(low_het$delta_rmse, 1), " Mg/ha ΔRMSE"))
  }
}

# Hex summary
if (!is.null(compiled$hex_rmse)) {
  n_nefin_better <- sum(compiled$hex_rmse$delta_rmse > 0, na.rm = TRUE)
  n_total <- sum(!is.na(compiled$hex_rmse$delta_rmse))
  pct_better <- 100 * n_nefin_better / n_total
  
  cat("  5. SPATIAL PATTERN:\n")
  cat("     Hexes where NEFIN better:", n_nefin_better, "/", n_total, 
      "(", round(pct_better, 1), "%)\n")
  cat("     Mean ΔRMSE:", round(mean(compiled$hex_rmse$delta_rmse, na.rm = TRUE), 1), "Mg/ha\n\n")
  
  summary_lines <- c(summary_lines,
                     paste0("NEFIN better in ", round(pct_better, 1), "% of hexes"))
}

# =============================================================================
# CREATE PUBLICATION FIGURE
# =============================================================================

cat("═══════════════════════════════════════════════════════════════════════════════\n")
cat("CREATING PUBLICATION FIGURES\n")
cat("═══════════════════════════════════════════════════════════════════════════════\n\n")

# Figure: Combined summary
if (!is.null(compiled$holdout) && nrow(compiled$holdout) > 0) {
  
  # Panel A: Holdout RMSE
  fig_a <- ggplot(compiled$holdout, aes(x = model, y = holdout_rmse, fill = trained_on)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f", holdout_rmse)), 
              position = position_dodge(0.9), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("FIA (fuzzed)" = "#e74c3c", "NEFIN (true)" = "#27ae60")) +
    labs(title = "A) Holdout RMSE", 
         subtitle = "Lower is better",
         x = "", y = "RMSE (Mg/ha)", fill = "") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  
  # Panel B: Sensor resolution (if available)
  if (!is.null(compiled$sensor)) {
    sensor_sd <- compiled$sensor %>% filter(metric == "NDVI_SD")
    
    fig_b <- ggplot(sensor_sd, aes(x = sensor, y = value, fill = sensor)) +
      geom_col(alpha = 0.8) +
      geom_text(aes(label = sprintf("%.4f", value)), vjust = -0.5, size = 3.5) +
      scale_fill_manual(values = c("MODIS" = "#3498db", "Sentinel-2" = "#e74c3c")) +
      labs(title = "B) NDVI Uncertainty by Resolution",
           subtitle = "SD from ±1.6km fuzzing",
           x = "", y = "NDVI SD") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "none")
  } else {
    fig_b <- ggplot() + theme_void() + 
      labs(title = "B) Sensor comparison not available")
  }
  
  # Panel C: Heterogeneity effect (if available)
  if (!is.null(compiled$heterogeneity)) {
    fig_c <- ggplot(compiled$heterogeneity, aes(x = het_class, y = delta_rmse, fill = het_class)) +
      geom_col(alpha = 0.8) +
      geom_text(aes(label = sprintf("+%.1f", delta_rmse)), vjust = -0.5, size = 3.5) +
      scale_fill_viridis_d(option = "D") +
      labs(title = "C) ΔRMSE by Heterogeneity",
           subtitle = "Higher = fuzzing hurts more",
           x = "Landscape Heterogeneity", y = "ΔRMSE (Mg/ha)") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "none")
  } else {
    fig_c <- ggplot() + theme_void() +
      labs(title = "C) Heterogeneity analysis not available")
  }
  
  # Panel D: Hex distribution
  if (!is.null(compiled$hex_rmse)) {
    fig_d <- ggplot(compiled$hex_rmse, aes(x = delta_rmse)) +
      geom_histogram(bins = 30, fill = "#9b59b6", alpha = 0.7) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
      geom_vline(xintercept = mean(compiled$hex_rmse$delta_rmse, na.rm = TRUE),
                 linetype = "solid", color = "darkgreen", linewidth = 1) +
      labs(title = "D) Hex-Level ΔRMSE Distribution",
           subtitle = "Positive = NEFIN better",
           x = "ΔRMSE (Mg/ha)", y = "Count") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))
  } else {
    fig_d <- ggplot() + theme_void() +
      labs(title = "D) Hex analysis not available")
  }
  
  # Combine
  fig_combined <- (fig_a | fig_b) / (fig_c | fig_d) +
    plot_annotation(
      title = "Impact of Coordinate Fuzzing on Forest Biomass Prediction",
      subtitle = "NEFIN true coordinates vs FIA fuzzed coordinates (±1.6km)",
      theme = theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12)
      )
    )
  
  ggsave(fs::path(output_dir, "figures", "fig_summary_4panel.png"), fig_combined,
         width = 14, height = 10, dpi = 300)
  cat("  ✓ fig_summary_4panel.png\n")
}

# =============================================================================
# SAVE MASTER SUMMARY
# =============================================================================

# Create summary dataframe
if (length(summary_lines) > 0) {
  summary_df <- data.frame(
    finding = summary_lines
  )
  write_csv(summary_df, fs::path(output_dir, "key_findings.csv"))
  cat("  ✓ key_findings.csv\n")
}

# Save compiled results
saveRDS(compiled, fs::path(output_dir, "compiled_results.rds"))
cat("  ✓ compiled_results.rds\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
cat("║  PHASE 2 ANALYSIS COMPLETE                                                   ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("Output directory:", output_dir, "\n\n")

cat("KEY FINDINGS:\n")
cat("─────────────────────────────────────────────────────────────────────────────────\n")
for (line in summary_lines) {
  cat("  •", line, "\n")
}
cat("─────────────────────────────────────────────────────────────────────────────────\n")

cat("\nAnalysis status:\n")
for (name in names(results)) {
  status <- if (isTRUE(results[[name]])) "✓" else "✗"
  cat("  ", status, name, "\n")
}

cat("\n")
