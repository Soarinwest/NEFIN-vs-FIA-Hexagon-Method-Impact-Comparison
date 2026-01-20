#!/usr/bin/env Rscript
# =============================================================================
# analyze_prediction_uncertainty.R
# 
# Analyze WHERE coordinate fuzzing impacts predictions most
# Connect uncertainty to landscape characteristics
#
# Author: Soren Donisvitch
# Date: December 2024
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(fs)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

config <- list(
  # Prediction rasters
  raster_dir = "runs/fuzzing_prediction_comparison/rasters",
  
  # Hex grids for aggregation
  hex_grids = list(
    `10kha` = "data/hex/hex_grid_10kha.geojson",
    `50kha` = "data/hex/hex_grid_50kha.geojson"
  ),
  
  # State boundaries
  states_path = "data/boundaries/states_5070.geojson",
  
  # Source NDVI for heterogeneity calc
  ndvi_path = "data/processed/ndvi/modis/MODIS_NDVI_5yr_blocked_2020_2024.tif",
  
  # Output
  output_dir = "runs/fuzzing_prediction_comparison/analysis"
)

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

analyze_uncertainty_patterns <- function(cfg = config) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  Analyzing Prediction Uncertainty Patterns                                   ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  fs::dir_create(cfg$output_dir, recurse = TRUE)
  fs::dir_create(fs::path(cfg$output_dir, "figures"), recurse = TRUE)
  
  # ===========================================================================
  # 1. LOAD RASTERS
  # ===========================================================================
  
  cat("Loading prediction rasters...\n")
  
  pred_nefin <- rast(fs::path(cfg$raster_dir, "biomass_pred_nefin.tif"))
  pred_fia <- rast(fs::path(cfg$raster_dir, "biomass_pred_fia.tif"))
  mc_sd <- rast(fs::path(cfg$raster_dir, "biomass_pred_mc_sd.tif"))
  mc_mean <- rast(fs::path(cfg$raster_dir, "biomass_pred_mc_mean.tif"))
  diff_map <- rast(fs::path(cfg$raster_dir, "diff_nefin_minus_fia.tif"))
  cv_map <- rast(fs::path(cfg$raster_dir, "prediction_cv_percent.tif"))
  
  cat("  ✓ Loaded 6 rasters\n")
  
  # ===========================================================================
  # 2. COMPUTE NDVI HETEROGENEITY (local SD)
  # ===========================================================================
  
  cat("\nComputing NDVI heterogeneity (local standard deviation)...\n")
  
  if (fs::file_exists(cfg$ndvi_path)) {
    ndvi <- rast(cfg$ndvi_path)
    
    # Compute local SD using a 5x5 window (~1.25km at 250m resolution)
    ndvi_sd <- focal(ndvi, w = 5, fun = "sd", na.rm = TRUE)
    names(ndvi_sd) <- "ndvi_heterogeneity"
    
    # Resample to match prediction rasters if needed
    if (!compareGeom(ndvi_sd, mc_sd, stopOnError = FALSE)) {
      ndvi_sd <- resample(ndvi_sd, mc_sd, method = "bilinear")
    }
    
    cat("  ✓ Computed NDVI heterogeneity\n")
    
    # Save
    writeRaster(ndvi_sd, fs::path(cfg$output_dir, "ndvi_heterogeneity.tif"), 
                overwrite = TRUE)
  } else {
    cat("  ⚠ NDVI raster not found, skipping heterogeneity\n")
    ndvi_sd <- NULL
  }
  
  # ===========================================================================
  # 3. CORRELATION: UNCERTAINTY vs HETEROGENEITY
  # ===========================================================================
  
  cat("\nAnalyzing uncertainty drivers...\n")
  
  # Sample points for correlation analysis
  set.seed(42)
  n_sample <- 50000
  
  # Stack relevant layers
  if (!is.null(ndvi_sd)) {
    analysis_stack <- c(mc_sd, pred_nefin, diff_map, cv_map, ndvi_sd)
    names(analysis_stack) <- c("mc_sd", "biomass", "diff", "cv", "ndvi_het")
  } else {
    analysis_stack <- c(mc_sd, pred_nefin, diff_map, cv_map)
    names(analysis_stack) <- c("mc_sd", "biomass", "diff", "cv")
  }
  
  # Sample random points
  sample_cells <- sample(which(!is.na(values(mc_sd))), min(n_sample, sum(!is.na(values(mc_sd)))))
  sample_vals <- analysis_stack[sample_cells]
  sample_df <- as.data.frame(sample_vals)
  sample_df <- sample_df[complete.cases(sample_df), ]
  
  cat("  Sampled", nrow(sample_df), "pixels for analysis\n")
  
  # Correlations
  cat("\n  Correlations with MC uncertainty (SD):\n")
  cor_biomass <- cor(sample_df$mc_sd, sample_df$biomass, use = "complete.obs")
  cat("    Biomass:           r =", round(cor_biomass, 3), "\n")
  
  if ("ndvi_het" %in% names(sample_df)) {
    cor_het <- cor(sample_df$mc_sd, sample_df$ndvi_het, use = "complete.obs")
    cat("    NDVI heterogeneity: r =", round(cor_het, 3), "\n")
  }
  
  # ===========================================================================
  # 4. CREATE UNCERTAINTY HOTSPOT MAP
  # ===========================================================================
  
  cat("\nIdentifying uncertainty hotspots...\n")
  
  # Define hotspot as MC_SD > 75th percentile
  sd_75 <- global(mc_sd, fun = function(x) quantile(x, 0.75, na.rm = TRUE))[[1]]
  sd_90 <- global(mc_sd, fun = function(x) quantile(x, 0.90, na.rm = TRUE))[[1]]
  
  cat("  75th percentile SD:", round(sd_75, 1), "Mg/ha\n")
  cat("  90th percentile SD:", round(sd_90, 1), "Mg/ha\n")
  
  # Create hotspot classification
  hotspot_class <- classify(mc_sd, 
    rcl = matrix(c(
      -Inf, sd_75, 1,   # Low uncertainty
      sd_75, sd_90, 2,  # Moderate uncertainty
      sd_90, Inf, 3     # High uncertainty (hotspot)
    ), ncol = 3, byrow = TRUE)
  )
  names(hotspot_class) <- "uncertainty_class"
  
  writeRaster(hotspot_class, fs::path(cfg$output_dir, "uncertainty_hotspots.tif"),
              overwrite = TRUE)
  cat("  ✓ Saved uncertainty_hotspots.tif\n")
  
  # ===========================================================================
  # 5. AGGREGATE TO HEX SCALES
  # ===========================================================================
  
  cat("\nAggregating uncertainty to hex scales...\n")
  
  hex_results <- list()
  
  for (scale_name in names(cfg$hex_grids)) {
    hex_path <- cfg$hex_grids[[scale_name]]
    
    if (!fs::file_exists(hex_path)) {
      cat("  ⚠", scale_name, "hex grid not found\n")
      next
    }
    
    cat("  Processing", scale_name, "...\n")
    
    hex <- st_read(hex_path, quiet = TRUE)
    
    # Transform to match raster CRS
    hex <- st_transform(hex, crs(mc_sd))
    
    # Extract zonal statistics
    # Ensure unique names
    extract_stack <- c(mc_sd, pred_nefin, pred_fia, diff_map)
    names(extract_stack) <- c("mc_sd", "pred_nefin", "pred_fia", "diff")
    
    hex_stats <- exact_extract(
      extract_stack,
      hex,
      fun = c("mean", "stdev", "median"),
      append_cols = "hex_id"
    )
    
    hex_results[[scale_name]] <- hex_stats
    
    # Save
    write_csv(hex_stats, fs::path(cfg$output_dir, paste0("hex_uncertainty_", scale_name, ".csv")))
    cat("    ✓ Saved hex_uncertainty_", scale_name, ".csv\n")
  }
  
  # ===========================================================================
  # 6. STATE-LEVEL SUMMARY
  # ===========================================================================
  
  cat("\nComputing state-level summaries...\n")
  
  if (fs::file_exists(cfg$states_path)) {
    states <- st_read(cfg$states_path, quiet = TRUE)
    states <- st_transform(states, crs(mc_sd))
    
    # Only NE states
    ne_states <- c("CT", "MA", "ME", "NH", "NJ", "NY", "PA", "RI", "VT")
    if ("STUSPS" %in% names(states)) {
      states <- states %>% filter(STUSPS %in% ne_states)
    }
    
    # Ensure unique names for extraction
    state_stack <- c(mc_sd, pred_nefin, pred_fia)
    names(state_stack) <- c("mc_sd", "pred_nefin", "pred_fia")
    
    state_stats <- exact_extract(
      state_stack,
      states,
      fun = c("mean", "stdev", "median"),
      append_cols = c("STUSPS", "NAME")
    )
    
    cat("\n  State-level uncertainty (MC SD):\n")
    cat("  ─────────────────────────────────────────\n")
    state_summary <- state_stats %>%
      arrange(desc(mean.mc_sd)) %>%
      select(STUSPS, mean_mc_sd = mean.mc_sd,
             mean_nefin = mean.pred_nefin,
             mean_fia = mean.pred_fia)
    
    print(state_summary)
    
    write_csv(state_stats, fs::path(cfg$output_dir, "state_uncertainty_summary.csv"))
    cat("\n  ✓ Saved state_uncertainty_summary.csv\n")
  }
  
  # ===========================================================================
  # 7. CREATE FIGURES
  # ===========================================================================
  
  cat("\nCreating figures...\n")
  
  # Figure 1: Uncertainty vs Biomass scatter
  p1 <- ggplot(sample_df, aes(x = biomass, y = mc_sd)) +
    geom_hex(bins = 50) +
    scale_fill_viridis_c(option = "magma", trans = "log10") +
    geom_smooth(method = "loess", color = "white", se = FALSE) +
    labs(
      title = "Prediction Uncertainty vs Biomass",
      subtitle = sprintf("r = %.3f | Higher biomass = more fuzzing impact", cor_biomass),
      x = "Predicted Biomass (Mg/ha)",
      y = "MC Prediction SD (Mg/ha)",
      fill = "Count"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(fs::path(cfg$output_dir, "figures", "uncertainty_vs_biomass.png"),
         p1, width = 10, height = 7, dpi = 300)
  cat("  ✓ uncertainty_vs_biomass.png\n")
  
  # Figure 2: Uncertainty vs NDVI heterogeneity
  if ("ndvi_het" %in% names(sample_df)) {
    p2 <- ggplot(sample_df, aes(x = ndvi_het, y = mc_sd)) +
      geom_hex(bins = 50) +
      scale_fill_viridis_c(option = "viridis", trans = "log10") +
      geom_smooth(method = "loess", color = "red", se = FALSE) +
      labs(
        title = "Prediction Uncertainty vs Landscape Heterogeneity",
        subtitle = sprintf("r = %.3f | Heterogeneous landscapes = more fuzzing impact", cor_het),
        x = "NDVI Local SD (heterogeneity)",
        y = "MC Prediction SD (Mg/ha)",
        fill = "Count"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(fs::path(cfg$output_dir, "figures", "uncertainty_vs_heterogeneity.png"),
           p2, width = 10, height = 7, dpi = 300)
    cat("  ✓ uncertainty_vs_heterogeneity.png\n")
  }
  
  # Figure 3: Histogram of differences
  p3 <- ggplot(sample_df, aes(x = diff)) +
    geom_histogram(bins = 100, fill = "#3498db", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    geom_vline(xintercept = mean(sample_df$diff, na.rm = TRUE), 
               linetype = "solid", color = "darkred", linewidth = 1) +
    annotate("text", x = mean(sample_df$diff) + 5, y = Inf, 
             label = sprintf("Mean = %.1f Mg/ha", mean(sample_df$diff, na.rm = TRUE)),
             vjust = 2, hjust = 0, color = "darkred") +
    labs(
      title = "NEFIN - FIA Prediction Difference",
      subtitle = "Positive = NEFIN predicts higher biomass",
      x = "Prediction Difference (Mg/ha)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(fs::path(cfg$output_dir, "figures", "prediction_difference_hist.png"),
         p3, width = 10, height = 6, dpi = 300)
  cat("  ✓ prediction_difference_hist.png\n")
  
  # ===========================================================================
  # 8. SUMMARY
  # ===========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  ANALYSIS COMPLETE                                                           ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  cat("Outputs:", cfg$output_dir, "\n\n")
  
  cat("KEY FINDINGS:\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  cat("  Uncertainty-Biomass correlation:     r =", round(cor_biomass, 3), "\n")
  if (exists("cor_het")) {
    cat("  Uncertainty-Heterogeneity correlation: r =", round(cor_het, 3), "\n")
  }
  cat("  High uncertainty threshold (90th):   ", round(sd_90, 1), "Mg/ha\n")
  cat("─────────────────────────────────────────────────────────────────────────────\n")
  cat("\n")
  
  invisible(list(
    sample_data = sample_df,
    correlations = list(biomass = cor_biomass)
  ))
}

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  # Check for exactextractr
  if (!requireNamespace("exactextractr", quietly = TRUE)) {
    cat("Installing exactextractr for zonal statistics...\n")
    install.packages("exactextractr")
  }
  library(exactextractr)
  
  analyze_uncertainty_patterns()
}
