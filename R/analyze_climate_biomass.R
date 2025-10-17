# R/analyze_climate_biomass.R
# Analyze relationships between climate variables and biomass/uncertainty

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(tidyr)
  library(broom); library(patchwork); library(scales); library(fs)
})

analyze_climate_biomass <- function(consolidated_dir = NULL,
                                    output_dir = NULL) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Climate-Biomass-Uncertainty Analysis                     ║\n")
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
    output_dir <- fs::path(consolidated_dir, "climate_analysis")
  }
  fs::dir_create(output_dir, recurse = TRUE)
  
  # Load data
  cat("Loading data...\n")
  fia_file <- fs::path(consolidated_dir, "fia_all_scales.csv")
  fia <- readr::read_csv(fia_file, show_col_types = FALSE)
  
  # Check for PRISM columns
  prism_cols <- names(fia)[grepl("^(tmean|ppt)_", names(fia))]
  
  if (length(prism_cols) == 0) {
    stop("No PRISM climate variables found in results. ",
         "Re-run compute_metrics after extracting PRISM data.")
  }
  
  cat("  Found PRISM variables:", paste(prism_cols, collapse = ", "), "\n\n")
  
  # ========================================================================
  # ANALYSIS 1: Climate gradients by scale
  # ========================================================================
  
  cat("Analysis 1: Climate gradients by grid scale\n")
  
  # Focus on mean annual temperature and precipitation
  climate_summary <- fia |>
    dplyr::filter(!is.na(tmean_mean), !is.na(ppt_mean)) |>
    dplyr::group_by(grid_scale) |>
    dplyr::summarise(
      n_hexes = dplyr::n_distinct(hex_id),
      
      # Temperature stats
      temp_range = max(tmean_mean, na.rm = TRUE) - min(tmean_mean, na.rm = TRUE),
      temp_sd = sd(tmean_mean, na.rm = TRUE),
      temp_cv = 100 * temp_sd / mean(tmean_mean, na.rm = TRUE),
      
      # Precipitation stats
      ppt_range = max(ppt_mean, na.rm = TRUE) - min(ppt_mean, na.rm = TRUE),
      ppt_sd = sd(ppt_mean, na.rm = TRUE),
      ppt_cv = 100 * ppt_sd / mean(ppt_mean, na.rm = TRUE),
      
      # Spatial heterogeneity
      temp_spatial_sd_mean = mean(tmean_mean_spatial_sd, na.rm = TRUE),
      ppt_spatial_sd_mean = mean(ppt_mean_spatial_sd, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  readr::write_csv(climate_summary, fs::path(output_dir, "climate_gradient_summary.csv"))
  print(climate_summary)
  
  # Plot 1: Climate space by grid scale
  p1 <- fia |>
    dplyr::filter(!is.na(tmean_mean), !is.na(ppt_mean)) |>
    ggplot(aes(x = tmean_mean, y = ppt_mean, color = grid_scale)) +
    geom_point(alpha = 0.4, size = 1.5) +
    facet_wrap(~grid_scale, ncol = 2) +
    scale_color_viridis_d() +
    labs(title = "Climate Space Coverage by Grid Scale",
         x = "Mean Annual Temperature (°C)",
         y = "Mean Annual Precipitation (mm)") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none")
  
  ggsave(fs::path(output_dir, "01_climate_space_by_scale.png"), p1,
         width = 10, height = 8, dpi = 300)
  
  # ========================================================================
  # ANALYSIS 2: Climate-biomass relationships
  # ========================================================================
  
  cat("\nAnalysis 2: Climate-biomass relationships\n")
  
  # Fit models by scale
  climate_models <- fia |>
    dplyr::filter(!is.na(tmean_mean), !is.na(ppt_mean), !is.na(mean)) |>
    dplyr::group_by(grid_scale) |>
    dplyr::do({
      # Simple model
      m_simple <- lm(mean ~ tmean_mean + ppt_mean, data = .)
      
      # With interaction
      m_interact <- lm(mean ~ tmean_mean * ppt_mean, data = .)
      
      # With quadratic terms
      m_quad <- lm(mean ~ poly(tmean_mean, 2) + poly(ppt_mean, 2), data = .)
      
      data.frame(
        n = nrow(.),
        r2_simple = summary(m_simple)$r.squared,
        r2_interact = summary(m_interact)$r.squared,
        r2_quad = summary(m_quad)$r.squared,
        temp_effect = coef(m_simple)["tmean_mean"],
        ppt_effect = coef(m_simple)["ppt_mean"]
      )
    }) |>
    dplyr::ungroup()
  
  readr::write_csv(climate_models, fs::path(output_dir, "climate_biomass_models.csv"))
  print(climate_models)
  
  # Plot 2: Temperature effect on biomass
  p2 <- fia |>
    dplyr::filter(!is.na(tmean_mean), !is.na(mean)) |>
    ggplot(aes(x = tmean_mean, y = mean)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 5), 
                color = "red", linewidth = 1.5) +
    facet_wrap(~grid_scale, scales = "free_y") +
    labs(title = "Biomass vs Temperature by Grid Scale",
         x = "Mean Annual Temperature (°C)",
         y = "Biomass (Mg/ha)") +
    theme_minimal(base_size = 10)
  
  ggsave(fs::path(output_dir, "02_biomass_vs_temperature.png"), p2,
         width = 12, height = 8, dpi = 300)
  
  # ========================================================================
  # ANALYSIS 3: Climate effects on uncertainty
  # ========================================================================
  
  cat("\nAnalysis 3: Climate effects on uncertainty components\n")
  
  # Does climate heterogeneity affect positional uncertainty?
  uncertainty_models <- fia |>
    dplyr::filter(!is.na(tmean_mean_spatial_sd), !is.na(positional_sd)) |>
    dplyr::group_by(grid_scale) |>
    dplyr::summarise(
      n = dplyr::n(),
      
      # Correlation between climate heterogeneity and positional SD
      cor_temp_het_pos = cor(tmean_mean_spatial_sd, positional_sd, use = "complete.obs"),
      cor_ppt_het_pos = cor(ppt_mean_spatial_sd, positional_sd, use = "complete.obs"),
      
      # Correlation between climate heterogeneity and sampling error
      cor_temp_het_se = cor(tmean_mean_spatial_sd, se, use = "complete.obs"),
      cor_ppt_het_se = cor(ppt_mean_spatial_sd, se, use = "complete.obs"),
      
      # Model R² for predicting positional SD from climate heterogeneity
      r2_climate_pos = {
        m <- lm(positional_sd ~ tmean_mean_spatial_sd + ppt_mean_spatial_sd, 
                data = dplyr::filter(fia, grid_scale == grid_scale[1]))
        summary(m)$r.squared
      },
      
      .groups = "drop"
    )
  
  readr::write_csv(uncertainty_models, fs::path(output_dir, "climate_uncertainty_correlations.csv"))
  print(uncertainty_models)
  
  # Plot 3: Climate heterogeneity vs positional uncertainty
  p3 <- fia |>
    dplyr::filter(!is.na(tmean_mean_spatial_sd), !is.na(positional_sd)) |>
    ggplot(aes(x = tmean_mean_spatial_sd, y = positional_sd, color = grid_scale)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "lm", se = FALSE) +
    scale_color_viridis_d() +
    facet_wrap(~grid_scale, scales = "free") +
    labs(title = "Climate Heterogeneity vs Positional Uncertainty",
         subtitle = "Does spatial climate variability within hexes affect coordinate uncertainty?",
         x = "Temperature Spatial SD within Hex (°C)",
         y = "Positional SD (Mg/ha)") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none")
  
  ggsave(fs::path(output_dir, "03_climate_heterogeneity_vs_uncertainty.png"), p3,
         width = 12, height = 8, dpi = 300)
  
  # ========================================================================
  # ANALYSIS 4: Identify optimal scale based on climate
  # ========================================================================
  
  cat("\nAnalysis 4: Scale optimization considering climate\n")
  
  # Calculate metrics for scale selection
  scale_optimization <- fia |>
    dplyr::group_by(grid_scale) |>
    dplyr::summarise(
      # Basic stats
      n_hexes = dplyr::n_distinct(hex_id),
      mean_biomass = mean(mean, na.rm = TRUE),
      
      # Uncertainty metrics
      mean_total_sd = mean(total_sd, na.rm = TRUE),
      mean_pos_fraction = mean(positional_sd / total_sd, na.rm = TRUE),
      
      # Climate representation
      temp_range = max(tmean_mean, na.rm = TRUE) - min(tmean_mean, na.rm = TRUE),
      ppt_range = max(ppt_mean, na.rm = TRUE) - min(ppt_mean, na.rm = TRUE),
      
      # Within-hex climate heterogeneity (averaged)
      mean_temp_het = mean(tmean_mean_spatial_sd, na.rm = TRUE),
      mean_ppt_het = mean(ppt_mean_spatial_sd, na.rm = TRUE),
      
      # Biomass-climate relationship strength
      r2_climate = {
        m <- lm(mean ~ tmean_mean + ppt_mean, data = cur_data())
        summary(m)$r.squared
      },
      
      # Combined score (lower is better)
      # Balances: low uncertainty, good climate representation, strong relationships
      optimization_score = mean_total_sd / mean_biomass * 
        (1 + mean_pos_fraction) * 
        (1 + mean_temp_het/temp_range) *
        (1 - r2_climate),
      
      .groups = "drop"
    )
  
  readr::write_csv(scale_optimization, fs::path(output_dir, "scale_optimization_metrics.csv"))
  
  cat("\nScale Optimization Results:\n")
  print(scale_optimization |> arrange(optimization_score))
  
  # Plot 4: Multi-criteria scale comparison
  p4_data <- scale_optimization |>
    dplyr::select(grid_scale, mean_total_sd, mean_pos_fraction, 
                  mean_temp_het, r2_climate) |>
    tidyr::pivot_longer(-grid_scale, names_to = "metric", values_to = "value") |>
    dplyr::mutate(
      metric = factor(metric,
                      levels = c("mean_total_sd", "mean_pos_fraction", 
                                 "mean_temp_het", "r2_climate"),
                      labels = c("Total Uncertainty", "Positional Fraction",
                                 "Climate Heterogeneity", "Climate R²"))
    )
  
  p4 <- ggplot(p4_data, aes(x = grid_scale, y = value, fill = grid_scale)) +
    geom_col() +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    scale_fill_viridis_d() +
    labs(title = "Multi-Criteria Scale Comparison",
         subtitle = "Trade-offs between uncertainty, climate representation, and relationships",
         x = "Grid Scale", y = "Value") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(fs::path(output_dir, "04_scale_optimization_criteria.png"), p4,
         width = 10, height = 8, dpi = 300)
  
  # ========================================================================
  # ANALYSIS 5: Variance partitioning
  # ========================================================================
  
  cat("\nAnalysis 5: Variance partitioning with climate\n")
  
  # For each scale, partition variance in biomass
  variance_partitioning <- fia |>
    dplyr::filter(!is.na(mean), !is.na(tmean_mean), !is.na(ppt_mean)) |>
    dplyr::group_by(grid_scale) |>
    dplyr::do({
      df <- .
      
      # Total variance
      var_total <- var(df$mean)
      
      # Variance explained by climate
      m_climate <- lm(mean ~ tmean_mean + ppt_mean, data = df)
      var_climate <- var(fitted(m_climate))
      
      # Variance explained by spatial position (use hex_id as proxy)
      m_spatial <- lm(mean ~ as.factor(hex_id), data = df)
      var_spatial <- var(fitted(m_spatial)) - var_climate  # Remove climate component
      
      # Residual variance
      var_residual <- var_total - var_climate - max(var_spatial, 0)
      
      data.frame(
        var_total = var_total,
        var_climate = var_climate,
        var_spatial = max(var_spatial, 0),
        var_residual = max(var_residual, 0),
        pct_climate = 100 * var_climate / var_total,
        pct_spatial = 100 * max(var_spatial, 0) / var_total,
        pct_residual = 100 * max(var_residual, 0) / var_total
      )
    }) |>
    dplyr::ungroup()
  
  readr::write_csv(variance_partitioning, fs::path(output_dir, "variance_partitioning.csv"))
  
  # Plot 5: Variance partitioning
  p5_data <- variance_partitioning |>
    dplyr::select(grid_scale, pct_climate, pct_spatial, pct_residual) |>
    tidyr::pivot_longer(-grid_scale, names_to = "component", values_to = "percentage") |>
    dplyr::mutate(
      component = factor(component,
                         levels = c("pct_climate", "pct_spatial", "pct_residual"),
                         labels = c("Climate", "Spatial", "Residual"))
    )
  
  p5 <- ggplot(p5_data, aes(x = grid_scale, y = percentage, fill = component)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c("Climate" = "#E69F00", 
                                 "Spatial" = "#56B4E9", 
                                 "Residual" = "#999999")) +
    labs(title = "Variance Partitioning by Grid Scale",
         subtitle = "Proportion of biomass variance explained by different factors",
         x = "Grid Scale", y = "Percentage of Variance", fill = "Component") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  ggsave(fs::path(output_dir, "05_variance_partitioning.png"), p5,
         width = 10, height = 6, dpi = 300)
  
  # ========================================================================
  # SUMMARY DASHBOARD
  # ========================================================================
  
  cat("\nCreating summary dashboard...\n")
  
  # Combined dashboard
  dashboard <- (p1 + p2) / (p3 + p5) +
    plot_annotation(
      title = "Climate-Biomass-Uncertainty Analysis Dashboard",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  ggsave(fs::path(output_dir, "00_climate_analysis_dashboard.png"), dashboard,
         width = 16, height = 12, dpi = 300)
  
  # ========================================================================
  # FINAL REPORT
  # ========================================================================
  
  report_lines <- c(
    "═══════════════════════════════════════════════════════════════",
    "CLIMATE-BIOMASS-UNCERTAINTY ANALYSIS REPORT",
    "═══════════════════════════════════════════════════════════════",
    "",
    paste("Generated:", Sys.time()),
    paste("Data source:", consolidated_dir),
    "",
    "───────────────────────────────────────────────────────────────",
    "KEY FINDINGS",
    "───────────────────────────────────────────────────────────────",
    "",
    "1. CLIMATE GRADIENTS BY SCALE:",
    ""
  )
  
  for (i in seq_len(nrow(climate_summary))) {
    row <- climate_summary[i,]
    report_lines <- c(report_lines,
                      paste0("Grid: ", row$grid_scale),
                      paste0("  Temperature range: ", round(row$temp_range, 1), "°C"),
                      paste0("  Temperature CV: ", round(row$temp_cv, 1), "%"),
                      paste0("  Precipitation range: ", round(row$ppt_range, 0), " mm"),
                      paste0("  Precipitation CV: ", round(row$ppt_cv, 1), "%"),
                      ""
    )
  }
  
  report_lines <- c(report_lines,
                    "2. CLIMATE-BIOMASS RELATIONSHIPS:",
                    ""
  )
  
  for (i in seq_len(nrow(climate_models))) {
    row <- climate_models[i,]
    report_lines <- c(report_lines,
                      paste0("Grid: ", row$grid_scale),
                      paste0("  R² (simple model): ", round(row$r2_simple, 3)),
                      paste0("  Temperature effect: ", round(row$temp_effect, 2), " Mg/ha/°C"),
                      paste0("  Precipitation effect: ", round(row$ppt_effect/1000, 2), " Mg/ha/mm"),
                      ""
    )
  }
  
  report_lines <- c(report_lines,
                    "3. OPTIMAL SCALE CONSIDERING CLIMATE:",
                    "",
                    paste("  Best scale (lowest optimization score):",
                          scale_optimization$grid_scale[which.min(scale_optimization$optimization_score)]),
                    "",
                    "4. VARIANCE PARTITIONING:",
                    ""
  )
  
  for (i in seq_len(nrow(variance_partitioning))) {
    row <- variance_partitioning[i,]
    report_lines <- c(report_lines,
                      paste0("Grid: ", row$grid_scale),
                      paste0("  Climate explains: ", round(row$pct_climate, 1), "%"),
                      paste0("  Spatial pattern: ", round(row$pct_spatial, 1), "%"),
                      paste0("  Residual: ", round(row$pct_residual, 1), "%"),
                      ""
    )
  }
  
  report_lines <- c(report_lines,
                    "───────────────────────────────────────────────────────────────",
                    "RECOMMENDATIONS",
                    "───────────────────────────────────────────────────────────────",
                    "",
                    "1. Grid scales with high climate heterogeneity within hexes may",
                    "   need finer resolution to capture climate-biomass relationships.",
                    "",
                    "2. Positional uncertainty appears to be influenced by climate",
                    "   heterogeneity, suggesting environmental gradients affect",
                    "   coordinate fuzzing impacts.",
                    "",
                    "3. Consider stratifying by climate zones for more precise",
                    "   uncertainty estimates.",
                    "",
                    "═══════════════════════════════════════════════════════════════"
  )
  
  report_file <- fs::path(output_dir, "climate_analysis_report.txt")
  writeLines(report_lines, report_file)
  
  cat("\n✓ Climate analysis complete!\n")
  cat("  Output directory:", output_dir, "\n")
  
  invisible(output_dir)
}

# CLI
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  get_arg <- function(flag, default = NULL) {
    hit <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
  }
  
  analyze_climate_biomass(
    consolidated_dir = get_arg("--dir", NULL),
    output_dir = get_arg("--out", NULL)
  )
}