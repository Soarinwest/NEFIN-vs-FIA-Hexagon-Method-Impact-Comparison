#!/usr/bin/env Rscript
# R/create_prediction_maps.R
# Create MODIS-resolution biomass prediction maps and hex-level ΔRMSE maps
# Author: Soren Donisvitch
# Date: December 2024

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sf)
  library(terra)
  library(xgboost)
  library(patchwork)
  library(viridis)
  library(scales)
  library(fs)
})

# =============================================================================
# CREATE MODIS-RESOLUTION PREDICTION SURFACE
# =============================================================================

create_modis_prediction_map <- function(
    model,
    ndvi_path,
    tmean_path,
    ppt_path,
    states_path = "data/boundaries/states_5070.geojson",
    output_path,
    title = "Predicted Biomass",
    chunk_size = 500000  # Process in chunks to avoid memory issues
) {
  
  cat("  Loading rasters...\n")
  
  # Load rasters
  r_ndvi <- terra::rast(ndvi_path)
  r_tmean <- terra::rast(tmean_path)
  r_ppt <- terra::rast(ppt_path)
  
  cat("    NDVI res:", res(r_ndvi), "| dims:", dim(r_ndvi)[1:2], "\n")
  cat("    Climate res:", res(r_tmean), "\n")
  
  # Resample climate to MODIS grid (MODIS is finer resolution)
  cat("  Resampling climate to MODIS grid...\n")
  r_tmean_resamp <- terra::resample(r_tmean, r_ndvi, method = "bilinear")
  r_ppt_resamp <- terra::resample(r_ppt, r_ndvi, method = "bilinear")
  
  # Create output raster (same grid as NDVI)
  r_pred <- r_ndvi
  names(r_pred) <- "predicted_biomass"
  
  # Get total cells
  n_cells <- ncell(r_ndvi)
  cat("  Total pixels:", format(n_cells, big.mark = ","), "\n")
  
  # Process in chunks to manage memory
  cat("  Predicting biomass at each pixel...\n")
  
  n_chunks <- ceiling(n_cells / chunk_size)
  pb <- txtProgressBar(min = 0, max = n_chunks, style = 3)
  
  for (i in 1:n_chunks) {
    start_cell <- (i - 1) * chunk_size + 1
    end_cell <- min(i * chunk_size, n_cells)
    
    # Extract values for this chunk
    ndvi_vals <- values(r_ndvi)[start_cell:end_cell]
    tmean_vals <- values(r_tmean_resamp)[start_cell:end_cell]
    ppt_vals <- values(r_ppt_resamp)[start_cell:end_cell]
    
    # Find valid pixels (all three covariates present)
    valid <- !is.na(ndvi_vals) & !is.na(tmean_vals) & !is.na(ppt_vals)
    
    if (sum(valid) > 0) {
      # Build prediction matrix
      X_chunk <- cbind(
        ndvi_modis = ndvi_vals[valid],
        tmean = tmean_vals[valid],
        ppt = ppt_vals[valid]
      )
      
      # Predict
      preds <- predict(model, xgb.DMatrix(X_chunk))
      
      # Store predictions
      chunk_preds <- rep(NA, length(valid))
      chunk_preds[valid] <- preds
      
      # Write to raster
      values(r_pred)[start_cell:end_cell] <- chunk_preds
    }
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Count valid predictions
  valid_preds <- sum(!is.na(values(r_pred)))
  cat("  Valid predictions:", format(valid_preds, big.mark = ","), 
      "(", round(100 * valid_preds / n_cells, 1), "%)\n")
  
  # Save raster
  raster_path <- gsub("\\.png$", ".tif", output_path)
  terra::writeRaster(r_pred, raster_path, overwrite = TRUE)
  cat("  ✓ Saved GeoTIFF:", raster_path, "\n")
  
  # Create visualization
  cat("  Creating map visualization...\n")
  
  # Convert to data frame for plotting (subsample for speed)
  pred_df <- as.data.frame(r_pred, xy = TRUE, na.rm = TRUE)
  names(pred_df)[3] <- "biomass"
  
  set.seed(42)
  if (nrow(pred_df) > 300000) {
    plot_df <- pred_df %>% sample_n(300000)
  } else {
    plot_df <- pred_df
  }
  
  # Load states
  states <- if (fs::file_exists(states_path)) {
    st_read(states_path, quiet = TRUE) %>% st_transform(crs(r_pred))
  } else NULL
  
  # Summary stats for subtitle
  mean_bio <- round(mean(pred_df$biomass, na.rm = TRUE), 1)
  sd_bio <- round(sd(pred_df$biomass, na.rm = TRUE), 1)
  
  p <- ggplot() +
    geom_raster(data = plot_df, aes(x = x, y = y, fill = biomass)) +
    scale_fill_viridis_c(
      option = "G", 
      direction = -1,
      limits = c(0, 250),
      oob = scales::squish,
      na.value = "transparent"
    ) +
    labs(
      title = title,
      subtitle = paste0("MODIS resolution (~250m) | Mean: ", mean_bio, " ± ", sd_bio, " Mg/ha"),
      fill = "Biomass\n(Mg/ha)"
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    )
  
  if (!is.null(states)) {
    p <- p + geom_sf(data = states, fill = NA, color = "black", linewidth = 0.4)
  }
  
  ggsave(output_path, p, width = 12, height = 10, dpi = 300)
  cat("  ✓ Saved map:", output_path, "\n")
  
  return(r_pred)
}

# =============================================================================
# CREATE DIFFERENCE MAP (NEFIN - FIA PREDICTIONS)
# =============================================================================

create_prediction_difference_map <- function(
    r_pred_nefin,
    r_pred_fia,
    states_path = "data/boundaries/states_5070.geojson",
    output_path
) {
  
  cat("  Computing difference (NEFIN - FIA)...\n")
  
  r_diff <- r_pred_nefin - r_pred_fia
  names(r_diff) <- "difference"
  
  # Save raster
  raster_path <- gsub("\\.png$", ".tif", output_path)
  terra::writeRaster(r_diff, raster_path, overwrite = TRUE)
  cat("  ✓ Saved GeoTIFF:", raster_path, "\n")
  
  # Get statistics
  diff_vals <- values(r_diff)
  diff_vals <- diff_vals[!is.na(diff_vals)]
  
  cat("\n  Difference statistics (NEFIN - FIA predictions):\n")
  cat("    Mean difference:", round(mean(diff_vals), 2), "Mg/ha\n")
  cat("    Median difference:", round(median(diff_vals), 2), "Mg/ha\n")
  cat("    SD:", round(sd(diff_vals), 2), "Mg/ha\n")
  cat("    Range:", round(min(diff_vals), 1), "to", round(max(diff_vals), 1), "Mg/ha\n")
  cat("    Pixels where NEFIN > FIA:", round(100 * sum(diff_vals > 0) / length(diff_vals), 1), "%\n")
  cat("    Pixels where |diff| > 25 Mg/ha:", round(100 * sum(abs(diff_vals) > 25) / length(diff_vals), 1), "%\n")
  
  # Convert to data frame for plotting
  diff_df <- as.data.frame(r_diff, xy = TRUE, na.rm = TRUE)
  names(diff_df)[3] <- "difference"
  
  # Subsample for plotting
  set.seed(42)
  if (nrow(diff_df) > 300000) {
    plot_df <- diff_df %>% sample_n(300000)
  } else {
    plot_df <- diff_df
  }
  
  # Load states
  states <- if (fs::file_exists(states_path)) {
    st_read(states_path, quiet = TRUE) %>% st_transform(crs(r_diff))
  } else NULL
  
  # Create map
  p <- ggplot() +
    geom_raster(data = plot_df, aes(x = x, y = y, fill = difference)) +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0,
      limits = c(-75, 75),
      oob = scales::squish,
      na.value = "transparent"
    ) +
    labs(
      title = "Prediction Difference: NEFIN-trained − FIA-trained Model",
      subtitle = paste0("Mean diff: ", round(mean(diff_vals), 1), " Mg/ha | ",
                       "Red = NEFIN predicts higher | Blue = FIA predicts higher"),
      fill = "Δ Biomass\n(Mg/ha)"
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    )
  
  if (!is.null(states)) {
    p <- p + geom_sf(data = states, fill = NA, color = "black", linewidth = 0.4)
  }
  
  ggsave(output_path, p, width = 12, height = 10, dpi = 300)
  cat("  ✓ Saved map:", output_path, "\n")
  
  # Also create histogram of differences
  hist_path <- gsub("\\.png$", "_histogram.png", output_path)
  
  p_hist <- ggplot(data.frame(diff = diff_vals), aes(x = diff)) +
    geom_histogram(bins = 100, fill = "#3498db", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    geom_vline(xintercept = mean(diff_vals), linetype = "solid", color = "darkgreen", linewidth = 1) +
    annotate("text", x = mean(diff_vals), y = Inf, vjust = 2, hjust = -0.1,
             label = paste0("Mean: ", round(mean(diff_vals), 1)), color = "darkgreen") +
    labs(
      title = "Distribution of Prediction Differences (NEFIN - FIA)",
      subtitle = "At every MODIS pixel in the study region",
      x = "Biomass Difference (Mg/ha)",
      y = "Pixel Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(hist_path, p_hist, width = 10, height = 6, dpi = 300)
  cat("  ✓ Saved histogram:", hist_path, "\n")
  
  return(r_diff)
}

# =============================================================================
# CREATE HEX-LEVEL ΔRMSE MAP
# =============================================================================

create_hex_delta_rmse_map <- function(
    nefin_data,
    model_nefin,
    model_fia,
    hex_path = "data/hex/hex_grid_50kha.geojson",
    states_path = "data/boundaries/states_5070.geojson",
    output_path
) {
  
  cat("  Computing hex-level ΔRMSE...\n")
  
  # Load hex grid
  hex <- st_read(hex_path, quiet = TRUE)
  
  # Ensure hex_id column
  if (!"hex_id" %in% names(hex)) {
    hex$hex_id <- hex[[1]]
  }
  
  # Convert NEFIN to sf and spatial join with hex
  nefin_sf <- st_as_sf(nefin_data, coords = c("lon", "lat"), crs = 4326) %>%
    st_transform(st_crs(hex))
  
  nefin_with_hex <- st_join(nefin_sf, hex %>% select(hex_id))
  
  # Predict with both models
  X <- as.matrix(st_drop_geometry(nefin_with_hex) %>% select(ndvi_modis, tmean, ppt))
  
  nefin_with_hex$pred_nefin <- predict(model_nefin, xgb.DMatrix(X))
  nefin_with_hex$pred_fia <- predict(model_fia, xgb.DMatrix(X))
  nefin_with_hex$observed <- nefin_with_hex$aglb_Mg_per_ha
  
  # Compute hex-level RMSE
  hex_rmse <- nefin_with_hex %>%
    st_drop_geometry() %>%
    filter(!is.na(hex_id)) %>%
    group_by(hex_id) %>%
    summarise(
      n_plots = n(),
      rmse_nefin = sqrt(mean((observed - pred_nefin)^2)),
      rmse_fia = sqrt(mean((observed - pred_fia)^2)),
      delta_rmse = rmse_fia - rmse_nefin,  # Positive = NEFIN better
      mean_biomass = mean(observed),
      .groups = "drop"
    ) %>%
    filter(n_plots >= 3)  # Require at least 3 plots per hex
  
  cat("  Hexes with sufficient plots:", nrow(hex_rmse), "\n")
  
  # Join to hex geometry
  hex_with_rmse <- hex %>%
    left_join(hex_rmse, by = "hex_id")
  
  # Load states
  states <- if (fs::file_exists(states_path)) {
    st_read(states_path, quiet = TRUE)
  } else NULL
  
  # Create map
  p <- ggplot() +
    geom_sf(data = hex, fill = "gray95", color = "gray80", linewidth = 0.1) +
    geom_sf(data = hex_with_rmse %>% filter(!is.na(delta_rmse)), 
            aes(fill = delta_rmse), color = "gray50", linewidth = 0.1) +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0,
      limits = c(-30, 30),
      oob = scales::squish
    ) +
    labs(
      title = "Hex-Level ΔRMSE: FIA-trained − NEFIN-trained",
      subtitle = "Red = true coordinates improve RMSE | Blue = fuzzed coords better",
      fill = "ΔRMSE\n(Mg/ha)"
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "right"
    )
  
  if (!is.null(states)) {
    p <- p + geom_sf(data = states, fill = NA, color = "black", linewidth = 0.5)
  }
  
  ggsave(output_path, p, width = 12, height = 10, dpi = 300)
  cat("  ✓ Saved:", output_path, "\n")
  
  # Statistics
  cat("\n  ΔRMSE statistics:\n")
  cat("    Mean ΔRMSE:", round(mean(hex_rmse$delta_rmse, na.rm = TRUE), 2), "Mg/ha\n")
  cat("    Hexes where NEFIN better:", sum(hex_rmse$delta_rmse > 0, na.rm = TRUE), "\n")
  cat("    Hexes where FIA better:", sum(hex_rmse$delta_rmse < 0, na.rm = TRUE), "\n")
  
  # Save hex RMSE data
  csv_path <- gsub("\\.png$", ".csv", output_path)
  write_csv(hex_rmse, csv_path)
  
  return(hex_rmse)
}

# =============================================================================
# CREATE FUZZ PRESSURE VS ΔRMSE SCATTER
# =============================================================================

create_fuzz_pressure_plot <- function(
    hex_rmse,
    hex_path = "data/hex/hex_grid_50kha.geojson",
    output_path
) {
  
  cat("  Creating fuzz pressure vs ΔRMSE plot...\n")
  
  # Fuzz pressure could be related to:
  # - NDVI heterogeneity within hex
  # - Plot density
  # - Biomass variability
  
  # Use n_plots as proxy for sampling density (inverse of fuzz pressure impact)
  # More plots = less sensitive to individual fuzzing
  
  p <- ggplot(hex_rmse, aes(x = n_plots, y = delta_rmse)) +
    geom_point(aes(color = mean_biomass), alpha = 0.6, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_smooth(method = "loess", color = "#e74c3c", se = TRUE) +
    scale_color_viridis_c(option = "G", direction = -1) +
    labs(
      title = "Sample Density vs ΔRMSE",
      subtitle = "Does true coordinates help more in sparsely sampled hexes?",
      x = "Number of Plots per Hex",
      y = "ΔRMSE (FIA - NEFIN, Mg/ha)",
      color = "Mean Biomass\n(Mg/ha)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(output_path, p, width = 10, height = 7, dpi = 300)
  cat("  ✓ Saved:", output_path, "\n")
}

# =============================================================================
# MAIN FUNCTION
# =============================================================================

create_all_prediction_maps <- function(
    model_dir = "runs/spatial_model_comparison",
    output_dir = "runs/prediction_maps",
    window = "2020_2024"
) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  Creating Prediction Maps                                                 ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  fs::dir_create(output_dir, recurse = TRUE)
  
  # Paths
  ndvi_path <- paste0("data/processed/ndvi/modis/MODIS_NDVI_5yr_blocked_", window, ".tif")
  tmean_path <- paste0("data/processed/prism/prism_tmean_ne_", window, ".tif")
  ppt_path <- paste0("data/processed/prism/prism_ppt_ne_", window, ".tif")
  
  # Check files exist
  if (!fs::file_exists(ndvi_path)) {
    stop("NDVI raster not found: ", ndvi_path)
  }
  
  # ===========================================================================
  # 1. TRAIN MODELS (if not already done)
  # ===========================================================================
  
  cat("═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 1: Training Models\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  # Load NEFIN data
  nefin <- read_csv("data/processed/nefin_processed.csv", show_col_types = FALSE)
  
  # Load covariates
  nefin_climate <- read_csv("data/processed/climate_at_plots/nefin_climate.csv", 
                            show_col_types = FALSE) %>%
    distinct(CN, .keep_all = TRUE)
  nefin_ndvi <- read_csv("data/processed/ndvi_at_plots/nefin_ndvi.csv", 
                         show_col_types = FALSE) %>%
    distinct(CN, .keep_all = TRUE)
  
  nefin <- nefin %>%
    left_join(nefin_climate %>% select(CN, tmean, ppt), by = "CN") %>%
    left_join(nefin_ndvi %>% select(CN, ndvi_modis), by = "CN")
  
  # Filter complete cases
  lon_col <- names(nefin)[grepl("lon", names(nefin), ignore.case = TRUE)][1]
  lat_col <- names(nefin)[grepl("lat", names(nefin), ignore.case = TRUE)][1]
  
  nefin_complete <- nefin %>%
    filter(!is.na(aglb_Mg_per_ha), !is.na(ndvi_modis), !is.na(tmean), !is.na(ppt)) %>%
    rename(lon = !!lon_col, lat = !!lat_col)
  
  cat("  NEFIN complete cases:", nrow(nefin_complete), "\n")
  
  # Train NEFIN model
  X_nefin <- as.matrix(nefin_complete %>% select(ndvi_modis, tmean, ppt))
  y_nefin <- nefin_complete$aglb_Mg_per_ha
  
  model_nefin <- xgb.train(
    data = xgb.DMatrix(X_nefin, label = y_nefin),
    nrounds = 50,
    params = list(objective = "reg:squarederror", eta = 0.1, max_depth = 4),
    verbose = 0
  )
  cat("  ✓ NEFIN model trained\n")
  
  # Train FIA model (using one MC replicate)
  jitter_file <- "data/processed/mc_jitter_library/replicates/rep_0001.csv"
  if (fs::file_exists(jitter_file)) {
    jitter <- read_csv(jitter_file, show_col_types = FALSE)
    
    # Load FIA biomass
    tree_files <- list.files("data/interim/fia/states", pattern = "tree.csv",
                             recursive = TRUE, full.names = TRUE)
    if (length(tree_files) > 0) {
      trees <- lapply(tree_files, function(f) {
        read_csv(f, show_col_types = FALSE) %>%
          select(any_of(c("PLT_CN", "DRYBIO_AG", "TPA_UNADJ")))
      }) %>% bind_rows()
      
      plot_biomass <- trees %>%
        filter(!is.na(DRYBIO_AG), !is.na(TPA_UNADJ)) %>%
        group_by(PLT_CN) %>%
        summarise(aglb_Mg_per_ha = sum(DRYBIO_AG * TPA_UNADJ * 0.0004536 * 2.471), .groups = "drop") %>%
        rename(CN = PLT_CN)
      
      jitter <- jitter %>% inner_join(plot_biomass, by = "CN")
    }
    
    # Extract covariates at fuzzed locations
    jitter_sf <- st_as_sf(jitter, coords = c("lon_jittered", "lat_jittered"), crs = 4326)
    
    # Quick extraction
    r_ndvi <- terra::rast(ndvi_path)
    r_tmean <- terra::rast(tmean_path)
    r_ppt <- terra::rast(ppt_path)
    
    jitter_reproj <- st_transform(jitter_sf, crs(r_ndvi))
    jitter$ndvi_modis <- terra::extract(r_ndvi, terra::vect(jitter_reproj))[[2]]
    
    jitter_reproj <- st_transform(jitter_sf, crs(r_tmean))
    jitter$tmean <- terra::extract(r_tmean, terra::vect(jitter_reproj))[[2]]
    jitter$ppt <- terra::extract(r_ppt, terra::vect(jitter_reproj))[[2]]
    
    jitter_complete <- jitter %>%
      filter(!is.na(aglb_Mg_per_ha), !is.na(ndvi_modis), !is.na(tmean), !is.na(ppt))
    
    X_fia <- as.matrix(jitter_complete %>% select(ndvi_modis, tmean, ppt))
    y_fia <- jitter_complete$aglb_Mg_per_ha
    
    model_fia <- xgb.train(
      data = xgb.DMatrix(X_fia, label = y_fia),
      nrounds = 50,
      params = list(objective = "reg:squarederror", eta = 0.1, max_depth = 4),
      verbose = 0
    )
    cat("  ✓ FIA model trained\n")
  } else {
    model_fia <- model_nefin  # Fallback
    cat("  ⚠ Using NEFIN model as FIA fallback\n")
  }
  
  # ===========================================================================
  # 2. CREATE MODIS-RESOLUTION PREDICTION MAPS
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 2: Creating MODIS-Resolution Prediction Maps\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  cat("Creating NEFIN-trained prediction map...\n")
  r_pred_nefin <- create_modis_prediction_map(
    model = model_nefin,
    ndvi_path = ndvi_path,
    tmean_path = tmean_path,
    ppt_path = ppt_path,
    output_path = fs::path(output_dir, "modis_biomass_nefin.png"),
    title = "Predicted Biomass (NEFIN-trained model)"
  )
  
  cat("\nCreating FIA-trained prediction map...\n")
  r_pred_fia <- create_modis_prediction_map(
    model = model_fia,
    ndvi_path = ndvi_path,
    tmean_path = tmean_path,
    ppt_path = ppt_path,
    output_path = fs::path(output_dir, "modis_biomass_fia.png"),
    title = "Predicted Biomass (FIA-trained model)"
  )
  
  # ===========================================================================
  # 3. CREATE DIFFERENCE MAP
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 3: Creating Prediction Difference Map\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  r_diff <- create_prediction_difference_map(
    r_pred_nefin = r_pred_nefin,
    r_pred_fia = r_pred_fia,
    output_path = fs::path(output_dir, "modis_prediction_difference.png")
  )
  
  # ===========================================================================
  # 4. CREATE HEX-LEVEL ΔRMSE MAP
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 4: Creating Hex-Level ΔRMSE Map\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  hex_rmse <- create_hex_delta_rmse_map(
    nefin_data = nefin_complete,
    model_nefin = model_nefin,
    model_fia = model_fia,
    output_path = fs::path(output_dir, "hex_delta_rmse.png")
  )
  
  # ===========================================================================
  # 5. CREATE FUZZ PRESSURE PLOT
  # ===========================================================================
  
  cat("\n═══════════════════════════════════════════════════════════════════════════\n")
  cat("STEP 5: Creating Fuzz Pressure Plot\n")
  cat("═══════════════════════════════════════════════════════════════════════════\n\n")
  
  create_fuzz_pressure_plot(
    hex_rmse = hex_rmse,
    output_path = fs::path(output_dir, "fuzz_pressure_vs_delta_rmse.png")
  )
  
  # ===========================================================================
  # SUMMARY
  # ===========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║  MAP CREATION COMPLETE                                                    ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  cat("Output directory:", output_dir, "\n\n")
  
  cat("Maps created:\n")
  cat("  • modis_biomass_nefin.png/.tif - NEFIN-trained model prediction surface\n")
  cat("  • modis_biomass_fia.png/.tif - FIA-trained model prediction surface\n")
  cat("  • modis_prediction_difference.png/.tif - Difference (NEFIN - FIA)\n")
  cat("  • hex_delta_rmse.png/.csv - Hex-level ΔRMSE map\n")
  cat("  • fuzz_pressure_vs_delta_rmse.png - Sample density vs ΔRMSE\n")
  
  invisible(list(
    r_pred_nefin = r_pred_nefin,
    r_pred_fia = r_pred_fia,
    r_diff = r_diff,
    hex_rmse = hex_rmse
  ))
}

# =============================================================================
# CLI
# =============================================================================

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  create_all_prediction_maps()
}
