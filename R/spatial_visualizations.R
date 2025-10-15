# R/spatial_visualizations.R
# Create maps using hex grid GeoJSON files

suppressPackageStartupMessages({
  library(sf); library(readr); library(dplyr); library(ggplot2)
  library(viridis); library(patchwork); library(fs); library(scales)
})

spatial_visualizations <- function(consolidated_dir = NULL,
                                   output_dir = NULL,
                                   grid_name = "fia") {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  Spatial Visualizations with Hex Grids                   ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Find most recent consolidated directory
  if (is.null(consolidated_dir)) {
    all_dirs <- list.dirs("runs", recursive = FALSE)
    consol_dirs <- all_dirs[grepl("^runs/consolidated_", all_dirs)]
    if (length(consol_dirs) == 0) {
      stop("No consolidated results found.")
    }
    consol_dirs <- consol_dirs[order(file.mtime(consol_dirs), decreasing = TRUE)]
    consolidated_dir <- consol_dirs[1]
  }
  
  cat("Using results from:", consolidated_dir, "\n")
  
  if (is.null(output_dir)) {
    output_dir <- fs::path(consolidated_dir, "spatial_maps")
  }
  fs::dir_create(output_dir, recurse = TRUE)
  
  # Turn off spherical geometry for faster processing
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  # Load data
  cat("Loading data...\n")
  fia_file <- fs::path(consolidated_dir, "fia_all_scales.csv")
  nefin_file <- fs::path(consolidated_dir, "fia_nefin_comparison_all_scales.csv")
  
  fia <- readr::read_csv(fia_file, show_col_types = FALSE)
  
  has_nefin <- fs::file_exists(nefin_file)
  nefin <- if (has_nefin) {
    readr::read_csv(nefin_file, show_col_types = FALSE)
  } else NULL
  
  # Load hex grids
  cat("Loading hex grids...\n")
  
  hex_paths <- list(
    fia = "data/hex/hex_grid.geojson",
    "1.5k" = "data/hex/hex_grid_1.5kac.geojson",
    "3k" = "data/hex/hex_grid_3kac.geojson",
    "6k" = "data/hex/hex_grid_6kac.geojson"
  )
  
  # Try to find available grids
  available_grids <- list()
  for (name in names(hex_paths)) {
    path <- hex_paths[[name]]
    if (fs::file_exists(path)) {
      cat("  ✓ Loading", name, "grid from", path, "\n")
      grid <- sf::st_read(path, quiet = TRUE)
      
      # Ensure hex_id exists
      if (!("hex_id" %in% names(grid))) {
        if ("ID" %in% names(grid)) grid <- dplyr::rename(grid, hex_id = ID)
        else if ("OBJECTID" %in% names(grid)) grid <- dplyr::rename(grid, hex_id = OBJECTID)
        else grid$hex_id <- seq_len(nrow(grid))
      }
      grid$hex_id <- as.character(grid$hex_id)
      
      available_grids[[name]] <- grid
    } else {
      cat("  ⚠ Grid not found:", name, "at", path, "\n")
    }
  }
  
  if (length(available_grids) == 0) {
    stop("No hex grids found. Check that GeoJSON files exist in data/hex/")
  }
  
  # ========================================================================
  # MAP SET 1: FIA Error Maps by Grid Scale
  # ========================================================================
  
  cat("\nCreating FIA error maps...\n")
  
  for (grid_name in names(available_grids)) {
    cat("  → Maps for", grid_name, "grid\n")
    
    grid_sf <- available_grids[[grid_name]]
    
    # Filter FIA data for this grid
    fia_grid <- fia |>
      dplyr::filter(grid_scale == grid_name) |>
      dplyr::group_by(hex_id) |>
      dplyr::summarise(
        mean_biomass = mean(mean, na.rm = TRUE),
        mean_se = mean(se, na.rm = TRUE),
        mean_pos_sd = mean(positional_sd, na.rm = TRUE),
        mean_total_sd = mean(total_sd, na.rm = TRUE),
        pos_fraction = mean(positional_sd / total_sd, na.rm = TRUE),
        mean_n_plots = mean(n_plots, na.rm = TRUE),
        n_years = dplyr::n(),
        .groups = "drop"
      )
    
    # Join with spatial
    map_data <- grid_sf |>
      dplyr::left_join(fia_grid, by = "hex_id")
    
    # Map 1: Mean biomass
    p1 <- ggplot(map_data) +
      geom_sf(aes(fill = mean_biomass), color = NA, linewidth = 0) +
      scale_fill_viridis_c(option = "G", na.value = "grey90", 
                           labels = comma, name = "AGLB\n(Mg/ha)") +
      labs(title = paste0("Mean Aboveground Biomass - ", toupper(grid_name), " Grid"),
           subtitle = "Averaged across all years") +
      theme_void(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "right")
    
    ggsave(fs::path(output_dir, paste0("map_", grid_name, "_01_biomass.png")), p1,
           width = 10, height = 8, dpi = 300)
    
    # Map 2: Sampling error
    p2 <- ggplot(map_data) +
      geom_sf(aes(fill = mean_se), color = NA, linewidth = 0) +
      scale_fill_viridis_c(option = "B", na.value = "grey90",
                           name = "SE\n(Mg/ha)") +
      labs(title = paste0("Sampling Error - ", toupper(grid_name), " Grid"),
           subtitle = "Uncertainty from plot variability") +
      theme_void(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "right")
    
    ggsave(fs::path(output_dir, paste0("map_", grid_name, "_02_sampling_error.png")), p2,
           width = 10, height = 8, dpi = 300)
    
    # Map 3: Positional SD
    p3 <- ggplot(map_data) +
      geom_sf(aes(fill = mean_pos_sd), color = NA, linewidth = 0) +
      scale_fill_viridis_c(option = "E", na.value = "grey90",
                           name = "Pos SD\n(Mg/ha)") +
      labs(title = paste0("Positional Uncertainty - ", toupper(grid_name), " Grid"),
           subtitle = "Uncertainty from coordinate fuzzing") +
      theme_void(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "right")
    
    ggsave(fs::path(output_dir, paste0("map_", grid_name, "_03_positional_sd.png")), p3,
           width = 10, height = 8, dpi = 300)
    
    # Map 4: Total uncertainty
    p4 <- ggplot(map_data) +
      geom_sf(aes(fill = mean_total_sd), color = NA, linewidth = 0) +
      scale_fill_viridis_c(option = "A", na.value = "grey90",
                           name = "Total SD\n(Mg/ha)") +
      labs(title = paste0("Total Uncertainty - ", toupper(grid_name), " Grid"),
           subtitle = "Combined sampling + positional error") +
      theme_void(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "right")
    
    ggsave(fs::path(output_dir, paste0("map_", grid_name, "_04_total_uncertainty.png")), p4,
           width = 10, height = 8, dpi = 300)
    
    # Map 5: Positional error fraction
    p5 <- ggplot(map_data) +
      geom_sf(aes(fill = pos_fraction * 100), color = NA, linewidth = 0) +
      scale_fill_viridis_c(option = "H", na.value = "grey90",
                           limits = c(0, 100),
                           name = "Pos %\nof Total") +
      labs(title = paste0("Positional Error Contribution - ", toupper(grid_name), " Grid"),
           subtitle = "% of total uncertainty from coordinate fuzzing") +
      theme_void(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "right")
    
    ggsave(fs::path(output_dir, paste0("map_", grid_name, "_05_positional_fraction.png")), p5,
           width = 10, height = 8, dpi = 300)
    
    # Map 6: Sample size
    p6 <- ggplot(map_data) +
      geom_sf(aes(fill = mean_n_plots), color = NA, linewidth = 0) +
      scale_fill_viridis_c(option = "C", na.value = "grey90",
                           trans = "log10", labels = comma,
                           name = "Plots\nper Hex") +
      labs(title = paste0("Sample Size - ", toupper(grid_name), " Grid"),
           subtitle = "Mean number of FIA plots per hex") +
      theme_void(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "right")
    
    ggsave(fs::path(output_dir, paste0("map_", grid_name, "_06_sample_size.png")), p6,
           width = 10, height = 8, dpi = 300)
    
    # Multi-panel overview
    overview <- (p1 | p4) / (p2 | p3) +
      plot_annotation(
        title = paste0("FIA Error Analysis Overview - ", toupper(grid_name), " Grid"),
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
      )
    
    ggsave(fs::path(output_dir, paste0("map_", grid_name, "_00_overview.png")), overview,
           width = 16, height = 14, dpi = 300)
  }
  
  # ========================================================================
  # MAP SET 2: FIA vs NEFIN Comparison Maps
  # ========================================================================
  
  if (!is.null(nefin)) {
    cat("\nCreating FIA vs NEFIN comparison maps...\n")
    
    for (grid_name in names(available_grids)) {
      cat("  → Comparison maps for", grid_name, "grid\n")
      
      grid_sf <- available_grids[[grid_name]]
      
      # Filter and aggregate NEFIN comparison data
      nefin_grid <- nefin |>
        dplyr::filter(grid_scale == grid_name) |>
        dplyr::group_by(hex_id) |>
        dplyr::summarise(
          mean_fia = mean(mean_fia, na.rm = TRUE),
          mean_nefin = mean(mean_nefin, na.rm = TRUE),
          mean_diff = mean(diff, na.rm = TRUE),
          abs_diff = mean(abs(diff), na.rm = TRUE),
          has_both = any(has_both),
          has_fia = any(has_fia),
          has_nefin = any(has_nefin),
          n_plots_fia = mean(n_plots_fia, na.rm = TRUE),
          n_plots_nefin = mean(n_plots_nefin, na.rm = TRUE),
          .groups = "drop"
        )
      
      # Join with spatial
      comp_map_data <- grid_sf |>
        dplyr::left_join(nefin_grid, by = "hex_id")
      
      # Map 7: Coverage
      coverage_data <- comp_map_data |>
        dplyr::mutate(
          coverage = dplyr::case_when(
            has_both ~ "Both",
            has_fia ~ "FIA Only",
            has_nefin ~ "NEFIN Only",
            TRUE ~ "Neither"
          )
        )
      
      p7 <- ggplot(coverage_data) +
        geom_sf(aes(fill = coverage), color = "grey40", linewidth = 0.1) +
        scale_fill_manual(
          values = c("Both" = "#2E8B57", "FIA Only" = "#4682B4", 
                     "NEFIN Only" = "#DAA520", "Neither" = "grey90"),
          na.value = "grey90",
          name = "Data Coverage"
        ) +
        labs(title = paste0("Data Coverage - ", toupper(grid_name), " Grid"),
             subtitle = "Which hexes have FIA, NEFIN, or both") +
        theme_void(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              legend.position = "right")
      
      ggsave(fs::path(output_dir, paste0("map_", grid_name, "_07_coverage.png")), p7,
             width = 10, height = 8, dpi = 300)
      
      # Map 8: Bias (NEFIN - FIA)
      p8 <- ggplot(comp_map_data |> dplyr::filter(has_both)) +
        geom_sf(aes(fill = mean_diff), color = NA, linewidth = 0) +
        scale_fill_gradient2(
          low = "#2166AC", mid = "white", high = "#B2182B",
          midpoint = 0, na.value = "grey90",
          name = "Bias\n(Mg/ha)",
          labels = comma
        ) +
        labs(title = paste0("Bias (NEFIN - FIA) - ", toupper(grid_name), " Grid"),
             subtitle = "Positive = NEFIN higher, Negative = FIA higher") +
        theme_void(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              legend.position = "right")
      
      ggsave(fs::path(output_dir, paste0("map_", grid_name, "_08_bias.png")), p8,
             width = 10, height = 8, dpi = 300)
      
      # Map 9: Absolute difference
      p9 <- ggplot(comp_map_data |> dplyr::filter(has_both)) +
        geom_sf(aes(fill = abs_diff), color = NA, linewidth = 0) +
        scale_fill_viridis_c(option = "A", na.value = "grey90",
                             name = "Abs Diff\n(Mg/ha)") +
        labs(title = paste0("Disagreement Magnitude - ", toupper(grid_name), " Grid"),
             subtitle = "Absolute difference between FIA and NEFIN") +
        theme_void(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              legend.position = "right")
      
      ggsave(fs::path(output_dir, paste0("map_", grid_name, "_09_abs_diff.png")), p9,
             width = 10, height = 8, dpi = 300)
      
      # Map 10: Ratio of NEFIN to FIA plots
      comp_map_data_ratio <- comp_map_data |>
        dplyr::filter(has_both) |>
        dplyr::mutate(
          plot_ratio = n_plots_nefin / (n_plots_fia + n_plots_nefin),
          ratio_category = dplyr::case_when(
            plot_ratio < 0.2 ~ "FIA Dominant (>80%)",
            plot_ratio < 0.4 ~ "FIA Majority (60-80%)",
            plot_ratio < 0.6 ~ "Balanced (40-60%)",
            plot_ratio < 0.8 ~ "NEFIN Majority (60-80%)",
            TRUE ~ "NEFIN Dominant (>80%)"
          )
        )
      
      p10 <- ggplot(comp_map_data_ratio) +
        geom_sf(aes(fill = plot_ratio), color = NA, linewidth = 0) +
        scale_fill_gradient2(
          low = "#4575B4", mid = "white", high = "#D73027",
          midpoint = 0.5, na.value = "grey90",
          name = "NEFIN\nFraction",
          labels = percent
        ) +
        labs(title = paste0("Sample Balance - ", toupper(grid_name), " Grid"),
             subtitle = "NEFIN plots as fraction of total (FIA + NEFIN)") +
        theme_void(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              legend.position = "right")
      
      ggsave(fs::path(output_dir, paste0("map_", grid_name, "_10_sample_balance.png")), p10,
             width = 10, height = 8, dpi = 300)
      
      # Comparison overview
      comp_overview <- (p7 | p8) / (p9 | p10) +
        plot_annotation(
          title = paste0("FIA vs NEFIN Comparison - ", toupper(grid_name), " Grid"),
          theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
        )
      
      ggsave(fs::path(output_dir, paste0("map_", grid_name, "_00_comparison_overview.png")), 
             comp_overview, width = 16, height = 14, dpi = 300)
    }
  }
  
  # ========================================================================
  # MAP SET 3: Multi-scale comparison
  # ========================================================================
  
  if (length(available_grids) > 1) {
    cat("\nCreating multi-scale comparison maps...\n")
    
    # Create biomass comparison across scales
    biomass_maps <- list()
    
    for (i in seq_along(available_grids)) {
      grid_name <- names(available_grids)[i]
      grid_sf <- available_grids[[grid_name]]
      
      fia_grid <- fia |>
        dplyr::filter(grid_scale == grid_name) |>
        dplyr::group_by(hex_id) |>
        dplyr::summarise(mean_biomass = mean(mean, na.rm = TRUE), .groups = "drop")
      
      map_data <- grid_sf |> dplyr::left_join(fia_grid, by = "hex_id")
      
      p <- ggplot(map_data) +
        geom_sf(aes(fill = mean_biomass), color = NA) +
        scale_fill_viridis_c(option = "G", na.value = "grey90",
                             limits = range(fia$mean, na.rm = TRUE),
                             name = "AGLB\n(Mg/ha)") +
        labs(title = toupper(grid_name)) +
        theme_void(base_size = 10) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              legend.position = "bottom",
              legend.key.width = unit(1.5, "cm"))
      
      biomass_maps[[i]] <- p
    }
    
    # Arrange in grid
    multi_biomass <- wrap_plots(biomass_maps, ncol = 2) +
      plot_annotation(
        title = "Mean Biomass Across Grid Scales",
        subtitle = "Same color scale for direct comparison",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5))
      )
    
    ggsave(fs::path(output_dir, "map_multiscale_biomass.png"), multi_biomass,
           width = 14, height = 12, dpi = 300)
    
    # Positional SD comparison
    pos_maps <- list()
    
    for (i in seq_along(available_grids)) {
      grid_name <- names(available_grids)[i]
      grid_sf <- available_grids[[grid_name]]
      
      fia_grid <- fia |>
        dplyr::filter(grid_scale == grid_name) |>
        dplyr::group_by(hex_id) |>
        dplyr::summarise(mean_pos_sd = mean(positional_sd, na.rm = TRUE), .groups = "drop")
      
      map_data <- grid_sf |> dplyr::left_join(fia_grid, by = "hex_id")
      
      p <- ggplot(map_data) +
        geom_sf(aes(fill = mean_pos_sd), color = NA) +
        scale_fill_viridis_c(option = "E", na.value = "grey90",
                             limits = range(fia$positional_sd, na.rm = TRUE),
                             name = "Pos SD\n(Mg/ha)") +
        labs(title = toupper(grid_name)) +
        theme_void(base_size = 10) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              legend.position = "bottom",
              legend.key.width = unit(1.5, "cm"))
      
      pos_maps[[i]] <- p
    }
    
    multi_pos <- wrap_plots(pos_maps, ncol = 2) +
      plot_annotation(
        title = "Positional Uncertainty Across Grid Scales",
        subtitle = "Does coarser grid = more positional error?",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5))
      )
    
    ggsave(fs::path(output_dir, "map_multiscale_positional.png"), multi_pos,
           width = 14, height = 12, dpi = 300)
  }
  
  # ========================================================================
  # DONE
  # ========================================================================
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  SPATIAL VISUALIZATIONS COMPLETE                          ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  cat("Output directory:\n")
  cat("  ", output_dir, "\n")
  cat("\n")
  cat("Created maps:\n")
  
  map_files <- list.files(output_dir, pattern = "\\.png$")
  cat("  Total maps: ", length(map_files), "\n\n")
  
  for (grid_name in names(available_grids)) {
    grid_maps <- map_files[grepl(paste0("^map_", grid_name), map_files)]
    cat("  ", toupper(grid_name), "grid: ", length(grid_maps), " maps\n")
  }
  
  cat("\n")
  cat("Map types per grid:\n")
  cat("  00 - Multi-panel overview\n")
  cat("  01 - Mean biomass\n")
  cat("  02 - Sampling error\n")
  cat("  03 - Positional SD\n")
  cat("  04 - Total uncertainty\n")
  cat("  05 - Positional fraction\n")
  cat("  06 - Sample size\n")
  
  if (!is.null(nefin)) {
    cat("  07 - Data coverage (FIA/NEFIN/Both)\n")
    cat("  08 - Bias (NEFIN - FIA)\n")
    cat("  09 - Absolute difference\n")
    cat("  10 - Sample balance\n")
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
  
  spatial_visualizations(
    consolidated_dir = get_arg("--dir", NULL),
    output_dir = get_arg("--out", NULL),
    grid_name = get_arg("--grid", "fia")
  )
}