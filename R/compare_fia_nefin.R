# R/compare_fia_nefin.R
# Compare FIA and NEFIN biomass estimates by hex grid

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(sf); library(fs)
})

source("R/utils_metrics.R")

`%||%` <- function(a,b) if (!is.null(a)) a else b

compare_fia_nefin <- function(fia_results = NULL,
                              nefin_assignments = "data/processed/nefin_hex_assignments.csv",
                              hex_grid_name = "fia",
                              years = 2018:2020,
                              level_window = 3,
                              out_dir = "runs/fia_nefin_comparison") {
  
  message("\n═══════════════════════════════════════════════════════════")
  message("FIA vs NEFIN COMPARISON")
  message("═══════════════════════════════════════════════════════════\n")
  
  fs::dir_create(out_dir, recurse = TRUE)
  
  # Find FIA results if not provided
  if (is.null(fia_results)) {
    # Look for most recent FIA run
    runs_dir <- "runs"
    run_dirs <- list.dirs(runs_dir, recursive = FALSE)
    fia_runs <- run_dirs[grepl(paste0("_", hex_grid_name, "_"), run_dirs)]
    
    if (!length(fia_runs)) {
      stop("No FIA results found for grid: ", hex_grid_name,
           "\n  Run: Rscript run_pipeline.R --compute --grid=", hex_grid_name)
    }
    
    # Get most recent
    fia_runs <- fia_runs[order(file.mtime(fia_runs), decreasing = TRUE)]
    fia_run <- fia_runs[1]
    fia_results <- fs::path(fia_run, "hex_aglb_results.csv")
    
    message("→ Using FIA results: ", fia_results)
  }
  
  if (!fs::file_exists(fia_results)) {
    stop("FIA results not found: ", fia_results)
  }
  
  # Load data
  message("→ Loading FIA results...")
  fia <- readr::read_csv(fia_results, show_col_types = FALSE)
  message("  Hex-years: ", nrow(fia))
  
  message("→ Loading NEFIN assignments...")
  if (!fs::file_exists(nefin_assignments)) {
    stop("NEFIN assignments not found: ", nefin_assignments,
         "\n  Run: Rscript R/assign_nefin_to_hex.R")
  }
  
  nefin <- readr::read_csv(nefin_assignments, show_col_types = FALSE)
  message("  Plots: ", nrow(nefin))
  
  # Get hex column for this grid
  hex_col <- paste0("hex_id_", hex_grid_name)
  
  if (!(hex_col %in% names(nefin))) {
    stop("NEFIN assignments missing column: ", hex_col,
         "\n  Available grids: ", paste(gsub("hex_id_", "", names(nefin)[grepl("^hex_id_", names(nefin))]), collapse = ", "))
  }
  
  # Aggregate NEFIN to hex level (matching FIA structure)
  message("\n→ Aggregating NEFIN by hex...")
  
  nefin_hex <- nefin |>
    dplyr::rename(hex_id = !!rlang::sym(hex_col)) |>
    dplyr::filter(!is.na(hex_id), !is.na(aglb_Mg_per_ha))
  
  # Create comparable time windows
  nefin_by_year <- lapply(years, function(yr) {
    year_window <- (yr - floor(level_window/2)):(yr + floor(level_window/2))
    
    nefin_hex |>
      dplyr::filter(MEASYEAR %in% year_window) |>
      dplyr::group_by(hex_id) |>
      dplyr::summarise(
        mean_nefin = mean(aglb_Mg_per_ha, na.rm = TRUE),
        se_nefin = sd(aglb_Mg_per_ha, na.rm = TRUE) / sqrt(dplyr::n()),
        n_plots_nefin = dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        year_label = yr,
        window = paste0(level_window, "y"),
        source = "NEFIN"
      )
  })
  
  nefin_summary <- dplyr::bind_rows(nefin_by_year)
  
  message("  NEFIN hex-years: ", nrow(nefin_summary))
  
  # Join FIA and NEFIN
  message("\n→ Joining FIA and NEFIN by hex...")
  
  comparison <- fia |>
    dplyr::select(hex_id, year_label, window, 
                  mean_fia = mean, se_fia = se, n_plots_fia = n_plots,
                  positional_sd, total_sd) |>
    dplyr::full_join(
      nefin_summary |> dplyr::select(hex_id, year_label, window,
                                     mean_nefin, se_nefin, n_plots_nefin),
      by = c("hex_id", "year_label", "window")
    ) |>
    dplyr::mutate(
      # Calculate difference
      diff = mean_nefin - mean_fia,
      pct_diff = 100 * diff / mean_fia,
      
      # Data availability
      has_fia = !is.na(mean_fia),
      has_nefin = !is.na(mean_nefin),
      has_both = has_fia & has_nefin
    )
  
  # Summary stats
  message("\n=== Comparison Summary ===")
  message("  Total hex-years: ", nrow(comparison))
  message("  FIA only: ", sum(comparison$has_fia & !comparison$has_nefin))
  message("  NEFIN only: ", sum(!comparison$has_fia & comparison$has_nefin))
  message("  Both: ", sum(comparison$has_both))
  
  if (sum(comparison$has_both) > 0) {
    both_data <- comparison |> dplyr::filter(has_both)
    
    stats <- both_data |>
      dplyr::summarise(
        n_hexes = dplyr::n(),
        mean_fia = mean(mean_fia, na.rm = TRUE),
        mean_nefin = mean(mean_nefin, na.rm = TRUE),
        mean_diff = mean(diff, na.rm = TRUE),
        rmse = sqrt(mean(diff^2, na.rm = TRUE)),
        correlation = cor(mean_fia, mean_nefin, use = "complete.obs")
      )
    
    message("\n  Agreement statistics (hexes with both):")
    message("    Mean FIA: ", round(stats$mean_fia, 2), " Mg/ha")
    message("    Mean NEFIN: ", round(stats$mean_nefin, 2), " Mg/ha")
    message("    Mean difference: ", round(stats$mean_diff, 2), " Mg/ha")
    message("    RMSE: ", round(stats$rmse, 2), " Mg/ha")
    message("    Correlation: ", round(stats$correlation, 3))
    
    # Write outputs
    readr::write_csv(stats, fs::path(out_dir, "comparison_stats.csv"))
    readr::write_csv(comparison, fs::path(out_dir, "hex_comparison.csv"))
    
    # Create scatter plot
    p1 <- ggplot(both_data, aes(x = mean_fia, y = mean_nefin)) +
      geom_point(alpha = 0.5, color = "steelblue") +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      geom_smooth(method = "lm", color = "blue", se = FALSE) +
      labs(title = paste0("FIA vs NEFIN Biomass (", hex_grid_name, " grid)"),
           subtitle = paste0("Correlation: ", round(stats$correlation, 3), 
                             ", RMSE: ", round(stats$rmse, 1), " Mg/ha"),
           x = "FIA AGLB (Mg/ha)", y = "NEFIN AGLB (Mg/ha)") +
      theme_minimal(base_size = 14)
    
    ggsave(fs::path(out_dir, "scatter_plot.png"), p1, width = 8, height = 7, dpi = 300)
    
    # Create difference histogram
    p2 <- ggplot(both_data, aes(x = diff)) +
      geom_histogram(bins = 30, fill = "steelblue", color = "white") +
      geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
      geom_vline(xintercept = stats$mean_diff, color = "blue", linetype = "solid", linewidth = 1) +
      labs(title = "Distribution of Differences (NEFIN - FIA)",
           subtitle = paste0("Mean bias: ", round(stats$mean_diff, 1), " Mg/ha"),
           x = "Difference (Mg/ha)", y = "Count") +
      theme_minimal(base_size = 14)
    
    ggsave(fs::path(out_dir, "difference_hist.png"), p2, width = 8, height = 6, dpi = 300)
    
    message("\n✓ Created plots: scatter_plot.png, difference_hist.png")
  } else {
    message("\n⚠ No hexes with both FIA and NEFIN data - cannot compute statistics")
  }
  
  message("\n✓ Comparison complete!")
  message("  Output directory: ", out_dir)
  
  invisible(list(
    comparison = comparison,
    out_dir = out_dir
  ))
}

# CLI
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  get_arg <- function(flag, default = NULL) {
    hit <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
  }
  
  compare_fia_nefin(
    fia_results = get_arg("--fia", NULL),
    nefin_assignments = get_arg("--nefin", "data/processed/nefin_hex_assignments.csv"),
    hex_grid_name = get_arg("--grid", "fia"),
    years = eval(parse(text = get_arg("--years", "2018:2020"))),
    level_window = as.integer(get_arg("--window", "3")),
    out_dir = get_arg("--out", "runs/fia_nefin_comparison")
  )
}