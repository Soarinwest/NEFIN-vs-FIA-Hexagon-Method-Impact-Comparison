suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(sf)
  library(fs)
  library(purrr)
  library(stringr)
  library(yaml)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

# Source scale name utilities
if (file.exists("R/utils_scale_names.R")) {
  source("R/utils_scale_names.R")
} else {
  # Inline fallback

  standardize_scale_name <- function(x) {
    x <- as.character(x)
    x[tolower(x) == "fia"] <- "64kha"
    x
  }
}

# ------------------------------------------------------------
# Load process.yml
# ------------------------------------------------------------
load_config <- function(cfg_path = "configs/process.yml") {
  yaml::read_yaml(cfg_path)
}

# ------------------------------------------------------------
# Helper: which column in NEFIN to use for a given grid
# e.g. grid "100ha" -> "hex_id_100ha"
# Special case: FIA grid might just be "hex_id" or "hex_id_fia"
# We'll support both "64kha" and FIA if you eventually tag it.
# ------------------------------------------------------------
hex_col_for_grid <- function(grid_name, nefin_cols) {
  # preferred name:
  candidate <- paste0("hex_id_", grid_name)
  
  if (candidate %in% nefin_cols) return(candidate)
  
  # fallback for "fia" style grids
  if (grid_name %in% c("fia", "64kha")) {
    if ("hex_id_fia" %in% nefin_cols) return("hex_id_fia")
    if ("hex_id" %in% nefin_cols)     return("hex_id")
  }
  
  # If nothing matches, throw a useful error
  stop("Could not find matching NEFIN hex column for grid '",
       grid_name, "'. Tried: ", candidate,
       " and fia fallbacks. NEFIN columns available are: ",
       paste(nefin_cols, collapse=", "))
}

# ------------------------------------------------------------
# Helper: build FIA run dir for a grid
#
# We were previously manually passing:
#   runs/2025-10-31_aglb_100ha_W5y/hex_aglb_results.csv
#
# We can infer that pattern from cfg$metric ("aglb"), grid name, and window size.
# We'll assume:
#   runs/<DATE>_<metric>_<grid>_W<window>y/hex_aglb_results.csv
#
# BUT: We don't actually know <DATE> from config. You might have multiple runs.
# So: we'll search runs/ for directories that look like they match.
#
# Example dir regex:
#   "^.*/[0-9]{4}-[0-9]{2}-[0-9]{2}_.+_<grid>_W[0-9]+y$"
#
# We'll pick the newest by mtime.
#
# NOTE: Handles both "fia" and "64kha" naming - searches for both
# ------------------------------------------------------------
find_fia_results_for_grid <- function(grid_name, metric, level_window) {
  runs_dir <- "runs"
  run_dirs <- list.dirs(runs_dir, recursive = FALSE)
  
  # window like "W5y"
  window_tag <- paste0("W", level_window, "y")
  
  # Build list of names to search for (handle fia <-> 64kha mapping)
  search_names <- grid_name
  if (tolower(grid_name) == "fia") {
    search_names <- c("fia", "64kha")
  } else if (tolower(grid_name) == "64kha") {
    search_names <- c("64kha", "fia")
  }
  
  # Search for each possible name
  candidates <- character(0)
  for (name in search_names) {
    tag <- paste0("_", name, "_", window_tag)
    matches <- run_dirs[grepl(tag, run_dirs, ignore.case = TRUE)]
    candidates <- c(candidates, matches)
  }
  
  candidates <- unique(candidates)
  
  if (!length(candidates)) {
    stop("No FIA results directory matched grid='", grid_name,
         "' and window '", window_tag,
         "'. Looked for '*_", paste(search_names, collapse = "_' or '*_"), 
         "_", window_tag, "' in ", runs_dir)
  }
  
  # pick most recent
  candidates <- candidates[order(file.mtime(candidates), decreasing = TRUE)]
  fia_dir <- candidates[1]
  
  fia_results_path <- fs::path(fia_dir, "hex_aglb_results.csv")
  if (!fs::file_exists(fia_results_path)) {
    stop("hex_aglb_results.csv not found in ", fia_dir)
  }
  
  list(
    fia_dir = fia_dir,
    fia_results_path = fia_results_path
  )
}

# ------------------------------------------------------------
# summarize FIA temporal stability (unchanged logic)
# ------------------------------------------------------------
summarize_fia_stability <- function(fia_df, scale_name, out_dir) {
  
  stability <- fia_df %>%
    group_by(hex_id) %>%
    summarise(
      n_years_covered      = n_distinct(year_label),
      mean_FIA_overall     = mean(mean, na.rm = TRUE),
      sd_FIA_across_years  = sd(mean, na.rm = TRUE),
      mean_SE_FIA          = mean(se, na.rm = TRUE),
      mean_positional_sd   = mean(positional_sd, na.rm = TRUE),
      mean_total_sd        = mean(total_sd, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(stability, fs::path(out_dir, paste0("fia_stability_", scale_name, ".csv")))
  
  # histogram of FIA temporal SD
  p_sd <- ggplot(stability, aes(x = sd_FIA_across_years)) +
    geom_histogram(bins = 30, fill = "gray60", color = "white") +
    labs(
      title = paste0("FIA Temporal Variability (", scale_name, " hexes)"),
      x = "SD of FIA biomass across FIA year_labels (Mg/ha)",
      y = "Count of hexes"
    ) +
    theme_minimal(base_size = 14)
  
  ggsave(
    filename = fs::path(out_dir, paste0("fia_temporal_sd_hist_", scale_name, ".png")),
    plot = p_sd,
    width = 8, height = 6, dpi = 300
  )
  
  stability
}

# ------------------------------------------------------------
# main 5-year comparison with recommendations, fuzz, bias
# now takes cfg for metadata
# FIXED: Proper SE combination using inverse-variance weighting
# ------------------------------------------------------------
compare_5y_window <- function(fia_df,
                              nefin_df,
                              cfg,
                              grid_name,
                              hex_col,
                              out_dir) {
  
  focal_years <- cfg$years
  level_window <- cfg$level_window
  jitter_radius_m <- cfg$jitter_radius_m %||% cfg$mc_reps # fallback silly but safe
  mc_reps <- cfg$mc_reps %||% NA_real_
  
  target_n <- 3          # FIA plots in 5y considered "good coverage"
  big_shift_threshold <- 10      # Mg/ha considered a large shift
  bias_concentration_threshold <- 0.8
  se_drop_good <- 0.15   # relative SE drop considered meaningful
  se_drop_small <- 0.05  # relative SE drop considered trivial
  
  # Minimum SE to consider valid for ratio calculations
  min_se_for_ratio <- 0.1  # Avoid division by near-zero SE
  
  # --- Prep NEFIN for these years and this grid ---
  nefin_xy <- nefin_df %>%
    rename(hex_id = !!sym(hex_col)) %>%
    filter(
      !is.na(hex_id),
      !is.na(aglb_Mg_per_ha),
      MEASYEAR %in% focal_years
    )
  
  # source concentration summary (bias risk)
  nefin_source_summary <- nefin_xy %>%
    group_by(hex_id, source) %>%
    summarise(
      n_source_plots = n(),
      .groups = "drop"
    ) %>%
    group_by(hex_id) %>%
    mutate(
      total_plots_hex = sum(n_source_plots),
      source_frac = n_source_plots / total_plots_hex
    ) %>%
    summarise(
      n_nefin_sources = n_distinct(source),
      max_source_frac = max(source_frac, na.rm = TRUE),
      .groups = "drop"
    )
  
  # collapse FIA over focal years
  fia_5y <- fia_df %>%
    filter(year_label %in% focal_years) %>%
    group_by(hex_id) %>%
    summarise(
      fia_only_mean     = mean(mean, na.rm = TRUE),
      fia_only_se       = mean(se, na.rm = TRUE),
      fia_only_n        = sum(n_plots, na.rm = TRUE),
      positional_sd_5y  = mean(positional_sd, na.rm = TRUE),
      total_sd_5y       = mean(total_sd, na.rm = TRUE),
      .groups = "drop"
    )
  
  # collapse NEFIN over focal years
  nefin_5y <- nefin_xy %>%
    group_by(hex_id) %>%
    summarise(
      nefin_mean_5y = mean(aglb_Mg_per_ha, na.rm = TRUE),
      nefin_se_5y   = sd(aglb_Mg_per_ha, na.rm = TRUE) / sqrt(n()),
      nefin_n_5y    = n(),
      .groups = "drop"
    )
  
  # join all hex-level summaries
  joined_5y <- full_join(fia_5y, nefin_5y, by = "hex_id") %>%
    full_join(nefin_source_summary, by = "hex_id")
  
  # =========================================================================
  # FIXED: Row-wise calculation with PROPER inverse-variance SE weighting
  # =========================================================================
  safe_row_calc <- function(fia_only_mean, fia_only_se, fia_only_n,
                            nefin_mean_5y, nefin_se_5y, nefin_n_5y) {
    
    n_fia   <- ifelse(is.na(fia_only_n),   0, fia_only_n)
    n_nefin <- ifelse(is.na(nefin_n_5y),   0, nefin_n_5y)
    n_combined <- n_fia + n_nefin
    
    # Weighting for MEAN: FIA is anchor unless FIA is weak
    # (This controls how much NEFIN pulls the estimate)
    w_fia   <- min(1, n_fia / target_n)
    w_nefin <- 1 - w_fia
    
    # augmented mean: weighted blend (FIA-anchored)
    augmented_mean <- NA_real_
    if (!is.na(fia_only_mean) && !is.na(nefin_mean_5y)) {
      augmented_mean <- w_fia * fia_only_mean + w_nefin * nefin_mean_5y
    } else if (!is.na(fia_only_mean)) {
      augmented_mean <- fia_only_mean
    } else if (!is.na(nefin_mean_5y)) {
      augmented_mean <- nefin_mean_5y
    }
    
    # =========================================================================
    # FIXED: augmented SE using INVERSE-VARIANCE WEIGHTING
    # 
    # When combining independent estimates, the combined SE is:
    #   SE_combined = 1 / sqrt(1/SE_fia^2 + 1/SE_nefin^2)
    # 
    # This is ALWAYS lower than either individual SE (which is the point!)
    # The old weighted-blend approach was wrong because it would return
    # augmented_se = fia_only_se when FIA had good coverage.
    # =========================================================================
    augmented_se <- NA_real_
    
    if (!is.na(fia_only_se) && !is.na(nefin_se_5y) && 
        fia_only_se > 0 && nefin_se_5y > 0) {
      # Inverse-variance weighting: proper statistical combination
      # This gives the SE of the weighted mean of two independent estimates
      augmented_se <- 1 / sqrt(1/fia_only_se^2 + 1/nefin_se_5y^2)
    } else if (!is.na(fia_only_se) && fia_only_se > 0) {
      # Only FIA available
      augmented_se <- fia_only_se
    } else if (!is.na(nefin_se_5y) && nefin_se_5y > 0) {
      # Only NEFIN available
      augmented_se <- nefin_se_5y
    }
    
    # sensitivity of mean
    est_change <- augmented_mean - fia_only_mean
    est_change_pct <- if (!is.na(fia_only_mean) && fia_only_mean != 0) {
      100 * est_change / fia_only_mean
    } else {
      NA_real_
    }
    
    # uncertainty improvement (positive = NEFIN reduced uncertainty)
    uncertainty_improvement <- fia_only_se - augmented_se
    relative_uncertainty_drop <- if (!is.na(fia_only_se) && fia_only_se > 0) {
      uncertainty_improvement / fia_only_se
    } else {
      NA_real_
    }
    
    list(
      augmented_mean             = augmented_mean,
      augmented_se               = augmented_se,
      w_fia                      = w_fia,
      w_nefin                    = w_nefin,
      est_change                 = est_change,
      est_change_pct             = est_change_pct,
      uncertainty_improvement    = uncertainty_improvement,
      relative_uncertainty_drop  = relative_uncertainty_drop
    )
  }
  
  joined_5y <- joined_5y %>%
    mutate(
      calc = purrr::pmap(
        list(
          fia_only_mean, fia_only_se, fia_only_n,
          nefin_mean_5y, nefin_se_5y, nefin_n_5y
        ),
        safe_row_calc
      )
    ) %>%
    mutate(
      augmented_mean             = purrr::map_dbl(calc, `[[`, "augmented_mean"),
      augmented_se               = purrr::map_dbl(calc, `[[`, "augmented_se"),
      w_fia                      = purrr::map_dbl(calc, `[[`, "w_fia"),
      w_nefin                    = purrr::map_dbl(calc, `[[`, "w_nefin"),
      est_change                 = purrr::map_dbl(calc, `[[`, "est_change"),
      est_change_pct             = purrr::map_dbl(calc, `[[`, "est_change_pct"),
      uncertainty_improvement    = purrr::map_dbl(calc, `[[`, "uncertainty_improvement"),
      relative_uncertainty_drop  = purrr::map_dbl(calc, `[[`, "relative_uncertainty_drop")
    ) %>%
    select(-calc)
  
  # fuzz pressure: how much FIA uncertainty is coming from positional fuzz?
  joined_5y <- joined_5y %>%
    mutate(
      fuzz_pressure = dplyr::case_when(
        # Can't calculate if either value is NA
        is.na(positional_sd_5y) | is.na(fia_only_se) ~ NA_real_,
        # If SE is too small (including 0), ratio is undefined/infinite
        fia_only_se < min_se_for_ratio ~ NA_real_,
        # Normal case: calculate ratio
        TRUE ~ positional_sd_5y / fia_only_se
      )
    )
  
  # Standardize grid_name (fia → 64kha)
  grid_name_display <- standardize_scale_name(grid_name)
  
  # per-hex recommendation
  joined_5y <- joined_5y %>%
    mutate(
      recommendation = dplyr::case_when(
        # FIA already strong, NEFIN doesn't help much
        w_fia >= 0.99 &
          (is.na(relative_uncertainty_drop) | relative_uncertainty_drop <= se_drop_small) ~
          "USE_FIA_ONLY",
        
        # NEFIN is filling gaps where FIA is weak and actually lowers uncertainty,
        # and it's not obviously dominated by one source.
        w_fia < 0.67 &
          relative_uncertainty_drop > se_drop_good &
          (
            (is.na(max_source_frac) | max_source_frac < bias_concentration_threshold) |
              (!is.na(n_nefin_sources) & n_nefin_sources > 1)
          ) ~
          "USE_AUGMENTED",
        
        # Looks like NEFIN is dragging estimate around but might be biased
        w_fia < 0.99 &
          !is.na(max_source_frac) & max_source_frac >= bias_concentration_threshold &
          abs(est_change) > big_shift_threshold ~
          "NEEDS_REVIEW_FOR_BIAS",
        
        TRUE ~ "MIXED_SIGNAL"
      ),
      # Use standardized scale name (fia → 64kha) in output
      scale_name = grid_name_display,
      window_years = paste0(min(focal_years), "-", max(focal_years)),
      mc_reps = mc_reps,
      jitter_radius_m = cfg$jitter_radius_m,
      level_window = level_window
    )
  
  # plots for this grid - filter out NA values for plotting
  plot_data <- joined_5y %>%
    filter(!is.na(fuzz_pressure), !is.na(relative_uncertainty_drop),
           is.finite(fuzz_pressure), is.finite(relative_uncertainty_drop))
  
  if (nrow(plot_data) > 0) {
    p_fuzz_vs_drop <- ggplot(plot_data,
                             aes(x = fuzz_pressure,
                                 y = relative_uncertainty_drop,
                                 color = recommendation)) +
      geom_point(alpha = 0.6) +
      labs(
        title = paste0("Fuzz sensitivity vs NEFIN benefit (", grid_name, ", ", min(focal_years), "–", max(focal_years), ")"),
        subtitle = paste0("FIA positional uncertainty modeled with ", mc_reps, " jitter draws at ±", cfg$jitter_radius_m, " m"),
        x = "Fuzz pressure (positional_sd_5y / FIA-only SE)",
        y = "Relative SE reduction after adding NEFIN"
      ) +
      theme_minimal(base_size = 14)
    ggsave(
      fs::path(out_dir, paste0("fuzz_vs_uncertainty_drop_", grid_name, ".png")),
      p_fuzz_vs_drop,
      width = 8, height = 6, dpi = 300
    )
  }
  
  # Filter for uncertainty improvement plot
  plot_data_unc <- joined_5y %>%
    filter(!is.na(uncertainty_improvement), is.finite(uncertainty_improvement))
  
  if (nrow(plot_data_unc) > 0) {
    p_unc <- ggplot(plot_data_unc, aes(x = uncertainty_improvement)) +
      geom_histogram(bins = 30, fill = "darkgreen", color = "white") +
      geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
      labs(
        title = paste0("Uncertainty Change With NEFIN (", grid_name, ", ", min(focal_years), "–", max(focal_years), ")"),
        subtitle = "Positive = NEFIN lowered SE (good)",
        x = "FIA-only SE − Augmented SE (Mg/ha)",
        y = "Hex count"
      ) +
      theme_minimal(base_size = 14)
    ggsave(
      fs::path(out_dir, paste0("uncertainty_improvement_hist_5y_", grid_name, ".png")),
      p_unc,
      width = 8, height = 6, dpi = 300
    )
  }
  
  # Filter for estimate change plot
  plot_data_change <- joined_5y %>%
    filter(!is.na(est_change), is.finite(est_change))
  
  if (nrow(plot_data_change) > 0) {
    p_change <- ggplot(plot_data_change, aes(x = est_change)) +
      geom_histogram(bins = 30, fill = "steelblue", color = "white") +
      geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
      labs(
        title = paste0("Shift in Biomass After Adding NEFIN (", grid_name, ", ", min(focal_years), "–", max(focal_years), ")"),
        subtitle = "augmented_mean − FIA_only_mean (Mg/ha)",
        x = "Change in estimate (Mg/ha)",
        y = "Hex count"
      ) +
      theme_minimal(base_size = 14)
    ggsave(
      fs::path(out_dir, paste0("est_change_hist_5y_", grid_name, ".png")),
      p_change,
      width = 8, height = 6, dpi = 300
    )
  }
  
  # write per-hex CSV
  write_csv(
    joined_5y,
    fs::path(out_dir, paste0("summary_5y_", grid_name, ".csv"))
  )
  
  # guidance text per scale - handle NA in median calculations
  guidance_df <- joined_5y %>%
    group_by(recommendation) %>%
    summarise(
      n_hex                 = dplyr::n(),
      median_fuzz_pressure  = median(fuzz_pressure, na.rm = TRUE),
      median_rel_unc_drop   = median(relative_uncertainty_drop, na.rm = TRUE),
      median_est_change     = median(est_change, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(frac_hex = n_hex / sum(n_hex))
  
  # Format guidance with NA handling
  guidance_lines <- c(
    paste0("Scale/grid: ", grid_name),
    paste0("Years pooled: ", min(focal_years), "-", max(focal_years),
           " (window size ", level_window, "y)"),
    paste0("FIA fuzz modeled with ", mc_reps, " Monte Carlo draws at ±",
           cfg$jitter_radius_m, " m"),
    "",
    "Recommendation categories across hexes:",
    paste(
      sprintf(
        "- %s: %.1f%% of hexes | median fuzz_pressure=%s | median rel_uncertainty_drop=%s | median est_change=%s Mg/ha",
        guidance_df$recommendation,
        100 * guidance_df$frac_hex,
        ifelse(is.na(guidance_df$median_fuzz_pressure), "NA", sprintf("%.2f", guidance_df$median_fuzz_pressure)),
        ifelse(is.na(guidance_df$median_rel_unc_drop), "NA", sprintf("%.2f", guidance_df$median_rel_unc_drop)),
        ifelse(is.na(guidance_df$median_est_change), "NA", sprintf("%.2f", guidance_df$median_est_change))
      ),
      collapse = "\n"
    ),
    "",
    "USE_FIA_ONLY: FIA had adequate coverage; NEFIN did not materially improve uncertainty.",
    "USE_AUGMENTED: FIA coverage weak AND NEFIN reduced uncertainty without obvious single-source dominance.",
    "NEEDS_REVIEW_FOR_BIAS: NEFIN pulled estimates noticeably, but is dominated by one source in that hex (likely not landscape-representative).",
    "MIXED_SIGNAL: No clear call; consider reporting FIA but footnoting fuzz pressure."
  )
  
  writeLines(
    guidance_lines,
    fs::path(out_dir, paste0("scale_guidance_", grid_name, ".txt"))
  )
  
  joined_5y
}

# ------------------------------------------------------------
# driver over all grids from cfg$hex_grids
# ------------------------------------------------------------
run_all_scales <- function(cfg_path = "configs/process.yml",
                           nefin_assign_path = "data/processed/nefin_hex_assignments.csv") {
  
  cfg <- load_config(cfg_path)
  
  message("\nLoaded config:")
  message("  years: ", paste(cfg$years, collapse=", "))
  message("  level_window: ", cfg$level_window, "y")
  message("  mc_reps: ", cfg$mc_reps)
  message("  jitter_radius_m: ", cfg$jitter_radius_m, " m")
  
  # load NEFIN once
  nefin_df <- read_csv(nefin_assign_path, show_col_types = FALSE)
  nefin_cols <- colnames(nefin_df)
  
  results <- list()
  
  for (grid in cfg$hex_grids) {
    grid_name <- grid$name
    message("\n═══════════════════════════════════════════════════════════")
    message("Grid/scale: ", grid_name)
    message("Hex path:   ", grid$path)
    message("═══════════════════════════════════════════════════════════")
    
    # figure out which NEFIN column to use for this grid
    this_hex_col <- hex_col_for_grid(grid_name, nefin_cols)
    message("Using NEFIN column: ", this_hex_col)
    
    # find FIA results dir for this grid using cfg$metric and cfg$level_window
    fia_info <- find_fia_results_for_grid(
      grid_name = grid_name,
      metric = cfg$metric,
      level_window = cfg$level_window
    )
    message("FIA dir: ", fia_info$fia_dir)
    message("FIA csv: ", fia_info$fia_results_path)
    
    fia_df <- read_csv(fia_info$fia_results_path, show_col_types = FALSE)
    
    # make output dir for this scale
    out_dir <- fs::path("runs", paste0("fia_nefin_comparison_", grid_name))
    fs::dir_create(out_dir, recurse = TRUE)
    
    # 1) FIA stability across time
    stability_df <- summarize_fia_stability(fia_df, grid_name, out_dir)
    
    # 2) 5-year pooled comparison with recommendations, fuzz, bias
    summary_5y_df <- compare_5y_window(
      fia_df = fia_df,
      nefin_df = nefin_df,
      cfg = cfg,
      grid_name = grid_name,
      hex_col = this_hex_col,
      out_dir = out_dir
    )
    
    results[[grid_name]] <- list(
      stability = stability_df,
      summary_5y = summary_5y_df,
      out_dir = out_dir
    )
    
    message("\n✓ Finished grid: ", grid_name)
    message("  Output dir: ", out_dir)
  }
  
  invisible(results)
}

# ------------------------------------------------------------
# CLI entry point
# ------------------------------------------------------------
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  run_all_scales()
}
