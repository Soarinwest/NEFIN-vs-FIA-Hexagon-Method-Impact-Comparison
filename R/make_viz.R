#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(ggplot2)
})

# Robust reader for run outputs (CSV or Parquet)
read_run_tbl <- function(path_base) {
  csv <- paste0(path_base, ".csv")
  if (file.exists(csv)) return(readr::read_csv(csv, show_col_types = FALSE))
  pq  <- paste0(path_base, ".parquet")
  if (file.exists(pq)) {
    if (!requireNamespace("arrow", quietly = TRUE)) stop("Install 'arrow' to read Parquet.")
    return(arrow::read_parquet(pq))
  }
  stop("Missing run output: ", csv, " / ", pq)
}

# Make a simple choropleth with quantile breaks
plot_hex_var <- function(hx, df, var, title, out_png) {
  brks <- stats::quantile(df[[var]], probs = c(0, .2, .4, .6, .8, 1), na.rm = TRUE)
  df$cut <- cut(df[[var]], include.lowest = TRUE, breaks = unique(brks))
  g <- hx |>
    left_join(df, by = "hex_id") |>
    ggplot() +
    geom_sf(aes(fill = cut), linewidth = 0.1, color = "grey30") +
    ggtitle(title) +
    guides(fill = guide_legend(title = var)) +
    theme_minimal()
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_png, g, width = 9, height = 7, dpi = 200)
}

# Main entry: make a few PNGs per year
make_run_viz <- function(run_dir, hex_path, years = NULL) {
  stopifnot(dir.exists(run_dir), file.exists(hex_path))
  hx <- sf::st_read(hex_path, quiet = TRUE)
  if (!"hex_id" %in% names(hx)) {
    if ("ID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = ID)
    else if ("OBJECTID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = OBJECTID)
    else hx$hex_id <- seq_len(nrow(hx))
  }
  viz_dir <- file.path(run_dir, "viz"); dir.create(viz_dir, showWarnings = FALSE)
  
  # If years not supplied, infer from files like hex_year_YYYY.*
  if (is.null(years)) {
    yfiles <- list.files(run_dir, pattern = "^hex_year_\\d{4}\\.(csv|parquet)$", full.names = FALSE)
    years <- sort(unique(as.integer(gsub("^hex_year_(\\d{4})\\..*$", "\\1", yfiles))))
  }
  
  # n_plots + aglb_mean per year
  for (yr in years) {
    yr_base <- file.path(run_dir, sprintf("hex_year_%d", yr))
    hy <- read_run_tbl(yr_base)
    hy$n_plots <- hy$n_plots %||% NA_integer_
    
    plot_hex_var(hx, hy, "n_plots",
                 sprintf("FIA plots per hex (%d)", yr),
                 file.path(viz_dir, sprintf("plots_per_hex_%d.png", yr)))
    
    if ("aglb_mean" %in% names(hy)) {
      plot_hex_var(hx, hy, "aglb_mean",
                   sprintf("Aboveground live biomass (Mg/ha) — mean (%d)", yr),
                   file.path(viz_dir, sprintf("aglb_mean_%d.png", yr)))
    }
    
    # Positional SD (from MC), if present
    ps_base <- file.path(run_dir, sprintf("hex_positional_sd_%d", yr))
    if (file.exists(paste0(ps_base, ".csv")) || file.exists(paste0(ps_base, ".parquet"))) {
      ps <- read_run_tbl(ps_base)
      if ("positional_sd" %in% names(ps)) {
        plot_hex_var(hx, ps, "positional_sd",
                     sprintf("Positional SD from MC (%d)", yr),
                     file.path(viz_dir, sprintf("positional_sd_%d.png", yr)))
        # quick histogram too
        h <- ggplot(ps, aes(x = positional_sd)) + geom_histogram(bins = 30) +
          ggtitle(sprintf("Distribution of positional SD (%d)", yr)) + theme_minimal()
        ggsave(file.path(viz_dir, sprintf("positional_sd_hist_%d.png", yr)), h, width = 8, height = 5, dpi = 200)
      }
    }
  }
  
  message("✔ Visuals written to: ", viz_dir)
}

`%||%` <- function(a,b) if (!is.null(a)) a else b