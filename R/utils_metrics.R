# R/utils_metrics.R
# Metric computation utilities
# FIXED: Proper NA handling for single-observation hexes

suppressPackageStartupMessages({
  library(dplyr); library(readr)
})

acres_to_ha <- 0.40468564224
lbs_to_Mg   <- 0.00045359237

# Weighted mean and SE
# FIXED: Returns NA for SE when n < 2 (can't estimate variance from single observation)
w_mean_se <- function(y, w) {
  keep <- is.finite(y) & is.finite(w) & w > 0
  y <- y[keep]; w <- w[keep]
  if (!length(y)) return(list(mean=NA_real_, se=NA_real_, n=0L, sumw=0))
  
  p <- w/sum(w)
  mu <- sum(p*y)
  
  # FIX: SE is undefined (NA) with only 1 observation
  # You cannot estimate variance from a single data point

  if (length(y) < 2) {
    return(list(mean=mu, se=NA_real_, n=length(y), sumw=sum(w)))
  }
  
  n_eff <- (sum(w)^2) / sum(w^2)
  var_w <- sum(p * (y - mu)^2)
  
  # Apply Bessel correction for weighted variance (n/(n-1))
  var_w_corrected <- var_w * length(y) / (length(y) - 1)
  se <- sqrt(var_w_corrected / n_eff)
  
  list(mean=mu, se=se, n=length(y), sumw=sum(w))
}

# Build plot-level AGLB
build_plot_aglb <- function(tree_csv, plot_csv) {
  stopifnot(fs::file_exists(tree_csv), fs::file_exists(plot_csv))
  
  tr <- suppressMessages(readr::read_csv(tree_csv, guess_max = 1e6, show_col_types = FALSE))
  pl <- suppressMessages(readr::read_csv(plot_csv, guess_max = 1e6, show_col_types = FALSE))
  
  source("R/utils_spatial.R")
  pl <- normalize_plot_coords(pl)
  
  agg <- tr |>
    dplyr::mutate(bio_lb_per_ac = DRYBIO_AG * TPA_UNADJ) |>
    dplyr::group_by(STATECD, PLT_CN) |>
    dplyr::summarise(bio_lb_per_ac = sum(bio_lb_per_ac, na.rm=TRUE), .groups="drop")
  
  df <- pl |>
    dplyr::select(CN, STATECD, MEASYEAR, lat_public, lon_public, UNITCD, COUNTYCD, PLOT) |>
    dplyr::left_join(agg, by = c("STATECD" = "STATECD", "CN" = "PLT_CN")) |>
    dplyr::mutate(
      aglb_Mg_per_ac = bio_lb_per_ac * 0.00045359237,
      aglb_Mg_per_ha = aglb_Mg_per_ac / 0.40468564224
    )
  
  df
}

# Aggregate to hex in time window
aggregate_hex <- function(df, hex_col = "hex_id", year_label, window_years = 3, w_col = NULL, y_col) {
  if (!is.data.frame(df)) stop("aggregate_hex: df must be a data.frame")
  
  years <- (year_label - floor(window_years/2)):(year_label + floor(window_years/2))
  dfw <- dplyr::filter(df, .data$MEASYEAR %in% years)
  if (!nrow(dfw)) return(tibble::tibble())
  
  dfw$w <- if (!is.null(w_col) && (w_col %in% names(dfw))) dfw[[w_col]] else 1
  
  key <- as.character(dfw[[hex_col]])
  out <- dfw |>
    dplyr::mutate(.key = key) |>
    dplyr::group_by(.key) |>
    dplyr::reframe({
      m <- w_mean_se(.data[[y_col]], w = w)
      tibble::tibble(mean = m$mean, se = m$se, n_plots = m$n, sum_weight = m$sumw)
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      hex_id = as.character(.key),
      year_label = year_label,
      window = paste0(window_years, "y")
    ) |>
    dplyr::select(-.key)
  
  out
}
