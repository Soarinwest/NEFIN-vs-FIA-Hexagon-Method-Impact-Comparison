#!/usr/bin/env Rscript
# process_to_hex.R — Build hex-level products from FIA (and optional NEFIN),
# including FIA positional Monte Carlo; state-safe joins.

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(purrr)
  library(stringr); library(yaml); library(fs); library(glue); library(tidyr)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b
acres_to_ha <- 0.40468564224
lbs_to_Mg   <- 0.00045359237

read_config <- function(path="configs/process.yml") {
  if (fs::file_exists(path)) yaml::read_yaml(path) else list()
}
ensure_dirs <- function(paths) purrr::walk(paths, ~ if (!fs::dir_exists(.x)) fs::dir_create(.x, recurse=TRUE))

# --- helpers to normalize coord columns from PLOT --------------------------------
find_col <- function(df, candidates) {
  # return the first matching column name (case-insensitive), else NA
  up <- toupper(names(df))
  for (cand in toupper(candidates)) {
    hit <- which(up == cand)
    if (length(hit)) return(names(df)[hit[1]])
  }
  NA_character_
}

normalize_plot_coords <- function(pl) {
  # common variants we’ve seen across FIADB drops
  lat_candidates <- c("lat_public","LAT_PUBLIC","LAT","LATITUDE","DEC_LAT","DEC_LATITUDE","Y")
  lon_candidates <- c("lon_public","LON_PUBLIC","LON","LONGITUDE","DEC_LON","DEC_LONGITUDE","X")
  
  lat_col <- find_col(pl, lat_candidates)
  lon_col <- find_col(pl, lon_candidates)
  
  if (is.na(lat_col) || is.na(lon_col))
    stop("Could not find latitude/longitude columns in PLOT; looked for: ",
         paste(lat_candidates, collapse=", "), " / ",
         paste(lon_candidates, collapse=", "))
  
  pl <- pl |>
    dplyr::mutate(
      lat_public = suppressWarnings(as.numeric(.data[[lat_col]])),
      lon_public = suppressWarnings(as.numeric(.data[[lon_col]]))
    )
  
  # basic sanity: drop impossible values
  pl <- pl |>
    dplyr::mutate(
      lat_public = dplyr::if_else(lat_public >= -90 & lat_public <= 90, lat_public, NA_real_),
      lon_public = dplyr::if_else(lon_public >= -180 & lon_public <= 180, lon_public, NA_real_)
    )
  
  pl
}

# Fast, practical weighted mean + SE using effective sample size
w_mean_se <- function(y, w) {
  keep <- is.finite(y) & is.finite(w) & w > 0
  y <- y[keep]; w <- w[keep]
  if (!length(y)) return(list(mean=NA_real_, se=NA_real_, n=0L, sumw=0))
  p <- w/sum(w); mu <- sum(p*y)
  n_eff <- (sum(w)^2) / sum(w^2)
  var_w <- sum(p * (y - mu)^2)
  se <- sqrt(var_w / max(n_eff, 1))
  list(mean=mu, se=se, n=length(y), sumw=sum(w))
}

# Builds plot-level AGLB (Mg/ha); joins are STATE-safe and coords are normalized.
build_plot_aglb <- function(tree_csv, plot_csv) {
  stopifnot(fs::file_exists(tree_csv), fs::file_exists(plot_csv))
  tr <- suppressMessages(readr::read_csv(tree_csv, guess_max = 1e6, show_col_types = FALSE))
  pl <- suppressMessages(readr::read_csv(plot_csv, guess_max = 1e6, show_col_types = FALSE))
  
  need_tree <- c("PLT_CN","STATECD","DRYBIO_AG","TPA_UNADJ")
  need_plot <- c("CN","STATECD","MEASYEAR","UNITCD","COUNTYCD","PLOT")
  if (!all(need_tree %in% names(tr))) stop("TREE missing: ", paste(setdiff(need_tree, names(tr)), collapse=", "))
  if (!all(need_plot %in% names(pl))) stop("PLOT missing: ", paste(setdiff(need_plot, names(pl)), collapse=", "))
  
  # normalize coord columns to lat_public / lon_public (handles LAT/LON, etc.)
  pl <- normalize_plot_coords(pl)
  
  agg <- tr |>
    dplyr::mutate(bio_lb_per_ac = DRYBIO_AG * TPA_UNADJ) |>
    dplyr::group_by(STATECD, PLT_CN) |>
    dplyr::summarise(bio_lb_per_ac = sum(bio_lb_per_ac, na.rm=TRUE), .groups="drop")
  
  df <- pl |>
    dplyr::select(CN, STATECD, MEASYEAR, lat_public, lon_public, UNITCD, COUNTYCD, PLOT) |>
    dplyr::left_join(agg, by = c("STATECD" = "STATECD", "CN" = "PLT_CN")) |>
    dplyr::mutate(
      aglb_Mg_per_ac = bio_lb_per_ac * 0.00045359237,   # lb → Mg
      aglb_Mg_per_ha = aglb_Mg_per_ac / 0.40468564224   # ac → ha
    )
  
  df
}

# Assign plots to hex grid (expects a hex_id field). Reprojects to EPSG:5070 internally.
assign_plots_to_hex <- function(df_plots, hex_path, hex_layer = NULL) {
  if (!fs::file_exists(hex_path)) stop("Hex grid not found: ", hex_path)
  hx <- if (is.null(hex_layer)) sf::st_read(hex_path, quiet = TRUE) else sf::st_read(hex_path, layer = hex_layer, quiet = TRUE)
  if (!("hex_id" %in% names(hx))) { hx$hex_id <- seq_len(nrow(hx)); warning("hex_id missing; created sequential IDs") }
  hx <- st_make_valid(hx); hx_5070 <- st_transform(hx, 5070)
  pts <- st_as_sf(df_plots, coords=c("lon_public","lat_public"), crs=4326, remove=FALSE) |>
    st_transform(5070)
  j <- st_join(pts, hx_5070["hex_id"], left=TRUE)
  st_drop_geometry(j) |>
    as.data.frame()
}

# Aggregate to hex in a centered time window on year_label
aggregate_hex <- function(df, hex_col = "hex_id", year_label, window_years = 3, w_col = NULL, y_col) {
  years <- (year_label - floor(window_years/2)):(year_label + floor(window_years/2))
  dfw <- if ("MEASYEAR" %in% names(df)) dplyr::filter(df, .data$MEASYEAR %in% years) else dplyr::filter(df, .data$mid %in% years)
  if (!nrow(dfw)) return(tibble())
  dfw$w <- if (!is.null(w_col) && (w_col %in% names(dfw))) dfw[[w_col]] else 1
  df_group <- dfw |>
    group_by(.data[[hex_col]]) |>
    reframe({
      m <- w_mean_se(.data[[y_col]], w = w)
      tibble(mean = m$mean, se = m$se, n_plots = m$n, sum_weight = m$sumw)
    }) |>
    ungroup() |>
    mutate(hex_id = .data[[hex_col]], year_label = year_label, window = paste0(window_years, "y"))
  df_group
}

# Monte Carlo positional SD for FIA (jitter public coords)
positional_mc <- function(df_points, hex_path, hex_layer = NULL, year_label, window_years = 3, R = 100, radius_m = 1609.34) {
  if (!nrow(df_points)) return(list(replicates=tibble(), summary=tibble()))
  hx <- if (is.null(hex_layer)) sf::st_read(hex_path, quiet = TRUE) else sf::st_read(hex_path, layer = hex_layer, quiet = TRUE)
  hx_5070 <- sf::st_transform(sf::st_make_valid(hx), 5070)
  pts <- sf::st_as_sf(df_points, coords=c("lon_public","lat_public"), crs=4326, remove=FALSE) |>
    sf::st_transform(5070)
  res <- vector("list", R)
  years <- (year_label - floor(window_years/2)):(year_label + floor(window_years/2))
  for (r in seq_len(R)) {
    u <- runif(nrow(pts)); rr <- sqrt(u)*radius_m; theta <- runif(nrow(pts),0,2*pi)
    dx <- rr*cos(theta); dy <- rr*sin(theta)
    pts_j <- sf::st_set_geometry(pts, sf::st_geometry(pts) + cbind(dx,dy))
    j <- sf::st_join(pts_j, hx_5070["hex_id"], left=TRUE) |> sf::st_drop_geometry()
    g <- j |>
      dplyr::filter(.data$MEASYEAR %in% years) |>
      dplyr::group_by(hex_id) |>
      dplyr::summarise(mean_rep = mean(aglb_Mg_per_ha, na.rm=TRUE), .groups="drop")
    g$replicate_id <- r
    res[[r]] <- g
  }
  all <- dplyr::bind_rows(res)
  pos <- all |>
    dplyr::group_by(hex_id) |>
    dplyr::summarise(positional_sd = stats::sd(mean_rep, na.rm=TRUE), .groups="drop") |>
    dplyr::mutate(year_label = year_label, window = paste0(window_years, "y"))
  list(replicates = all, summary = pos)
}

process_to_hex <- function(project_dir = ".",
                           hex_path = "data/hex/hex_grid.geojson",
                           hex_layer = NULL,
                           metric = "aglb",
                           years = 2018:2020,
                           level_window = 3,
                           rate_window = 5,
                           mc_reps = 100,
                           jitter_radius_m = 1609.34,
                           run_id = NULL) {
  
  # Prefer region-wide slices if present (from state-by-state pull)
  fia_root <- if (fs::dir_exists(fs::path(project_dir, "data", "interim", "fia_region"))) {
    fs::path(project_dir, "data", "interim", "fia_region")
  } else {
    fs::path(project_dir, "data", "interim", "fia")
  }
  
  paths <- list(
    fia_dir  = fia_root,
    nefin_dir = fs::path(project_dir, "data", "interim", "nefin"),
    runs_dir  = fs::path(project_dir, "runs")
  )
  ensure_dirs(paths)
  
  # 1) FIA plot-level metric
  fia_plot <- build_plot_aglb(
    tree_csv = fs::path(paths$fia_dir, "tree.csv"),
    plot_csv = fs::path(paths$fia_dir, "plot.csv")
  ) |> dplyr::mutate(source = "fia")
  
  fia_plot_hex <- assign_plots_to_hex(fia_plot, hex_path, hex_layer)
  
  # 2) NEFIN optional
  nefin_plot_path <- fs::path(paths$nefin_dir, "nefin_plot.csv")
  if (fs::file_exists(nefin_plot_path)) {
    ne <- readr::read_csv(nefin_plot_path, show_col_types = FALSE)
    need <- c("plot_id","visit_year","lon","lat","aglb_Mg_per_ha")
    if (!all(need %in% names(ne))) stop("nefin_plot.csv missing columns: ", paste(need, collapse=", "))
    ne <- ne |>
      dplyr::transmute(
        CN = NA_integer_, STATECD = NA_integer_,
        MEASYEAR = visit_year, lat_public = lat, lon_public = lon,
        aglb_Mg_per_ha = aglb_Mg_per_ha,
        UNITCD = NA_integer_, COUNTYCD = NA_integer_, PLOT = plot_id,
        source = "nefin"
      )
    ne_hex <- assign_plots_to_hex(ne, hex_path, hex_layer)
  } else {
    ne_hex <- tibble::tibble()
  }
  
  out_design <- vector("list", length(years))
  out_pos_sd <- vector("list", length(years))
  out_reps <- vector("list", length(years))
  
  for (i in seq_along(years)) {
    yr <- years[i]
    message(glue("=== Year {yr} ==="))
    fia_hex_stats <- aggregate_hex(fia_plot_hex, year_label = yr, window_years = level_window, y_col = "aglb_Mg_per_ha") |>
      dplyr::mutate(source = "fia")
    ne_hex_stats <- if (nrow(ne_hex)) aggregate_hex(ne_hex, year_label = yr, window_years = level_window, y_col = "aglb_Mg_per_ha") |> dplyr::mutate(source="nefin") else tibble::tibble()
    fused_hex_stats <- dplyr::bind_rows(
      fia_plot_hex |> dplyr::mutate(source="fia"),
      ne_hex |> dplyr::mutate(source="nefin")
    ) |>
      aggregate_hex(year_label = yr, window_years = level_window, y_col = "aglb_Mg_per_ha") |>
      dplyr::mutate(source = "fused")
    mc <- positional_mc(fia_plot_hex, hex_path, hex_layer, year_label = yr, window_years = level_window,
                        R = mc_reps, radius_m = jitter_radius_m)
    out_pos_sd[[i]] <- mc$summary |> dplyr::mutate(metric = metric)
    out_reps[[i]]   <- mc$replicates |> dplyr::mutate(metric = metric, year_label = yr, window = paste0(level_window,"y"))
    out_design[[i]] <- dplyr::bind_rows(fia_hex_stats, ne_hex_stats, fused_hex_stats) |> dplyr::mutate(metric = metric)
  }
  
  design <- dplyr::bind_rows(out_design)
  pos_sd <- dplyr::bind_rows(out_pos_sd)
  reps   <- dplyr::bind_rows(out_reps)
  
  design_wide <- design |>
    dplyr::select(hex_id, year_label, window, source, metric, mean, se, n_plots, sum_weight) |>
    tidyr::pivot_wider(names_from = source, values_from = c(mean, se, n_plots), names_sep = "_")
  
  joined <- design_wide |>
    dplyr::left_join(pos_sd |> dplyr::select(hex_id, year_label, window, positional_sd), by = c("hex_id","year_label","window")) |>
    dplyr::mutate(
      total_sd_fia = sqrt((se_fia %||% NA_real_)^2 + (positional_sd %||% 0)^2),
      se_ratio = dplyr::if_else(is.finite(se_fia) & is.finite(se_fused) & se_fused > 0, se_fia / se_fused, NA_real_)
    )
  
  run_id <- run_id %||% paste0(format(Sys.Date(), "%Y-%m-%d"), "_", metric, "_W", level_window, "y")
  out_dir <- fs::path(paths$runs_dir, run_id); if (!fs::dir_exists(out_dir)) fs::dir_create(out_dir, recurse = TRUE)
  
  readr::write_csv(joined, fs::path(out_dir, "hex_joined.csv"))
  readr::write_csv(reps,   fs::path(out_dir, "fia_mc_replicates.csv"))
  
  message("✔ Wrote: ", fs::path(out_dir, "hex_joined.csv"))
  message("✔ Wrote: ", fs::path(out_dir, "fia_mc_replicates.csv"))
  invisible(list(hex_joined = fs::path(out_dir, "hex_joined.csv"),
                 fia_mc_replicates = fs::path(out_dir, "fia_mc_replicates.csv"),
                 run_id = run_id))
}

# CLI (only runs when called directly, not when sourced)
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  cfg <- read_config("configs/process.yml")
  process_to_hex(
    project_dir = cfg$project_dir %||% ".",
    hex_path = cfg$hex_path %||% "data/hex/hex_grid.geojson",
    hex_layer = cfg$hex_layer %||% NULL,
    metric = cfg$metric %||% "aglb",
    years = cfg$years %||% 2018:2020,
    level_window = cfg$level_window %||% 3,
    rate_window = cfg$rate_window %||% 5,
    mc_reps = cfg$mc_reps %||% 100,
    jitter_radius_m = cfg$jitter_radius_m %||% 1609.34,
    run_id = cfg$run_id %||% NULL
  )
}
