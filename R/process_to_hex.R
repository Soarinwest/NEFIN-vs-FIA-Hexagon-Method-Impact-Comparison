# process_to_hex.R 
### Foreword ---------------------------------------------------------------------------------
### Author: Soren Donisvitch
### Date: 10/02/2025
### Dependents: R (>= 3.5), DBI, RSQLite, fs, curl, readr, dplyr, glue, withr, jsonlite, rlang
### Foreword: The use or application of these code without permission of the author is prohibited.The author is not ->
###           liable for the use, modification, or any other application of this or other provided scripts.
### Build hex-level products from FIA (and optional NEFIN), including FIA positional Monte Carlo; state-safe joins.

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

# --- helpers --------------------------------
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

write_tbl <- function(df, path, fmt = "parquet", compression = "zstd") {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  fmt <- tolower(fmt)
  if (fmt == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) stop("Install 'arrow'")
    arrow::write_parquet(df, path, compression = compression)
  } else {
    readr::write_csv(df, path)
  }
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

build_hex_union_5070 <- function(hex_path, hex_layer = NULL) {
  hx <- if (is.null(hex_layer)) sf::st_read(hex_path, quiet = TRUE) else sf::st_read(hex_path, layer = hex_layer, quiet = TRUE)
  if (!("hex_id" %in% names(hx))) {
    if ("ID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = ID)
    else if ("OBJECTID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = OBJECTID)
    else hx$hex_id <- seq_len(nrow(hx))
  }
  hx <- sf::st_make_valid(sf::st_zm(hx, drop = TRUE))
  hx$hex_id <- as.character(hx$hex_id)
  sf::st_transform(hx, 5070) |>
    sf::st_union() |>
    sf::st_make_valid()
}

build_state_polys_5070 <- function(state_geo_path, state_field = "STATEFP") {
  st <- sf::st_read(state_geo_path, quiet = TRUE) |> sf::st_make_valid() |> sf::st_transform(5070)
  if (!state_field %in% names(st)) stop("state_field not found in state layer: ", state_field)
  vals <- st[[state_field]]
  if (is.character(vals)) st$STATECD <- suppressWarnings(as.integer(vals)) else st$STATECD <- as.integer(vals)
  st[, c("STATECD","geometry")]
}

options(dplyr.summarise.inform = TRUE)

dbg_hex <- function(df, tag) {
  nm <- names(df)
  message(sprintf("[DBG] %s: %d rows, %d cols", tag, nrow(df), length(nm)))
  if ("hex_id" %in% nm) {
    cls <- paste(class(df$hex_id), collapse = "|")
    sam <- tryCatch(paste(head(unique(df$hex_id), 5), collapse = ", "), error = function(e) "<err>")
    nas <- sum(is.na(df$hex_id))
    message(sprintf("       hex_id: class=%s, NA=%d, sample={%s}", cls, nas, sam))
  } else {
    message("       hex_id: <MISSING>")
    cand <- intersect(c("hex_id.x","hex_id.y","ID","OBJECTID","id","HEX_ID","hexid"), nm)
    if (length(cand)) message("       candidates present: ", paste(cand, collapse=", "))
  }
  invisible(df)
}

ensure_hex_id <- function(df) {
  nm <- names(df)
  if (!("hex_id" %in% nm)) {
    cand <- intersect(c("hex_id.y","hex_id.x","ID","OBJECTID","id","HEX_ID","hexid"), nm)
    if (length(cand)) {
      df$hex_id <- df[[cand[1]]]
    } else {
      stop("ensure_hex_id: no hex_id nor candidates. Columns: ", paste(nm, collapse=", "))
    }
  }
  df$hex_id <- as.character(df$hex_id)
  df
}

constrained_jitter_once <- function(pts_5070, radius_m, hex_union_5070,
                                    state_polys_5070 = NULL, max_reroll = 20) {
  n <- nrow(pts_5070); if (!n) return(pts_5070)
  
  # Ensure all inputs have identical CRS by forcing them to share one CRS object
  crs_ref <- sf::st_crs(5070)
  sf::st_crs(pts_5070) <- crs_ref
  sf::st_crs(hex_union_5070) <- crs_ref
  if (!is.null(state_polys_5070)) sf::st_crs(state_polys_5070) <- crs_ref
  
  # first proposal
  u <- runif(n); r <- sqrt(u) * radius_m; theta <- runif(n, 0, 2*pi)
  dx <- r * cos(theta); dy <- r * sin(theta)
  moved <- sf::st_set_geometry(pts_5070, sf::st_geometry(pts_5070) + cbind(dx, dy))
  sf::st_crs(moved) <- crs_ref  # Ensure moved also has same CRS
  
  inside <- suppressWarnings(sf::st_within(moved, hex_union_5070, sparse = FALSE)[,1])
  if (!is.null(state_polys_5070) && "STATECD" %in% names(pts_5070)) {
    st_lut <- split(state_polys_5070, state_polys_5070$STATECD)
    chk <- logical(n)
    for (i in seq_len(n)) {
      poly <- st_lut[[as.character(pts_5070$STATECD[i])]]
      if (!is.null(poly)) {
        sf::st_crs(poly) <- crs_ref  # Ensure poly has same CRS
        chk[i] <- as.logical(sf::st_within(moved[i,], poly, sparse = FALSE)[,1])
      } else chk[i] <- TRUE
    }
    inside <- inside & chk
  }
  
  tries <- 1L
  while (any(!inside) && tries < max_reroll) {
    idx <- which(!inside)
    u <- runif(length(idx)); r <- sqrt(u) * radius_m; theta <- runif(length(idx), 0, 2*pi)
    dx[idx] <- r * cos(theta); dy[idx] <- r * sin(theta)
    moved[idx, ] <- sf::st_set_geometry(pts_5070[idx, ], sf::st_geometry(pts_5070[idx, ]) + cbind(dx[idx], dy[idx]))
    sf::st_crs(moved) <- crs_ref  # Reset CRS after geometry changes
    
    inside[idx] <- suppressWarnings(sf::st_within(moved[idx,], hex_union_5070, sparse = FALSE)[,1])
    if (!is.null(state_polys_5070) && "STATECD" %in% names(pts_5070)) {
      for (j in seq_along(idx)) {
        i <- idx[j]
        poly <- st_lut[[as.character(pts_5070$STATECD[i])]]
        if (!is.null(poly)) {
          sf::st_crs(poly) <- crs_ref
          inside[i] <- inside[i] & as.logical(sf::st_within(moved[i,], poly, sparse = FALSE)[,1])
        }
      }
    }
    tries <- tries + 1L
  }
  
  # snap any stragglers to nearest boundary
  if (any(!inside)) {
    bad <- which(!inside)
    for (i in bad) {
      poly <- if (!is.null(state_polys_5070) && "STATECD" %in% names(pts_5070)) {
        p <- st_lut[[as.character(pts_5070$STATECD[i])]] |> sf::st_intersection(hex_union_5070)
        sf::st_crs(p) <- crs_ref
        p
      } else {
        hex_union_5070
      }
      nearest <- sf::st_nearest_points(moved[i,], poly)
      moved[i,] <- sf::st_set_geometry(moved[i,], sf::st_cast(nearest, "POINT")[2])
    }
    sf::st_crs(moved) <- crs_ref  # Final CRS assignment
  }
  moved
}

assign_plots_to_hex <- function(df_plots, hex_path, hex_layer = NULL) {
  if (!fs::file_exists(hex_path)) stop("Hex grid not found: ", hex_path)
  
  old_s2 <- sf::sf_use_s2(); on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)  # GEOS planar; we’ll join in EPSG:5070
  
  # ---- Read hex and ensure a 'hex_id' column ----
  hx <- if (is.null(hex_layer)) sf::st_read(hex_path, quiet = TRUE) else sf::st_read(hex_path, layer = hex_layer, quiet = TRUE)
  hx <- sf::st_make_valid(sf::st_zm(hx, drop = TRUE, what = "ZM"))
  if (is.na(sf::st_crs(hx))) sf::st_crs(hx) <- sf::st_crs(4326)
  
  # Guarantee a hex_id attribute
  if (!("hex_id" %in% names(hx))) {
    if ("ID" %in% names(hx)) {
      hx <- dplyr::rename(hx, hex_id = ID)
    } else if ("OBJECTID" %in% names(hx)) {
      hx <- dplyr::rename(hx, hex_id = OBJECTID)
    } else {
      hx$hex_id <- seq_len(nrow(hx))
      warning("hex_id missing; created sequential IDs")
    }
  }
  hx$hex_id <- as.character(hx$hex_id)
  
  # ---- Build points from coords ----
  dfp <- df_plots |>
    dplyr::filter(is.finite(lat_public), is.finite(lon_public),
                  lat_public >= -90, lat_public <= 90,
                  lon_public >= -180, lon_public <= 180)
  if (!nrow(dfp)) {
    warning("No valid coordinates found in df_plots")
    return(dplyr::mutate(df_plots, hex_id = NA_character_))
  }
  
  pts <- sf::st_as_sf(dfp, coords = c("lon_public","lat_public"), crs = 4326, remove = FALSE)
  pts <- sf::st_make_valid(sf::st_zm(pts, drop = TRUE, what = "ZM"))
  
  # If points already carry a hex_id from an earlier pass, drop it to avoid .x/.y suffixing
  if ("hex_id" %in% names(pts)) pts$hex_id <- NULL
  
  # ---- Project both to 5070; drop CRS to bypass strict equality test ----
  crs5070 <- sf::st_crs(5070)
  hx_5070  <- sf::st_transform(hx,  crs5070)
  pts_5070 <- sf::st_transform(pts, crs5070)

  # ---- Join ----
  joined <- sf::st_join(pts_5070, hx_5070["hex_id"], left = TRUE, join = sf::st_intersects)
  out <- sf::st_drop_geometry(joined) |> as.data.frame()
  
  # ---- Normalize hex_id column name if suffixed or differently named ----
  if (!("hex_id" %in% names(out))) {
    cand <- intersect(c("hex_id.y","hex_id.x","ID","OBJECTID"), names(out))
    if (length(cand)) {
      out$hex_id <- out[[cand[1]]]
      out[cand] <- NULL
    } else {
      stop("Join ran but no hex_id column found in result; check hex layer attributes.")
    }
  }
  
  # quick summary for sanity
  miss <- sum(is.na(out$hex_id))
  message(sprintf("• assign_plots_to_hex: matched %s / %s points (%.1f%%), %s NA",
                  format(nrow(out) - miss, big.mark = ","),
                  format(nrow(out), big.mark = ","),
                  100 * (nrow(out) - miss) / nrow(out),
                  miss))
  # ensure canonical type
  out$hex_id <- as.character(out$hex_id)
  return(out)
}

# Aggregate to hex in a centered time window on year_label
aggregate_hex <- function(df, hex_col = "hex_id", year_label, window_years = 3, w_col = NULL, y_col) {
  if (!is.data.frame(df)) stop("aggregate_hex: df must be a data.frame")
  if (!(hex_col %in% names(df))) df <- ensure_hex_id(df)
  
  years <- (year_label - floor(window_years/2)):(year_label + floor(window_years/2))
  dfw <- if ("MEASYEAR" %in% names(df)) dplyr::filter(df, .data$MEASYEAR %in% years) else dplyr::filter(df, .data$mid %in% years)
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
# Monte Carlo positional SD for FIA (jitter public coords)
positional_mc <- function(df_points, hex_path, hex_layer = NULL,
                          year_label, window_years = 3,
                          R = 100, radius_m = 1609.34,
                          persist_replicates = FALSE, thin_every = 0L,
                          out_dir = NULL) {
  if (!nrow(df_points)) return(list(replicates = dplyr::tibble(), summary = dplyr::tibble()))
  
  cfg <- read_config("configs/process.yml")
  use_hex_union <- isTRUE(cfg$mask$use_hex_union)
  use_state     <- isTRUE(cfg$mask$use_state_constraint)
  
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  # Create single CRS reference object
  crs_ref <- sf::st_crs(5070)
  
  # Build hex union and state polys if needed
  hex_union_5070  <- if (use_hex_union) {
    hu <- build_hex_union_5070(hex_path, hex_layer)
    sf::st_crs(hu) <- crs_ref
    hu
  } else NULL
  
  state_polys_5070 <- if (use_state) {
    sp <- build_state_polys_5070(cfg$mask$state_geo_path %||% "", cfg$mask$state_field %||% "STATEFP")
    sf::st_crs(sp) <- crs_ref
    sp
  } else NULL
  
  old_s2 <- sf::sf_use_s2(); on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  # Hex
  hx <- if (is.null(hex_layer)) sf::st_read(hex_path, quiet = TRUE) else sf::st_read(hex_path, layer = hex_layer, quiet = TRUE)
  if (is.na(sf::st_crs(hx))) sf::st_crs(hx) <- sf::st_crs(4326)
  hx <- sf::st_make_valid(sf::st_zm(hx, drop = TRUE, what = "ZM"))
  
  if (!("hex_id" %in% names(hx))) {
    if ("ID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = ID)
    else if ("OBJECTID" %in% names(hx)) hx <- dplyr::rename(hx, hex_id = OBJECTID)
    else hx$hex_id <- seq_len(nrow(hx))
  }
  hx$hex_id <- as.character(hx$hex_id)
  
  # Points
  pts <- sf::st_as_sf(df_points, coords = c("lon_public","lat_public"), crs = 4326, remove = FALSE)
  pts <- sf::st_make_valid(sf::st_zm(pts, drop = TRUE, what = "ZM"))
  
  # Project both to 5070 and drop CRS to avoid equality check
  crs5070 <- sf::st_crs(5070)
  hx_5070  <- sf::st_transform(hx,  crs5070)
  pts_5070 <- sf::st_transform(pts, crs5070)

  
  if (thin_every > 0L && !is.null(out_dir)) {
    dir.create(file.path(out_dir, "replicates"), recursive = TRUE, showWarnings = FALSE)
  }
  
  years <- (year_label - floor(window_years/2)):(year_label + floor(window_years/2))
  
  # If not persisting, keep only streaming stats: n, sum, sumsq
  if (!persist_replicates) {
    acc <- dplyr::tibble(hex_id = character(), n = integer(), sum = double(), sumsq = double())
    for (r in seq_len(R)) {
      u <- runif(nrow(pts_5070)); rr <- sqrt(u) * radius_m; theta <- runif(nrow(pts_5070), 0, 2*pi)
      dx <- rr * cos(theta); dy <- rr * sin(theta)
      if (!is.null(hex_union_5070) || !is.null(state_polys_5070)) {
        pts_j <- constrained_jitter_once(pts_5070, radius_m, hex_union_5070, state_polys_5070, max_reroll = as.integer(cfg$mask$max_reroll %||% 20))
      } else {
        u <- runif(nrow(pts_5070)); rr <- sqrt(u) * radius_m; theta <- runif(nrow(pts_5070), 0, 2*pi)
        dx <- rr * cos(theta); dy <- rr * sin(theta)
        pts_j <- sf::st_set_geometry(pts_5070, sf::st_geometry(pts_5070) + cbind(dx, dy))
      }
      
      j <- sf::st_join(pts_j, hx_5070["hex_id"], left = TRUE, join = sf::st_intersects) |> sf::st_drop_geometry()
      if (!("hex_id" %in% names(j))) {
        cand <- intersect(c("hex_id.y","hex_id.x","ID","OBJECTID"), names(j))
        if (length(cand)) { j$hex_id <- j[[cand[1]]]; j[cand] <- NULL }
      }
      g <- j |>
        dplyr::filter(.data$MEASYEAR %in% years) |>
        dplyr::group_by(hex_id) |>
        dplyr::summarise(mean_rep = mean(aglb_Mg_per_ha, na.rm = TRUE), .groups = "drop")
      
      g$hex_id <- as.character(g$hex_id)
      
      acc <- dplyr::full_join(acc, g, by = "hex_id") |>
        dplyr::mutate(
          n     = dplyr::coalesce(n, 0L) + 1L,
          sum   = dplyr::coalesce(sum, 0) + dplyr::coalesce(mean_rep, 0),
          sumsq = dplyr::coalesce(sumsq, 0) + dplyr::coalesce(mean_rep^2, 0)
        ) |>
        dplyr::select(hex_id, n, sum, sumsq)
      
      # (optional) write a thinned replicate for QC
      if (thin_every > 0L && (r %% thin_every == 0L) && !is.null(out_dir)) {
        g$replicate_id <- r
        readr::write_csv(g, file.path(out_dir, "replicates", sprintf("rep_%03d.csv", r)))
      }
    }
    
    out <- acc |>
      dplyr::mutate(
        positional_sd = sqrt(pmax(0, sumsq/n - (sum/n)^2)),
        year_label = year_label,
        window = paste0(window_years, "y")
      ) |>
      dplyr::select(hex_id, positional_sd, year_label, window)
    
    return(list(replicates = dplyr::tibble(), summary = out))
  }
  
  # legacy path (persist all replicates in memory)
  res <- vector("list", R)
  for (r in seq_len(R)) {
    u <- runif(nrow(pts_5070)); rr <- sqrt(u) * radius_m; theta <- runif(nrow(pts_5070), 0, 2*pi)
    dx <- rr * cos(theta); dy <- rr * sin(theta)
    pts_j <- sf::st_set_geometry(pts_5070, sf::st_geometry(pts_5070) + cbind(dx, dy))
    
    j <- sf::st_join(pts_j, hx_5070["hex_id"], left = TRUE, join = sf::st_intersects) |> sf::st_drop_geometry()
    if (!("hex_id" %in% names(j))) {
      cand <- intersect(c("hex_id.y","hex_id.x","ID","OBJECTID"), names(j))
      if (length(cand)) { j$hex_id <- j[[cand[1]]]; j[cand] <- NULL }
    }
    g <- j |>
      dplyr::filter(.data$MEASYEAR %in% years) |>
      dplyr::group_by(hex_id) |>
      dplyr::summarise(mean_rep = mean(aglb_Mg_per_ha, na.rm = TRUE), .groups = "drop")
    g$replicate_id <- r
    res[[r]] <- g
    
    if (thin_every > 0L && (r %% thin_every == 0L) && !is.null(out_dir)) {
      readr::write_csv(g, file.path(out_dir, "replicates", sprintf("rep_%03d.csv", r)))
    }
  }
  
  all <- dplyr::bind_rows(res)
  pos <- all |>
    dplyr::group_by(hex_id) |>
    dplyr::summarise(positional_sd = stats::sd(mean_rep, na.rm = TRUE), .groups = "drop") |>
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
  
  # read optional processing config for thinning/format
  cfg <- read_config("configs/process.yml")
  persist_reps <- isTRUE(cfg$persist_replicates %||% FALSE)
  thin_every   <- as.integer(cfg$thin_every %||% 0L)
  
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
  
  # establish run_id/out_dir BEFORE we call positional_mc()
  run_id <- run_id %||% cfg$run_id %||%
    paste0(format(Sys.Date(), "%Y-%m-%d"), "_", metric, "_W", level_window, "y")
  out_dir <- fs::path(paths$runs_dir, run_id)
  if (!fs::dir_exists(out_dir)) fs::dir_create(out_dir, recurse = TRUE)
  message("📁 Output dir: ", out_dir)
  
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
    mc <- positional_mc(
      fia_plot_hex, hex_path, hex_layer,
      year_label = yr, window_years = level_window,
      R = mc_reps, radius_m = jitter_radius_m,
      persist_replicates = persist_reps,
      thin_every = thin_every,
      out_dir = out_dir
    )
    out_pos_sd[[i]] <- mc$summary |> dplyr::mutate(metric = metric)
    out_reps[[i]]   <- mc$replicates |> dplyr::mutate(metric = metric, year_label = yr, window = paste0(level_window,"y"))
    out_design[[i]] <- dplyr::bind_rows(fia_hex_stats, ne_hex_stats, fused_hex_stats) |> dplyr::mutate(metric = metric)
  }
  
  design <- dplyr::bind_rows(out_design) |> dplyr::mutate(hex_id = as.character(hex_id))
  pos_sd <- dplyr::bind_rows(out_pos_sd) |> dplyr::mutate(hex_id = as.character(hex_id))
  reps   <- dplyr::bind_rows(out_reps)   |> dplyr::mutate(hex_id = as.character(hex_id))
  
  design_wide <- design |>
    dplyr::mutate(hex_id = as.character(hex_id)) |>
    dplyr::select(hex_id, year_label, window, source, metric, mean, se, n_plots, sum_weight) |>
    tidyr::pivot_wider(names_from = source, values_from = c(mean, se, n_plots), names_sep = "_")
  
  joined <- design_wide |>
    dplyr::mutate(hex_id = as.character(hex_id)) |>
    dplyr::left_join(pos_sd |> dplyr::mutate(hex_id = as.character(hex_id)) |>
                       dplyr::select(hex_id, year_label, window, positional_sd),
                     by = c("hex_id","year_label","window")) |>
    dplyr::mutate(
      total_sd_fia = sqrt((se_fia %||% NA_real_)^2 + (positional_sd %||% 0)^2),
      se_ratio = dplyr::if_else(is.finite(se_fia) & is.finite(se_fused) & se_fused > 0, se_fia / se_fused, NA_real_)
    )
  
  
  readr::write_csv(joined, fs::path(out_dir, "hex_joined.csv"))
  message("✔ Wrote: ", fs::path(out_dir, "hex_joined.csv"))
  
  if (nrow(reps)) {
    readr::write_csv(reps, fs::path(out_dir, "fia_mc_replicates.csv"))
    message("✔ Wrote: ", fs::path(out_dir, "fia_mc_replicates.csv"))
  } else {
    message("✔ MC replicates not persisted (see 'replicates/' for thinned samples, if enabled).")
  }
  
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
