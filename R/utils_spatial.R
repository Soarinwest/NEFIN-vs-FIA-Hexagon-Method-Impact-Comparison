# R/utils_spatial.R
# Shared spatial utilities extracted from process_to_hex.R

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(fs)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

# Find column by candidate names (case-insensitive)
find_col <- function(df, candidates) {
  up <- toupper(names(df))
  for (cand in toupper(candidates)) {
    hit <- which(up == cand)
    if (length(hit)) return(names(df)[hit[1]])
  }
  NA_character_
}

# Normalize plot coordinate columns to lat_public/lon_public
normalize_plot_coords <- function(pl) {
  lat_candidates <- c("lat_public","LAT_PUBLIC","LAT","LATITUDE","DEC_LAT","DEC_LATITUDE","Y")
  lon_candidates <- c("lon_public","LON_PUBLIC","LON","LONGITUDE","DEC_LON","DEC_LONGITUDE","X")
  
  lat_col <- find_col(pl, lat_candidates)
  lon_col <- find_col(pl, lon_candidates)
  
  if (is.na(lat_col) || is.na(lon_col))
    stop("Could not find latitude/longitude columns in PLOT")
  
  pl <- pl |>
    dplyr::mutate(
      lat_public = suppressWarnings(as.numeric(.data[[lat_col]])),
      lon_public = suppressWarnings(as.numeric(.data[[lon_col]]))
    )
  
  pl <- pl |>
    dplyr::mutate(
      lat_public = dplyr::if_else(lat_public >= -90 & lat_public <= 90, lat_public, NA_real_),
      lon_public = dplyr::if_else(lon_public >= -180 & lon_public <= 180, lon_public, NA_real_)
    )
  
  pl
}

# Assign plots to hexes using spatial join
assign_plots_to_hex <- function(df_plots, hex_path, hex_layer = NULL) {
  if (!fs::file_exists(hex_path)) stop("Hex grid not found: ", hex_path)
  
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  hx <- if (is.null(hex_layer)) sf::st_read(hex_path, quiet = TRUE) else sf::st_read(hex_path, layer = hex_layer, quiet = TRUE)
  hx <- sf::st_make_valid(sf::st_zm(hx, drop = TRUE, what = "ZM"))
  if (is.na(sf::st_crs(hx))) sf::st_crs(hx) <- sf::st_crs(4326)
  
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
  
  if ("hex_id" %in% names(pts)) pts$hex_id <- NULL
  
  crs5070 <- sf::st_crs(5070)
  hx_5070  <- sf::st_transform(hx, crs5070)
  pts_5070 <- sf::st_transform(pts, crs5070)
  
  joined <- sf::st_join(pts_5070, hx_5070["hex_id"], left = TRUE, join = sf::st_intersects)
  out <- sf::st_drop_geometry(joined) |> as.data.frame()
  
  if (!("hex_id" %in% names(out))) {
    cand <- intersect(c("hex_id.y","hex_id.x","ID","OBJECTID"), names(out))
    if (length(cand)) {
      out$hex_id <- out[[cand[1]]]
      out[cand] <- NULL
    }
  }
  
  miss <- sum(is.na(out$hex_id))
  message(sprintf("  Assigned %s / %s plots (%.1f%%), %s NA",
                  format(nrow(out) - miss, big.mark = ","),
                  format(nrow(out), big.mark = ","),
                  100 * (nrow(out) - miss) / nrow(out),
                  miss))
  
  out$hex_id <- as.character(out$hex_id)
  return(out)
}

# Build hex union for jitter constraints
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

# Build state polygons for jitter constraints
build_state_polys_5070 <- function(state_geo_path, state_field = "STATEFP") {
  st <- sf::st_read(state_geo_path, quiet = TRUE) |> sf::st_make_valid() |> sf::st_transform(5070)
  if (!state_field %in% names(st)) stop("state_field not found: ", state_field)
  vals <- st[[state_field]]
  if (is.character(vals)) st$STATECD <- suppressWarnings(as.integer(vals)) else st$STATECD <- as.integer(vals)
  st[, c("STATECD","geometry")]
}

# Constrained jitter (keeps points inside hex union and state boundaries)
constrained_jitter_once <- function(pts_5070, radius_m, hex_union_5070,
                                    state_polys_5070 = NULL, max_reroll = 20) {
  n <- nrow(pts_5070)
  if (!n) return(pts_5070)
  
  crs_ref <- sf::st_crs(5070)
  sf::st_crs(pts_5070) <- crs_ref
  sf::st_crs(hex_union_5070) <- crs_ref
  if (!is.null(state_polys_5070)) sf::st_crs(state_polys_5070) <- crs_ref
  
  u <- runif(n); r <- sqrt(u) * radius_m; theta <- runif(n, 0, 2*pi)
  dx <- r * cos(theta); dy <- r * sin(theta)
  moved <- sf::st_set_geometry(pts_5070, sf::st_geometry(pts_5070) + cbind(dx, dy))
  sf::st_crs(moved) <- crs_ref
  
  inside <- suppressWarnings(sf::st_within(moved, hex_union_5070, sparse = FALSE)[,1])
  
  if (!is.null(state_polys_5070) && "STATECD" %in% names(pts_5070)) {
    st_lut <- split(state_polys_5070, state_polys_5070$STATECD)
    chk <- logical(n)
    for (i in seq_len(n)) {
      poly <- st_lut[[as.character(pts_5070$STATECD[i])]]
      if (!is.null(poly) && nrow(poly) > 0) {
        sf::st_crs(poly) <- crs_ref
        chk[i] <- as.logical(sf::st_within(moved[i,], poly, sparse = FALSE)[1,1])
      } else chk[i] <- TRUE
    }
    inside <- inside & chk
  }
  
  tries <- 1L
  while (any(!inside) && tries < max_reroll) {
    idx <- which(!inside)
    u <- runif(length(idx)); r <- sqrt(u) * radius_m; theta <- runif(length(idx), 0, 2*pi)
    dx[idx] <- r * cos(theta); dy[idx] <- r * sin(theta)
    new_geom <- sf::st_geometry(pts_5070[idx, ]) + cbind(dx[idx], dy[idx])
    moved[idx, ] <- sf::st_set_geometry(pts_5070[idx, ], new_geom)
    sf::st_crs(moved) <- crs_ref
    
    inside[idx] <- suppressWarnings(sf::st_within(moved[idx,], hex_union_5070, sparse = FALSE)[,1])
    
    if (!is.null(state_polys_5070) && "STATECD" %in% names(pts_5070)) {
      for (j in seq_along(idx)) {
        i <- idx[j]
        poly <- st_lut[[as.character(pts_5070$STATECD[i])]]
        if (!is.null(poly) && nrow(poly) > 0) {
          sf::st_crs(poly) <- crs_ref
          inside[i] <- inside[i] & as.logical(sf::st_within(moved[i,], poly, sparse = FALSE)[1,1])
        }
      }
    }
    tries <- tries + 1L
  }
  
  if (any(!inside)) {
    bad <- which(!inside)
    for (i in bad) {
      poly <- if (!is.null(state_polys_5070) && "STATECD" %in% names(pts_5070)) {
        p <- st_lut[[as.character(pts_5070$STATECD[i])]]
        if (!is.null(p) && nrow(p) > 0) {
          p_int <- sf::st_intersection(p, hex_union_5070)
          sf::st_crs(p_int) <- crs_ref
          p_int
        } else {
          hex_union_5070
        }
      } else {
        hex_union_5070
      }
      nearest <- sf::st_nearest_points(moved[i,], poly)
      moved[i,] <- sf::st_set_geometry(moved[i,], sf::st_cast(nearest, "POINT")[2])
    }
    sf::st_crs(moved) <- crs_ref
  }
  
  moved
}