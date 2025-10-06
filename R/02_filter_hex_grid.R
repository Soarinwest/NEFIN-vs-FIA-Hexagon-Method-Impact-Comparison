#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(sf); library(dplyr); library(stringr); library(fs); library(glue)
})
`%||%` <- function(a,b) if (!is.null(a)) a else b
canon_states <- function(states) {
  states <- toupper(trimws(states))
  
  # 50 states (name/abbr) + add DC once with correct FIPS
  base <- data.frame(
    name = toupper(state.name),
    abbr = toupper(state.abb),
    fips = c(
      "01","02","04","05","06","08","09","10","12","13",
      "15","16","17","18","19","20","21","22","23","24",
      "25","26","27","28","29","30","31","32","33","34",
      "35","36","37","38","39","40","41","42","44","45",
      "46","47","48","49","50","51","53","54","55","56"
    ),
    stringsAsFactors = FALSE
  )
  dc <- data.frame(name = "DISTRICT OF COLUMBIA", abbr = "DC", fips = "11", stringsAsFactors = FALSE)
  map <- unique(rbind(base, dc))
  
  hit <- subset(map, name %in% states | abbr %in% states | fips %in% states)
  if (!nrow(hit)) stop("No valid states parsed from: ", paste(states, collapse=", "))
  list(names = hit$name, abbr = hit$abbr, fips = hit$fips)
}
filter_by_attribute <- function(g_hex, states) {
  st <- canon_states(states)
  cols_up <- toupper(names(g_hex))
  fld_state_name <- names(g_hex)[match("STATE", cols_up)]
  fld_state_fips <- names(g_hex)[match(c("STATE_FIPS","STATEFP","STATEFP10","STATEFP20"), cols_up)]
  fld_state_fips <- fld_state_fips[!is.na(fld_state_fips)][1] %||% NA_character_
  if (!is.na(fld_state_name) && is.character(g_hex[[fld_state_name]])) {
    g_hex |> dplyr::filter(toupper(.data[[fld_state_name]]) %in% st$names)
  } else if (!is.na(fld_state_fips)) {
    val <- g_hex[[fld_state_fips]]
    val <- if (is.numeric(val)) sprintf("%02d", as.integer(val)) else stringr::str_pad(val, 2, pad = "0")
    g_hex |> dplyr::mutate(.fips = val) |> dplyr::filter(.fips %in% st$fips) |> dplyr::select(-.fips)
  } else {
    stop("No STATE or STATE_FIPS-like field found for attribute filtering.")
  }
}
write_geojson_4326 <- function(g, out_path) {
  g4326 <- sf::st_transform(g, 4326)
  fs::dir_create(fs::path_dir(out_path))
  sf::st_write(g4326, out_path, delete_dsn = TRUE, quiet = TRUE)
  message("✔ Wrote: ", out_path, " (", nrow(g4326), " hexes)")
}
filter_hex_grid <- function(in_path, states, out_path = "data/hex/hex_grid.geojson", method = c("attribute"), ensure_hex_id = TRUE) {
  method <- match.arg(method)
  if (!fs::file_exists(in_path)) stop("Input not found: ", in_path)
  message("→ Reading national grid: ", in_path)
  g <- sf::st_read(in_path, quiet = TRUE) |> sf::st_make_valid()
  if (method == "attribute") {
    g_sub <- filter_by_attribute(g, states)
  } else {
    stop("Only 'attribute' method is enabled in this script.")
  }
  if (!nrow(g_sub)) stop("No hex cells matched the requested states.")
  if (ensure_hex_id) {
    if (!("hex_id" %in% names(g_sub))) {
      if ("ID" %in% names(g_sub)) g_sub$hex_id <- as.character(g_sub$ID) else g_sub$hex_id <- sprintf("H%06d", seq_len(nrow(g_sub)))
    }
  }
  write_geojson_4326(g_sub, out_path)
  invisible(out_path)
}
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  get <- function(flag, default=NULL) { hit <- grep(paste0("^",flag,"="), args, value=TRUE); if (length(hit)) sub(paste0("^",flag,"="),"",hit[1]) else default }
  in_path <- get("--in", stop("Missing --in=path/to/national_hex.geojson"))
  states  <- strsplit(get("--states", stop("Missing --states=VT,NH,...")), ",")[[1]]
  out_path <- get("--out", "data/hex/hex_grid.geojson")
  method <- get("--method", "attribute")
  filter_hex_grid(in_path, states, out_path, method = method)
}
