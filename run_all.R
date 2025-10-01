#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(fs); library(yaml) })
`%||%` <- function(a,b) if (!is.null(a)) a else b
source("R/fiaDataPull.R")
source("R/process_to_hex.R")
source("R/filter_hex_grid.R")

paths <- init_project(".")
cfg_hex <- yaml::read_yaml("configs/hex_filter.yml")
cfg_fia <- if (fs::file_exists("configs/fia_pull.yml")) yaml::read_yaml("configs/fia_pull.yml") else list()
cfg_proc <- if (fs::file_exists("configs/process.yml")) yaml::read_yaml("configs/process.yml") else list()

hex_out <- cfg_hex$out_path %||% "data/hex/hex_grid.geojson"
if (!fs::file_exists(hex_out) || isTRUE(cfg_hex$overwrite)) {
  fs::dir_create(fs::path_dir(hex_out))
  message("▶ Filtering hex grid to states: ", paste(cfg_hex$states, collapse=", "))
  filter_hex_grid(cfg_hex$national_hex_path, cfg_hex$states, hex_out, method = cfg_hex$method %||% "attribute")
} else {
  message("✔ Using existing hex grid: ", hex_out)
}

DEFAULT_SQLITE_ZIP <- "https://apps.fs.usda.gov/fia/datamart/Databases/SQLite_FIADB_ENTIRE.zip"
fia_pull(
  project_dir = ".",
  states2 = cfg_fia$states2 %||% cfg_hex$states %||% c("VT","NH","ME","MA","CT","RI","NY"),
  years   = if (!is.null(cfg_fia$years) && length(cfg_fia$years) == 2) seq(cfg_fia$years[[1]], cfg_fia$years[[2]]) else 2008:2022,
  eval_type = cfg_fia$eval_type %||% "EXPN",
  db_zip_url = cfg_fia$db_zip_url %||% DEFAULT_SQLITE_ZIP,
  overwrite_zip = isTRUE(cfg_fia$overwrite_zip %||% FALSE),
  overwrite_unzip = isTRUE(cfg_fia$overwrite_unzip %||% FALSE),
  keep_zip = isTRUE(cfg_fia$keep_zip %||% FALSE)
)

process_to_hex(
  project_dir = cfg_proc$project_dir %||% ".",
  hex_path = hex_out,
  hex_layer = cfg_proc$hex_layer %||% NULL,
  metric = cfg_proc$metric %||% "aglb",
  years = cfg_proc$years %||% 2018:2020,
  level_window = cfg_proc$level_window %||% 3,
  rate_window = cfg_proc$rate_window %||% 5,
  mc_reps = cfg_proc$mc_reps %||% 25,
  jitter_radius_m = cfg_proc$jitter_radius_m %||% 1609.34,
  run_id = cfg_proc$run_id %||% NULL
)
cat("✔ Master run complete. See 'runs/' for outputs.\n")
