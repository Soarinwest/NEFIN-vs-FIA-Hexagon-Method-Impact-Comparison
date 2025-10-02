### run_all.R — master runner.R ###
### Foreword ---------------------------------------------------------------------------------
### Title: run_all.R
### Author: Soren Donisvitch
### Date: 10/02/2025
### Dependents: R (>= 3.5)
### Description: run_all.R — master runner: init - hex filter - FIA pull (state-by-state or single) - hex processing
### Foreword: The use or application of these code without permission of the author is prohibited.The author is not ->
###           liable for the use, modification, or any other application of this or other provided scripts.

suppressPackageStartupMessages({
  library(fs)
  library(yaml)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

# -----------------------------
# 0) Load component scripts
# -----------------------------
safe_source <- function(path) {
  if (!file.exists(path)) stop("Required script missing: ", path, call. = FALSE)
  source(path, local = FALSE)
}

safe_source("R/fiaDataPull.R")        # init_project(), ensure_fiadb_sqlite(), fia_pull(), helpers
# Optional state-by-state pull (only used if configs/fia_states.yml says so)
if (file.exists("R/fiaDataPull_states.R")) safe_source("R/fiaDataPull_states.R")
safe_source("R/filter_hex_grid.R")    # filter_hex_grid()
safe_source("R/process_to_hex.R")     # process_to_hex(); should read from fia_region if present

# -----------------------------
# 1) Initialize project layout
# -----------------------------
paths <- init_project(".")
message("✔ Project initialized at: ", normalizePath("."))

# -----------------------------
# 2) Read configs
# -----------------------------
cfg_hex_path  <- "configs/hex_filter.yml"
cfg_fia_path  <- "configs/fia_pull.yml"     # used for single-bundle mode
cfg_states_path <- "configs/fia_states.yml" # used for state-by-state mode
cfg_proc_path <- "configs/process.yml"

cfg_hex   <- if (file.exists(cfg_hex_path))   yaml::read_yaml(cfg_hex_path)   else list()
cfg_fia   <- if (file.exists(cfg_fia_path))   yaml::read_yaml(cfg_fia_path)   else list()
cfg_states<- if (file.exists(cfg_states_path)) yaml::read_yaml(cfg_states_path) else list()
cfg_proc  <- if (file.exists(cfg_proc_path))  yaml::read_yaml(cfg_proc_path)  else list()

# -----------------------------
# 3) Filter hex grid to states
# -----------------------------
hex_out <- cfg_hex$out_path %||% "data/hex/hex_grid.geojson"
if (!fs::file_exists(hex_out) || isTRUE(cfg_hex$overwrite %||% FALSE)) {
  if (is.null(cfg_hex$national_hex_path)) {
    stop("configs/hex_filter.yml must set 'national_hex_path' to your national hex GeoJSON (EPSG:3857 per your file).")
  }
  dir_create(fs::path_dir(hex_out))
  message("▶ Filtering hex grid to states: ", paste(cfg_hex$states %||% character(), collapse = ", "))
  filter_hex_grid(
    in_path = cfg_hex$national_hex_path,
    states  = cfg_hex$states %||% character(),
    out_path = hex_out,
    method  = cfg_hex$method %||% "attribute"
  )
} else {
  message("✔ Using existing hex grid: ", hex_out)
}

# -----------------------------
# 4) FIA pull
#    Prefer state-by-state if configs/fia_states.yml says use_states: true
# -----------------------------
ran_state_mode <- FALSE
if (length(cfg_states) && isTRUE(cfg_states$use_states)) {
  if (!file.exists("R/fiaDataPull_states.R")) {
    stop("configs/fia_states.yml requests state-by-state pull, but R/fiaDataPull_states.R is missing.")
  }
  message("▶ Running state-by-state FIA pull …")
  # This writes per-state slices under data/interim/fia/states/<ABBR>/ and
  # combines them to data/interim/fia_region/
  fia_pull_states(project_dir = ".", cfg_path = cfg_states_path)
  ran_state_mode <- TRUE
} else {
  message("▶ Running single-bundle FIA pull …")
  DEFAULT_SQLITE_ZIP <- "https://apps.fs.usda.gov/fia/datamart/Databases/SQLite_FIADB_ENTIRE.zip"
  # Years handling: allow [start, end] or an explicit vector
  years_vec <- {
    y <- cfg_fia$years
    if (is.null(y)) 2008:2022 else {
      if (is.list(y) && length(y) == 2) seq(y[[1]], y[[2]]) else unlist(y)
    }
  }
  # States for single-bundle: fall back to hex config states
  states2 <- cfg_fia$states2 %||% (cfg_hex$states %||% c("VT","NH","ME","MA","CT","RI","NY"))
  fia_pull(
    project_dir    = ".",
    states2        = states2,
    years          = years_vec,
    eval_type      = cfg_fia$eval_type %||% "EXPN",
    db_zip_url     = cfg_fia$db_zip_url %||% DEFAULT_SQLITE_ZIP,
    overwrite_zip   = isTRUE(cfg_fia$overwrite_zip %||% FALSE),
    overwrite_unzip = isTRUE(cfg_fia$overwrite_unzip %||% FALSE),
    keep_zip        = isTRUE(cfg_fia$keep_zip %||% FALSE)
  )
}

# -----------------------------
# 5) Process to hex
#    (process_to_hex should auto-prefer data/interim/fia_region if it exists)
# -----------------------------
# Years handling for processing
proc_years <- {
  y <- cfg_proc$years
  if (is.null(y)) 2018:2020 else {
    if (is.list(y) && length(y) == 2) seq(y[[1]], y[[2]]) else unlist(y)
  }
}

message("▶ Processing to hex …")
res <- process_to_hex(
  project_dir     = cfg_proc$project_dir %||% ".",
  hex_path        = hex_out,
  hex_layer       = cfg_proc$hex_layer %||% NULL,
  metric          = cfg_proc$metric %||% "aglb",
  years           = proc_years,
  level_window    = cfg_proc$level_window %||% 3,
  rate_window     = cfg_proc$rate_window %||% 5,
  mc_reps         = cfg_proc$mc_reps %||% 25,
  jitter_radius_m = cfg_proc$jitter_radius_m %||% 1609.34,
  run_id          = cfg_proc$run_id %||% NULL
)

# --- quick visuals for a smoke test ---
if (file.exists("R/make_viz.R")) source("R/make_viz.R")
if (exists("make_run_viz")) {
  make_run_viz(run_dir = res$run_dir, hex_path = hex_out, years = cfg_proc$years %||% NULL)
}

# -----------------------------
# 6) Epilog
# -----------------------------
root <- "data/interim"
fia_root_used <- if (fs::dir_exists(file.path(root, "fia_region"))) file.path(root, "fia_region") else file.path(root, "fia")
cat("✔ Master run complete.\n",
    "  • FIA source: ", fia_root_used, if (ran_state_mode) " (combined state-by-state)" else " (single-bundle)", "\n",
    "  • Hex grid:   ", hex_out, "\n",
    "  • Outputs:    ", dirname(res$hex_joined), "\n", sep = "")

cat("✔ Master run complete. See 'runs/' for outputs.\n")
