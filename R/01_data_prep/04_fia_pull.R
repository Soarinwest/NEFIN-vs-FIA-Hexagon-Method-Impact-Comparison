### Foreword ---------------------------------------------------------------------------------
### Author: Soren Donisvitch
### Date: 10/02/2025
### Dependents: R (>= 3.5), DBI, RSQLite, fs, curl, readr, dplyr, glue, withr, jsonlite, rlang
### Foreword: The use or application of these code without permission of the author is prohibited.The author is not ->
###           liable for the use, modification, or any other application of this or other provided scripts.
### fiaDataPull.R — Download FIADB (zipped SQLite), reuse existing unzipped DB if present,
### unzip then remove the ZIP, and export analysis-ready CSV slices for core tables.
### Can be sourced by a master script (exports init_project() and fia_pull()) or run via CLI.

suppressPackageStartupMessages({
  library(DBI); library(RSQLite); library(fs); library(curl)
  library(readr); library(dplyr); library(glue); library(withr)
  library(jsonlite); library(rlang)
})

# --------------------------
# Defaults & helpers
# --------------------------
DEFAULT_SQLITE_ZIP <- "https://apps.fs.usda.gov/fia/datamart/Databases/SQLite_FIADB_ENTIRE.zip"

`%||%` <- function(a,b) if (!is.null(a)) a else b

abbr_to_statecd <- function(abbr) {
  lut <- c(
    AL=1, AZ=4, AR=5, CA=6, CO=8, CT=9, DE=10, FL=12, GA=13, ID=16, IL=17, IN=18,
    IA=19, KS=20, KY=21, LA=22, ME=23, MD=24, MA=25, MI=26, MN=27, MS=28, MO=29,
    MT=30, NE=31, NV=32, NH=33, NJ=34, NM=35, NY=36, NC=37, ND=38, OH=39, OK=40,
    OR=41, PA=42, RI=44, SC=45, SD=46, TN=47, TX=48, UT=49, VT=50, VA=51, WA=53,
    WV=54, WI=55, WY=56, DC=11, PR=72, VI=78, AK=2, HI=15
  )
  unname(lut[toupper(abbr)])
}

dir_create_safe <- function(path) if (!fs::dir_exists(path)) fs::dir_create(path, recurse = TRUE)

download_zip <- function(url, dest, overwrite = FALSE) {
  if (fs::file_exists(dest) && !overwrite) {
    message("✔ Using existing zip: ", dest)
    return(dest)
  }
  dir_create_safe(fs::path_dir(dest))
  message("↓ Downloading: ", url)
  curl::curl_download(url, destfile = dest, mode = "wb")
  message("✔ Downloaded to: ", dest)
  dest
}

# Prefer existing unzipped DB; unzip zip → then delete zip (by default)
ensure_fiadb_sqlite <- function(db_zip_url,
                                cache_dir,
                                overwrite_unzip = FALSE,
                                overwrite_zip   = FALSE,
                                keep_zip        = FALSE) {
  dir_create_safe(cache_dir)
  unzip_dir <- fs::path(cache_dir, "unzipped")
  dir_create_safe(unzip_dir)

  # 1) Reuse unzipped DB if present (unless overwrite_unzip=TRUE)
  existing <- c(
    fs::dir_ls(unzip_dir, recurse = TRUE, glob="*.db"),
    fs::dir_ls(unzip_dir, recurse = TRUE, glob="*.sqlite")
  )
  if (length(existing) && !overwrite_unzip) {
    db_path <- existing[[1]]
    message("✔ Reusing existing unzipped SQLite: ", db_path)
    return(db_path)
  }

  # 2) If overwriting, clear unzipped dir
  if (overwrite_unzip && fs::dir_exists(unzip_dir)) {
    message("↯ Clearing previous unzipped contents: ", unzip_dir)
    fs::dir_delete(unzip_dir); fs::dir_create(unzip_dir)
  }

  # 3) Ensure ZIP exists, then unzip
  zip_path <- fs::path(cache_dir, basename(db_zip_url))
  download_zip(db_zip_url, zip_path, overwrite = overwrite_zip)

  message("↯ Unzipping to: ", unzip_dir)
  utils::unzip(zipfile = zip_path, exdir = unzip_dir)

  # 4) Locate DB
  cands <- c(
    fs::dir_ls(unzip_dir, recurse = TRUE, glob="*.db"),
    fs::dir_ls(unzip_dir, recurse = TRUE, glob="*.sqlite")
  )
  if (!length(cands)) stop("No .db/.sqlite found after unzip in: ", unzip_dir)
  db_path <- cands[[1]]
  message("✔ SQLite DB ready: ", db_path)

  # 5) Remove the ZIP to save space (unless keep_zip=TRUE)
  if (!keep_zip && fs::file_exists(zip_path)) {
    fs::file_delete(zip_path)
    message("🧹 Removed ZIP to save space: ", zip_path)
  }

  db_path
}

coalesce_latlon <- function(df) {
  lat_cols <- intersect(c("LAT","LATITUDE"), names(df))
  lon_cols <- intersect(c("LON","LONGITUDE"), names(df))
  df |>
    mutate(
      lat_public = dplyr::coalesce(!!!rlang::syms(lat_cols)),
      lon_public = dplyr::coalesce(!!!rlang::syms(lon_cols))
    )
}

# --------------------------
# Project init
# --------------------------
init_project <- function(project_dir = ".") {
  paths <- list(
    project_dir = project_dir,
    data_raw    = fs::path(project_dir, "data", "raw"),
    data_raw_fia = fs::path(project_dir, "data", "raw", "fia_sqlite"),
    data_interim = fs::path(project_dir, "data", "interim"),
    data_interim_fia = fs::path(project_dir, "data", "interim", "fia"),
    data_processed = fs::path(project_dir, "data", "processed"),
    outputs_fig = fs::path(project_dir, "outputs", "figures"),
    outputs_tab = fs::path(project_dir, "outputs", "tables"),
    outputs_maps = fs::path(project_dir, "outputs", "maps"),
    logs = fs::path(project_dir, "logs"),
    configs = fs::path(project_dir, "configs")
  )
  lapply(paths, dir_create_safe)
  # seed a YAML-like text config for reference
  cfg_path <- fs::path(paths$configs, "fia_pull.yml")
  if (!fs::file_exists(cfg_path)) {
    yaml <- paste(
      "states2: [VT, NH, ME, MA, CT, RI, NY]",
      "years: [2008, 2022]",
      "eval_type: EXPN",
      sprintf("db_zip_url: %s", DEFAULT_SQLITE_ZIP),
      "overwrite_zip: false",
      "overwrite_unzip: false",
      "keep_zip: false",
      sep = "
"
    )
    writeLines(yaml, cfg_path)
  }
  paths
}

get_ppsa_plots <- function(con, evalids) {
  if (!nrow(evalids)) return(data.frame())
  tbls <- DBI::dbListTables(con)
  stopifnot("POP_PLOT_STRATUM_ASSGN" %in% tbls)
  cols <- DBI::dbListFields(con, "POP_PLOT_STRATUM_ASSGN")
  pltcol <- if ("PLT_CN" %in% cols) "PLT_CN" else if ("PLOT_CN" %in% cols) "PLOT_CN" else stop("No PLT_CN/PLOT_CN in PPSA")
  sql <- glue::glue("
    SELECT DISTINCT EVALID, {pltcol} AS PLT_CN, STRATUM_CN
    FROM POP_PLOT_STRATUM_ASSGN
    WHERE EVALID IN ({paste(unique(evalids$EVALID), collapse=',')})
  ")
  DBI::dbGetQuery(con, sql)
}

# --------------------------
# FIADB routines
# --------------------------
# get_evalids(): robust across state DB variants
get_evalids <- function(con, states2, years, eval_type = "EXPN") {
  statecds <- abbr_to_statecd(states2); statecds <- statecds[!is.na(statecds)]
  if (!length(statecds)) stop("No valid STATECDs parsed from: ", paste(states2, collapse=","))
  
  yrs <- sort(unique(as.integer(years)))
  tbls <- DBI::dbListTables(con)
  if (!"POP_EVAL" %in% tbls) stop("POP_EVAL not found")
  
  pe_cols <- DBI::dbListFields(con, "POP_EVAL")
  state_col <- if ("STATECD" %in% pe_cols) "STATECD" else stop("STATECD not found in POP_EVAL")
  
  # year columns present in your screenshot
  have_report <- "REPORT_YEAR_NM" %in% pe_cols
  have_end    <- "END_INVYR" %in% pe_cols
  
  # year filter
  yr_filter <- if (have_report) {
    glue::glue("AND CAST(e.REPORT_YEAR_NM AS INTEGER) IN ({paste(yrs, collapse=',')})")
  } else if (have_end) {
    glue::glue("AND e.END_INVYR IN ({paste(yrs, collapse=',')})")
  } else {
    "" # last resort: no year filter here
  }
  
  sel_rpt  <- if (have_report) "CAST(e.REPORT_YEAR_NM AS INTEGER) AS RPT_YR," else "NULL AS RPT_YR,"
  sel_estn <- if ("ESTN_YR" %in% pe_cols) "e.ESTN_YR AS ESTN_YR," else if (have_end) "e.END_INVYR AS ESTN_YR," else "NULL AS ESTN_YR,"
  
  # base (no type filter)
  base_sql <- glue::glue("
    SELECT DISTINCT e.EVALID, e.EVAL_GRP_CN,
           {sel_rpt}
           {sel_estn}
           e.{state_col} AS STATECD
    FROM POP_EVAL e
    WHERE e.{state_col} IN ({paste(statecds, collapse=',')})
      {yr_filter}
  ")
  
  # try to apply eval_type if possible
  if (!is.null(eval_type) && nzchar(eval_type) && "POP_EVAL_TYP" %in% tbls) {
    pet_cols <- DBI::dbListFields(con, "POP_EVAL_TYP")
    typ_col  <- intersect(c("EVAL_TYP","EVALTYPE","EVAL_TYP_CD","EVAL_TYPE","EVALTYP","EVAL_TYP_NM"), pet_cols)
    typ_col  <- if (length(typ_col)) typ_col[1] else NA_character_
    
    # usable join?
    join_on_evalid <- ("EVALID" %in% pet_cols) && ("EVALID" %in% pe_cols)
    join_on_grp    <- ("EVAL_GRP_CN" %in% pet_cols) && ("EVAL_GRP_CN" %in% pe_cols)
    
    if (!is.na(typ_col) && (join_on_evalid || join_on_grp)) {
      sql <- if (join_on_evalid) {
        glue::glue("
          SELECT DISTINCT e.EVALID, e.EVAL_GRP_CN,
                 {sel_rpt}
                 {sel_estn}
                 e.{state_col} AS STATECD
          FROM POP_EVAL e
          JOIN POP_EVAL_TYP t ON e.EVALID = t.EVALID
          WHERE e.{state_col} IN ({paste(statecds, collapse=',')})
            {yr_filter}
            AND t.{typ_col} = {DBI::dbQuoteString(con, eval_type)}
        ")
      } else {
        glue::glue("
          SELECT DISTINCT e.EVALID, e.EVAL_GRP_CN,
                 {sel_rpt}
                 {sel_estn}
                 e.{state_col} AS STATECD
          FROM POP_EVAL e
          JOIN POP_EVAL_TYP t ON e.EVAL_GRP_CN = t.EVAL_GRP_CN
          WHERE e.{state_col} IN ({paste(statecds, collapse=',')})
            {yr_filter}
            AND t.{typ_col} = {DBI::dbQuoteString(con, eval_type)}
        ")
      }
      out <- DBI::dbGetQuery(con, sql)
      if (nrow(out)) return(out)
      message("   • No rows after eval_type join; falling back to no-type filter.")
    } else {
      message("   • POP_EVAL_TYP present but no usable join/type col; skipping type filter.")
    }
  } else if (!is.null(eval_type) && nzchar(eval_type)) {
    message("   • POP_EVAL_TYP not present; skipping type filter.")
  }
  
  DBI::dbGetQuery(con, base_sql)
}

export_core_slices <- function(con, evalids_df, out_dir) {
  dir_create_safe(out_dir)
  readr::write_csv(evalids_df, fs::path(out_dir, "pop_eval_filtered.csv"))
  readr::write_csv(DBI::dbGetQuery(con, "SELECT * FROM POP_STRATUM"),
                   fs::path(out_dir, "pop_stratum.csv"))
  readr::write_csv(DBI::dbGetQuery(con, "SELECT * FROM POP_ESTN_UNIT"),
                   fs::path(out_dir, "pop_estn_unit.csv"))

  evalid_list <- paste(unique(evalids_df$EVALID), collapse = ",")
  pps <- DBI::dbGetQuery(con, glue("
    SELECT * FROM POP_PLOT_STRATUM_ASSGN
    WHERE EVALID IN ({evalid_list})
  "))
  readr::write_csv(pps, fs::path(out_dir, "pop_plot_stratum_assgn.csv"))

  # temp table of PLT_CN filtered by our EVALIDs
  DBI::dbExecute(con, "DROP TABLE IF EXISTS __PLTS__")
  DBI::dbExecute(con, "CREATE TEMP TABLE __PLTS__ (PLT_CN INTEGER PRIMARY KEY)")
  DBI::dbExecute(con, glue("
    INSERT INTO __PLTS__
    SELECT DISTINCT PLT_CN FROM POP_PLOT_STRATUM_ASSGN WHERE EVALID IN ({evalid_list})
  "))

  plot_df <- DBI::dbGetQuery(con, "SELECT * FROM PLOT WHERE PLT_CN IN (SELECT PLT_CN FROM __PLTS__)")
  plot_df <- coalesce_latlon(plot_df)
  readr::write_csv(plot_df, fs::path(out_dir, "plot.csv"))

  cond_df <- DBI::dbGetQuery(con, "SELECT * FROM COND WHERE PLT_CN IN (SELECT PLT_CN FROM __PLTS__)")
  readr::write_csv(cond_df, fs::path(out_dir, "cond.csv"))

  tree_df <- DBI::dbGetQuery(con, "SELECT * FROM TREE WHERE PLT_CN IN (SELECT PLT_CN FROM __PLTS__)")
  readr::write_csv(tree_df, fs::path(out_dir, "tree.csv"))

  if ("SEEDLING" %in% DBI::dbListTables(con)) {
    seed_df <- DBI::dbGetQuery(con, "SELECT * FROM SEEDLING WHERE PLT_CN IN (SELECT PLT_CN FROM __PLTS__)")
    readr::write_csv(seed_df, fs::path(out_dir, "seedling.csv"))
  }
  message("✔ Core FIADB slices written to: ", out_dir)
  invisible(TRUE)
}

# --------------------------
# Orchestrator
# --------------------------
fia_pull <- function(project_dir = ".",
                     states2 = c("VT","NH","ME","MA","CT","RI","NY"),
                     years   = 2008:2022,
                     eval_type = "EXPN",
                     db_zip_url = DEFAULT_SQLITE_ZIP,
                     overwrite_zip = FALSE,
                     overwrite_unzip = FALSE,
                     keep_zip = FALSE) {

  paths <- init_project(project_dir)
  log_file <- fs::path(paths$logs, sprintf("fia_pull_%s.log", format(Sys.time(), "%Y%m%d-%H%M%S")))
  sink(log_file, split = TRUE); on.exit(sink(), add = TRUE)

  message("▶ FIADB pull starting …")
  message("  project_dir: ", normalizePath(project_dir))

  # 1) Ensure/reuse SQLite
  db_path <- ensure_fiadb_sqlite(
    db_zip_url      = db_zip_url,
    cache_dir       = paths$data_raw_fia,
    overwrite_unzip = overwrite_unzip,
    overwrite_zip   = overwrite_zip,
    keep_zip        = keep_zip
  )

  # 2) Connect + checks
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path, loadable.extensions = FALSE)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  needed <- c("PLOT","COND","TREE","POP_EVAL","POP_EVAL_TYP",
              "POP_PLOT_STRATUM_ASSGN","POP_STRATUM","POP_ESTN_UNIT")
  missing <- setdiff(needed, DBI::dbListTables(con))
  if (length(missing)) stop("Missing FIADB tables: ", paste(missing, collapse=", "))

  # 3) Filter EVALIDs
  evals <- get_evalids(con, states2, years, eval_type)
  if (nrow(evals) == 0) stop("No POP_EVAL rows for states=", paste(states2, collapse=","),
                             " years=", paste(range(years), collapse=":"))

  # 4) Export slices
  export_core_slices(con, evals, out_dir = paths$data_interim_fia)

  # 5) Manifest
  manifest <- list(
    states = states2, years = years, eval_type = eval_type,
    db_zip_url = db_zip_url, db_file = db_path,
    outputs = list(
      dir = paths$data_interim_fia,
      files = c("plot.csv","cond.csv","tree.csv","seedling.csv",
                "pop_eval_filtered.csv","pop_plot_stratum_assgn.csv",
                "pop_stratum.csv","pop_estn_unit.csv")
    ),
    created_at = as.character(Sys.time())
  )
  writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE),
             fs::path(paths$data_interim_fia, "manifest.json"))

  message("✔ Done. Manifest: ", fs::path(paths$data_interim_fia, "manifest.json"))
  invisible(list(paths = paths, manifest = manifest))
}

# --------------------------
# CLI entrypoint
# --------------------------
parse_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  get <- function(flag, default=NULL) { hit <- grep(paste0("^",flag,"="), args, value=TRUE); if (length(hit)) sub(paste0("^",flag,"="),"",hit[1]) else default }
  states2 <- strsplit(get("--states","VT,NH,ME,MA,CT,RI,NY"), ",")[[1]] |> toupper() |> trimws()
  years_s <- get("--years","2008:2022"); years_rng <- as.integer(strsplit(years_s,":")[[1]]); years <- seq(min(years_rng), max(years_rng))
  list(
    project_dir = get("--project","."),
    states2 = states2,
    years = years,
    eval_type = get("--eval","EXPN"),
    db_zip_url = get("--url", DEFAULT_SQLITE_ZIP),
    overwrite_zip   = isTRUE(as.logical(get("--overwrite_zip","FALSE"))),
    overwrite_unzip = isTRUE(as.logical(get("--overwrite_unzip","FALSE"))),
    keep_zip        = isTRUE(as.logical(get("--keep_zip","FALSE")))
  )
}

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  cli <- parse_cli()
  fia_pull(
    project_dir = cli$project_dir,
    states2 = cli$states2,
    years = cli$years,
    eval_type = cli$eval_type,
    db_zip_url = cli$db_zip_url,
    overwrite_zip = cli$overwrite_zip,
    overwrite_unzip = cli$overwrite_unzip,
    keep_zip = cli$keep_zip
  )
}
