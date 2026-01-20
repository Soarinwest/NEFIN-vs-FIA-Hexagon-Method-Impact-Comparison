# fiaDataPull_states.R — state-by-state FIADB export with column-aware joins and SQL debug logs
### Foreword ---------------------------------------------------------------------------------
### Author: Soren Donisvitch
### Date: 10/02/2025
### Dependents: R (>= 3.5), DBI, RSQLite, fs, curl, readr, dplyr, glue, withr, jsonlite, rlang
### Foreword: The use or application of these code without permission of the author is prohibited.The author is not ->
###           liable for the use, modification, or any other application of this or other provided scripts.

suppressPackageStartupMessages({
  library(fs); library(yaml); library(dplyr); library(readr); library(glue); library(DBI); library(RSQLite)
})

source("R/03_fia_pull.R")

`%||%` <- function(a,b) if (!is.null(a)) a else b
dir_create_safe <- function(p) if (!fs::dir_exists(p)) fs::dir_create(p, recurse=TRUE)

# ----------------------------- DEBUG LOG HELPERS --------------------------------
sql_log_dir <- function(state) {
  d <- fs::path("logs","sql", state); dir_create_safe(d); d
}
log_sql <- function(state, name, sql) {
  fp <- fs::path(sql_log_dir(state), paste0(name, ".sql"))
  writeLines(sql, fp); invisible(fp)
}
log_head <- function(state, name, df, n = 50) {
  fp <- fs::path(sql_log_dir(state), paste0(name, "_head.csv"))
  suppressMessages(readr::write_csv(head(df, n), fp)); invisible(fp)
}
log_cols <- function(state, table, cols) {
  fp <- fs::path(sql_log_dir(state), paste0("cols_", table, ".txt"))
  writeLines(paste(cols, collapse = "\n"), fp); invisible(fp)
}

# ----------------------------- ROBUST EVALID DISCOVERY --------------------------
get_evalids <- function(con, states2, years, eval_type = "EXPN") {
  statecds <- abbr_to_statecd(states2); statecds <- statecds[!is.na(statecds)]
  if (!length(statecds)) stop("No valid STATECDs parsed from: ", paste(states2, collapse=","))
  
  yrs <- sort(unique(as.integer(years)))
  tbls <- DBI::dbListTables(con)
  if (!"POP_EVAL" %in% tbls) stop("POP_EVAL not found")
  
  pe_cols <- DBI::dbListFields(con, "POP_EVAL")
  state_col <- if ("STATECD" %in% pe_cols) "STATECD" else stop("STATECD not found in POP_EVAL")
  
  have_report <- "REPORT_YEAR_NM" %in% pe_cols
  have_end    <- "END_INVYR" %in% pe_cols
  yr_filter <- if (have_report) {
    glue("AND CAST(e.REPORT_YEAR_NM AS INTEGER) IN ({paste(yrs, collapse=',')})")
  } else if (have_end) {
    glue("AND e.END_INVYR IN ({paste(yrs, collapse=',')})")
  } else {
    ""
  }
  sel_rpt  <- if (have_report) "CAST(e.REPORT_YEAR_NM AS INTEGER) AS RPT_YR," else "NULL AS RPT_YR,"
  sel_estn <- if ("ESTN_YR" %in% pe_cols) "e.ESTN_YR AS ESTN_YR," else if (have_end) "e.END_INVYR AS ESTN_YR," else "NULL AS ESTN_YR,"
  
  base_sql <- glue("
    SELECT DISTINCT e.EVALID, e.EVAL_GRP_CN,
           {sel_rpt}
           {sel_estn}
           e.{state_col} AS STATECD
    FROM POP_EVAL e
    WHERE e.{state_col} IN ({paste(statecds, collapse=',')})
      {yr_filter}
  ")
  
  # optional eval_type filter via POP_EVAL_TYP
  if (!is.null(eval_type) && nzchar(eval_type) && "POP_EVAL_TYP" %in% tbls) {
    pet_cols <- DBI::dbListFields(con, "POP_EVAL_TYP")
    typ_col  <- intersect(c("EVAL_TYP","EVALTYPE","EVAL_TYP_CD","EVAL_TYPE","EVALTYP","EVAL_TYP_NM"), pet_cols)
    typ_col  <- if (length(typ_col)) typ_col[1] else NA_character_
    
    join_on_evalid <- ("EVALID" %in% pet_cols) && ("EVALID" %in% pe_cols)
    join_on_grp    <- ("EVAL_GRP_CN" %in% pet_cols) && ("EVAL_GRP_CN" %in% pe_cols)
    
    if (!is.na(typ_col) && (join_on_evalid || join_on_grp)) {
      sql <- if (join_on_evalid) {
        glue("
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
        glue("
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

# ----------------------------- PPSA (plots per EVALID) -------------------------
get_ppsa_plots <- function(con, evalids, state_abbr, debug = TRUE) {
  if (!nrow(evalids)) return(tibble())
  tbls <- DBI::dbListTables(con)
  if (!"POP_PLOT_STRATUM_ASSGN" %in% tbls) stop("POP_PLOT_STRATUM_ASSGN not found")
  
  cols <- DBI::dbListFields(con, "POP_PLOT_STRATUM_ASSGN")
  if (debug) log_cols(state_abbr, "POP_PLOT_STRATUM_ASSGN", cols)
  
  pltcol <- if ("PLT_CN" %in% cols) "PLT_CN" else if ("PLOT_CN" %in% cols) "PLOT_CN" else stop("No PLT_CN/PLOT_CN in PPSA")
  sql <- glue("
    WITH ppsa AS (
      SELECT DISTINCT EVALID, {pltcol} AS PLT_CN, STRATUM_CN
      FROM POP_PLOT_STRATUM_ASSGN
      WHERE EVALID IN ({paste(unique(evalids$EVALID), collapse=',')})
    )
    SELECT * FROM ppsa
  ")
  if (debug) log_sql(state_abbr, "01_ppsa", sql)
  out <- DBI::dbGetQuery(con, sql)
  if (debug) log_head(state_abbr, "01_ppsa", out)
  out
}

# ----------------------------- EXPORT SLICES (STATE) ---------------------------
export_core_slices_state <- function(con, evalids, out_dir, state_abbr, debug = TRUE) {
  dir_create_safe(out_dir)
  
  # 1) ppsa (EVALID, PLT_CN, STRATUM_CN)
  ppsa <- get_ppsa_plots(con, evalids, state_abbr, debug = debug)
  if (!nrow(ppsa)) { warning("No PPSA rows for ", state_abbr); return(invisible(NULL)) }
  
  # 2) PLOT subset via join to ppsa (PLOT.CN matches PLT_CN)
  sql_plot <- "
    WITH ppsa AS (
      SELECT DISTINCT EVALID, PLT_CN, STRATUM_CN FROM POP_PLOT_STRATUM_ASSGN
    )
    SELECT p.*
    FROM PLOT p
    JOIN (SELECT DISTINCT PLT_CN FROM ppsa) x ON p.CN = x.PLT_CN
  "
  # but the CTE above assumes we already aliased PLT_CN; we’ll rebuild with the right source column:
  tbl_cols <- DBI::dbListFields(con, "POP_PLOT_STRATUM_ASSGN")
  src_plt <- if ("PLT_CN" %in% tbl_cols) "PLT_CN" else "PLOT_CN"
  sql_plot <- glue("
    WITH ppsa AS (
      SELECT DISTINCT EVALID, {src_plt} AS PLT_CN, STRATUM_CN
      FROM POP_PLOT_STRATUM_ASSGN
      WHERE EVALID IN ({paste(unique(evalids$EVALID), collapse=',')})
    )
    SELECT p.*
    FROM PLOT p
    JOIN (SELECT DISTINCT PLT_CN FROM ppsa) x ON p.CN = x.PLT_CN
  ")
  if (debug) log_sql(state_abbr, "02_plot", sql_plot)
  plot_df <- DBI::dbGetQuery(con, sql_plot)
  if (debug) log_head(state_abbr, "02_plot", plot_df)
  if (!nrow(plot_df)) warning("PLOT subset came back empty for ", state_abbr)
  
  # 3) TREE subset via join to ppsa (TREE.PLT_CN matches ppsa.PLT_CN)
  sql_tree <- glue("
    WITH ppsa AS (
      SELECT DISTINCT EVALID, {src_plt} AS PLT_CN, STRATUM_CN
      FROM POP_PLOT_STRATUM_ASSGN
      WHERE EVALID IN ({paste(unique(evalids$EVALID), collapse=',')})
    )
    SELECT t.*
    FROM TREE t
    JOIN (SELECT DISTINCT PLT_CN FROM ppsa) x ON t.PLT_CN = x.PLT_CN
  ")
  if (debug) log_sql(state_abbr, "03_tree", sql_tree)
  tree_df <- DBI::dbGetQuery(con, sql_tree)
  if (debug) log_head(state_abbr, "03_tree", tree_df)
  
  # 4) COND subset (optional but useful later)
  have_cond <- "COND" %in% DBI::dbListTables(con)
  if (have_cond) {
    sql_cond <- glue("
      WITH ppsa AS (
        SELECT DISTINCT EVALID, {src_plt} AS PLT_CN, STRATUM_CN
        FROM POP_PLOT_STRATUM_ASSGN
        WHERE EVALID IN ({paste(unique(evalids$EVALID), collapse=',')})
      )
      SELECT c.*
      FROM COND c
      JOIN (SELECT DISTINCT PLT_CN FROM ppsa) x ON c.PLT_CN = x.PLT_CN
    ")
    if (debug) log_sql(state_abbr, "04_cond", sql_cond)
    cond_df <- DBI::dbGetQuery(con, sql_cond)
    if (debug) log_head(state_abbr, "04_cond", cond_df)
  } else {
    cond_df <- tibble()
  }
  
  # 5) SEEDLING subset (if exists)
  have_seed <- "SEEDLING" %in% DBI::dbListTables(con)
  if (have_seed) {
    sql_seed <- glue("
      WITH ppsa AS (
        SELECT DISTINCT EVALID, {src_plt} AS PLT_CN, STRATUM_CN
        FROM POP_PLOT_STRATUM_ASSGN
        WHERE EVALID IN ({paste(unique(evalids$EVALID), collapse=',')})
      )
      SELECT s.*
      FROM SEEDLING s
      JOIN (SELECT DISTINCT PLT_CN FROM ppsa) x ON s.PLT_CN = x.PLT_CN
    ")
    if (debug) log_sql(state_abbr, "05_seedling", sql_seed)
    seed_df <- DBI::dbGetQuery(con, sql_seed)
    if (debug) log_head(state_abbr, "05_seedling", seed_df)
  } else {
    seed_df <- tibble()
  }
  
  # 6) POP tables (filtered)
  # Column introspection for POP_STRATUM
  ps_cols <- DBI::dbListFields(con, "POP_STRATUM")
  if (debug) log_cols(state_abbr, "POP_STRATUM", ps_cols)
  
  # decide the key column in POP_STRATUM (varies by state drops)
  ps_key <- if ("STRATUM_CN" %in% ps_cols) {
    "STRATUM_CN"
  } else if ("CN" %in% ps_cols) {
    "CN"             # many DBs use CN as the key here
  } else {
    stop("POP_STRATUM has neither STRATUM_CN nor CN — can’t join.")
  }
  
  # Build POP_STRATUM SQL joining PPSA’s STRATUM_CN to the detected key
  sql_stratum <- glue("
  WITH ppsa AS (
    SELECT DISTINCT EVALID, {src_plt} AS PLT_CN, STRATUM_CN
    FROM POP_PLOT_STRATUM_ASSGN
    WHERE EVALID IN ({paste(unique(evalids$EVALID), collapse=',')})
  )
  SELECT s.*
  FROM POP_STRATUM s
  JOIN (SELECT DISTINCT STRATUM_CN FROM ppsa) x ON s.{ps_key} = x.STRATUM_CN
")
  if (debug) log_sql(state_abbr, "06_pop_stratum", sql_stratum)
  pop_stratum <- DBI::dbGetQuery(con, sql_stratum)
  if (debug) log_head(state_abbr, "06_pop_stratum", pop_stratum)
  
  # POP_ESTN_UNIT: join via the correct key
  eu_cols <- DBI::dbListFields(con, "POP_ESTN_UNIT")
  if (debug) log_cols(state_abbr, "POP_ESTN_UNIT", eu_cols)
  
  eu_key <- if ("CN" %in% eu_cols) {
    "CN"               # common: EU table keyed by CN
  } else if ("ESTN_UNIT_CN" %in% eu_cols) {
    "ESTN_UNIT_CN"
  } else {
    stop("POP_ESTN_UNIT has neither CN nor ESTN_UNIT_CN — can’t join.")
  }
  
  if (nrow(pop_stratum) && "ESTN_UNIT_CN" %in% names(pop_stratum)) {
    eus <- unique(pop_stratum$ESTN_UNIT_CN)
    eus <- eus[is.finite(eus)]
    if (length(eus)) {
      sql_eu <- glue("SELECT * FROM POP_ESTN_UNIT WHERE {eu_key} IN ({paste(eus, collapse=',')})")
      if (debug) log_sql(state_abbr, "07_pop_estn_unit", sql_eu)
      pop_eu <- DBI::dbGetQuery(con, sql_eu)
      if (debug) log_head(state_abbr, "07_pop_estn_unit", pop_eu)
    } else {
      pop_eu <- tibble()
    }
  } else {
    pop_eu <- tibble()
  }
  
  # POP_EVAL (filtered by EVALIDs we used)
  evalids_in <- paste(unique(evalids$EVALID), collapse=",")
  sql_eval <- glue("SELECT * FROM POP_EVAL WHERE EVALID IN ({evalids_in})")
  if (debug) log_sql(state_abbr, "08_pop_eval", sql_eval)
  pop_eval <- DBI::dbGetQuery(con, sql_eval)
  if (debug) log_head(state_abbr, "08_pop_eval", pop_eval)
  
  # POP_EVAL_TYP if present (rows for these EVALIDs or eval groups)
  if ("POP_EVAL_TYP" %in% DBI::dbListTables(con)) {
    pet_cols <- DBI::dbListFields(con, "POP_EVAL_TYP")
    if ("EVALID" %in% pet_cols) {
      sql_pet <- glue("SELECT * FROM POP_EVAL_TYP WHERE EVALID IN ({evalids_in})")
    } else if ("EVAL_GRP_CN" %in% pet_cols && "EVAL_GRP_CN" %in% names(pop_eval)) {
      grp_in <- paste(unique(pop_eval$EVAL_GRP_CN), collapse=",")
      sql_pet <- glue("SELECT * FROM POP_EVAL_TYP WHERE EVAL_GRP_CN IN ({grp_in})")
    } else {
      sql_pet <- NULL
    }
    if (!is.null(sql_pet)) {
      if (debug) log_sql(state_abbr, "09_pop_eval_typ", sql_pet)
      pop_eval_typ <- DBI::dbGetQuery(con, sql_pet)
      if (debug) log_head(state_abbr, "09_pop_eval_typ", pop_eval_typ)
    } else {
      pop_eval_typ <- tibble()
    }
  } else {
    pop_eval_typ <- tibble()
  }
  
  
  # 7) Write CSVs
  if (nrow(plot_df))     readr::write_csv(plot_df, fs::path(out_dir, "plot.csv"))
  if (nrow(tree_df))     readr::write_csv(tree_df, fs::path(out_dir, "tree.csv"))
  if (nrow(cond_df))     readr::write_csv(cond_df, fs::path(out_dir, "cond.csv"))
  if (nrow(seed_df))     readr::write_csv(seed_df, fs::path(out_dir, "seedling.csv"))
  if (nrow(pop_eval))    readr::write_csv(pop_eval, fs::path(out_dir, "pop_eval_filtered.csv"))
  if (nrow(pop_eval_typ))readr::write_csv(pop_eval_typ, fs::path(out_dir, "pop_eval_typ.csv"))
  if (nrow(pop_stratum)) readr::write_csv(pop_stratum, fs::path(out_dir, "pop_stratum.csv"))
  if (nrow(pop_eu))      readr::write_csv(pop_eu, fs::path(out_dir, "pop_estn_unit.csv"))
  
  invisible(TRUE)
}

# ----------------------------- PER-STATE DRIVER --------------------------------
pull_one_state <- function(project_dir, abbr, url, years, eval_type, overwrite_zip, overwrite_unzip, keep_zip) {
  paths <- list(
    raw_state   = fs::path(project_dir, "data", "raw", "fia_sqlite", abbr),
    out_state   = fs::path(project_dir, "data", "interim", "fia", "states", abbr)
  )
  lapply(paths, dir_create_safe)
  
  db_path <- ensure_fiadb_sqlite(
    db_zip_url      = url,
    cache_dir       = paths$raw_state,
    overwrite_unzip = overwrite_unzip,
    overwrite_zip   = overwrite_zip,
    keep_zip        = keep_zip
  )
  
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path, loadable.extensions = FALSE)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  evals <- get_evalids(con, states2 = abbr, years = years, eval_type = eval_type)
  if (nrow(evals) == 0) {
    message("⚠ No evaluations matched for ", abbr, " in years ", paste(range(years), collapse=":"))
    return(invisible(NULL))
  }
  
  # NEW: export state slices with robust joins + SQL logs
  export_core_slices_state(con, evals, out_dir = paths$out_state, state_abbr = abbr, debug = TRUE)
  
  # annotate source state (for downstream troubleshooting)
  for (f in c("plot.csv","cond.csv","tree.csv","seedling.csv",
              "pop_eval_filtered.csv","pop_eval_typ.csv","pop_stratum.csv","pop_estn_unit.csv")) {
    fp <- fs::path(paths$out_state, f)
    if (fs::file_exists(fp)) {
      df <- readr::read_csv(fp, show_col_types = FALSE)
      if (!("SRC_STATE" %in% names(df))) df$SRC_STATE <- abbr
      readr::write_csv(df, fp)
    }
  }
  
  invisible(paths$out_state)
}

# ----------------------------- COMBINE STATES ----------------------------------
combine_states_to_region <- function(project_dir) {
  root_states <- fs::path(project_dir, "data", "interim", "fia", "states")
  out_region  <- fs::path(project_dir, "data", "interim", "fia_region")
  if (!fs::dir_exists(root_states)) stop("No per-state slices found at: ", root_states)
  if (!fs::dir_exists(out_region)) fs::dir_create(out_region, recurse = TRUE)
  
  # --- normalizers --------------------------------------------------------
  norm_identity <- function(df) df
  
  norm_pop_stratum <- function(df) {
    # Some DBs use CN as the stratum key; normalize to STRATUM_CN
    if (!("STRATUM_CN" %in% names(df)) && "CN" %in% names(df)) {
      df <- dplyr::rename(df, STRATUM_CN = CN)
    }
    # Some DBs may name EU key as CN; normalize to ESTN_UNIT_CN
    if (!("ESTN_UNIT_CN" %in% names(df)) && "ESTN_UNIT" %in% names(df)) {
      df <- dplyr::rename(df, ESTN_UNIT_CN = ESTN_UNIT)
    }
    df
  }
  
  norm_pop_estn_unit <- function(df) {
    # Normalize EU key to ESTN_UNIT_CN
    if (!("ESTN_UNIT_CN" %in% names(df)) && "CN" %in% names(df)) {
      df <- dplyr::rename(df, ESTN_UNIT_CN = CN)
    }
    df
  }
  
  # robust binder with optional normalizer and optional distinct on available keys
  bind_many <- function(fname, .distinct = FALSE, .keys = NULL, .norm = norm_identity) {
    files <- tryCatch(fs::dir_ls(root_states, recurse = TRUE, type = "file", glob = paste0("*/", fname)),
                      error = function(e) character(0))
    if (!length(files)) return(invisible(NULL))
    
    dfs <- lapply(files, function(p) suppressMessages(readr::read_csv(p, show_col_types = FALSE)))
    dfs <- lapply(dfs, .norm)
    all <- dplyr::bind_rows(dfs)
    
    if (.distinct && !is.null(.keys)) {
      keys_in <- intersect(.keys, names(all))
      if (length(keys_in)) {
        all <- dplyr::distinct(all, dplyr::across(dplyr::all_of(keys_in)), .keep_all = TRUE)
      } else {
        message("▶ combine_states_to_region: skipping distinct for ", fname, " — none of keys present: ", paste(.keys, collapse=", "))
      }
    }
    out <- fs::path(out_region, fname)
    readr::write_csv(all, out)
    out
  }
  
  # --- combine ------------------------------------------------------------
  bind_many("plot.csv")
  bind_many("cond.csv")
  bind_many("tree.csv")
  if (length(fs::dir_ls(root_states, glob = "*/seedling.csv"))) bind_many("seedling.csv")
  
  bind_many("pop_eval_filtered.csv", .distinct = TRUE, .keys = c("EVALID"))
  bind_many("pop_eval_typ.csv",      .distinct = TRUE, .keys = NULL)  # schema varies; keep all
  
  # PPSA is only written if you chose to output it; harmless if missing
  bind_many("pop_plot_stratum_assgn.csv", .distinct = TRUE, .keys = c("EVALID","PLT_CN","STRATUM_CN"))
  
  # Normalize keys then de-dup
  bind_many("pop_stratum.csv",    .distinct = TRUE, .keys = c("STRATUM_CN"),   .norm = norm_pop_stratum)
  bind_many("pop_estn_unit.csv",  .distinct = TRUE, .keys = c("ESTN_UNIT_CN"), .norm = norm_pop_estn_unit)
  
  message("✔ Region-wide FIA slices written to: ", out_region)
  invisible(out_region)
}

# ----------------------------- PUBLIC ENTRY ------------------------------------
fia_pull_states <- function(project_dir = ".",
                            cfg_path = "configs/fia_states.yml") {
  if (!fs::file_exists(cfg_path)) stop("Missing config: ", cfg_path)
  cfg <- yaml::read_yaml(cfg_path)
  if (!isTRUE(cfg$use_states)) stop("configs/fia_states.yml exists but 'use_states: true' not set.")
  stopifnot(length(cfg$states) > 0)
  
  years <- if (!is.null(cfg$years) && length(cfg$years) == 2) seq(cfg$years[[1]], cfg$years[[2]]) else 2008:2022
  eval_type <- cfg$eval_type %||% "EXPN"
  
  for (rec in cfg$states) {
    abbr <- rec$abbr; url <- rec$url
    message(glue("▶ Processing {abbr} …"))
    pull_one_state(project_dir, abbr, url, years, eval_type,
                   overwrite_zip = isTRUE(cfg$overwrite_zip %||% FALSE),
                   overwrite_unzip = isTRUE(cfg$overwrite_unzip %||% FALSE),
                   keep_zip        = isTRUE(cfg$keep_zip %||% FALSE))
  }
  
  combine_states_to_region(project_dir)
}

# CLI
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  fia_pull_states(".","configs/fia_states.yml")
}