# R/01_init_project.R
# Initialize project directory structure

suppressPackageStartupMessages({ library(fs) })

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
    configs = fs::path(project_dir, "configs"),
    runs = fs::path(project_dir, "runs")
  )
  lapply(paths, function(p) if (!fs::dir_exists(p)) fs::dir_create(p, recurse=TRUE))
  invisible(paths)
}

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  dir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
  init_project(dir)
  cat("✔ Project initialized at:", normalizePath(dir), "\n")
}