# =============================================================================
# NEFIN vs FIA Dashboard - Main Application
# VERSION: 4.0 - 2024-12-17
# Updated to use actual project data
# =============================================================================

library(shiny)
library(bslib)
library(leaflet)
library(leaflet.extras)
library(plotly)
library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(viridis)
library(scales)
library(DT)
library(htmltools)
library(shinycssloaders)

# Source modules
source("R/utils.R")
source("R/mod_overview.R")
source("R/mod_data.R")
source("R/mod_fuzzing.R")
source("R/mod_models.R")
source("R/mod_maps.R")

# =============================================================================
# Configuration - UPDATE THESE PATHS to match your local setup
# =============================================================================
# Base path - when running from dashboard folder, we need to go up one level
# to reach the project root where data/ and runs/ are located
BASE_PATH <- ".."  # Goes up from dashboard/ to project root

config <- list(
  # Core plot data
  fia_complete = file.path(BASE_PATH, "data/processed/fia_complete.csv"),
  nefin_processed = file.path(BASE_PATH, "data/processed/nefin_processed.csv"),
  
  # Climate/NDVI at plots
  fia_climate = file.path(BASE_PATH, "data/processed/climate_at_plots/fia_climate.csv"),
  nefin_climate = file.path(BASE_PATH, "data/processed/climate_at_plots/nefin_climate.csv"),
  fia_ndvi = file.path(BASE_PATH, "data/processed/ndvi_at_plots/fia_ndvi.csv"),
  nefin_ndvi = file.path(BASE_PATH, "data/processed/ndvi_at_plots/nefin_ndvi.csv"),
  
  # Boundaries
  states = file.path(BASE_PATH, "data/boundaries/states_5070.geojson"),
  
  # Hex grids
  hex_grids = list(
    `100ha` = file.path(BASE_PATH, "data/hex/hex_grid_100ha.geojson"),
    `500ha` = file.path(BASE_PATH, "data/hex/hex_grid_500ha.geojson"),
    `1kha` = file.path(BASE_PATH, "data/hex/hex_grid_1kha.geojson"),
    `5kha` = file.path(BASE_PATH, "data/hex/hex_grid_5kha.geojson"),
    `10kha` = file.path(BASE_PATH, "data/hex/hex_grid_10kha.geojson"),
    `50kha` = file.path(BASE_PATH, "data/hex/hex_grid_50kha.geojson"),
    `64kha` = file.path(BASE_PATH, "data/hex/hex_grid_6kac.geojson"),
    `100kha` = file.path(BASE_PATH, "data/hex/hex_grid_100kha.geojson")
  ),
  
  # Hex biomass results (from runs/)
  hex_results = list(
    `100ha` = file.path(BASE_PATH, "runs/2025-12-15_aglb_100ha_W5y/hex_aglb_results.csv"),
    `500ha` = file.path(BASE_PATH, "runs/2025-12-15_aglb_500ha_W5y/hex_aglb_results.csv"),
    `1kha` = file.path(BASE_PATH, "runs/2025-12-15_aglb_1kha_W5y/hex_aglb_results.csv"),
    `5kha` = file.path(BASE_PATH, "runs/2025-12-15_aglb_5kha_W5y/hex_aglb_results.csv"),
    `10kha` = file.path(BASE_PATH, "runs/2025-12-15_aglb_10kha_W5y/hex_aglb_results.csv"),
    `50kha` = file.path(BASE_PATH, "runs/2025-12-15_aglb_50kha_W5y/hex_aglb_results.csv"),
    `64kha` = file.path(BASE_PATH, "runs/2025-12-15_aglb_64kha_W5y/hex_aglb_results.csv"),
    `100kha` = file.path(BASE_PATH, "runs/2025-12-15_aglb_100kha_W5y/hex_aglb_results.csv")
  ),
  
  # Consolidated results
  consolidated = file.path(BASE_PATH, "runs/consolidated_20251211_120654"),
  fia_all_scales = file.path(BASE_PATH, "runs/consolidated_20251211_120654/fia_all_scales.csv"),
  fia_nefin_comparison = file.path(BASE_PATH, "runs/consolidated_20251211_120654/fia_nefin_comparison_all_scales.csv"),
  fia_summary_by_scale = file.path(BASE_PATH, "runs/consolidated_20251211_120654/fia_summary_by_scale.csv"),
  
  # Model comparison results
  model_comparison = file.path(BASE_PATH, "runs/spatial_model_comparison/model_comparison_summary.csv"),
  holdout_results = file.path(BASE_PATH, "runs/spatial_model_comparison/holdout_prediction_results.csv"),
  
  # Fuzzing effect analysis
  covariate_uncertainty = file.path(BASE_PATH, "runs/fuzzing_effect_analysis/covariate_uncertainty_summary.csv"),
  prediction_uncertainty = file.path(BASE_PATH, "runs/fuzzing_effect_analysis/prediction_uncertainty_by_plot.csv"),
  
  # Prediction rasters
  pred_fia_tif = file.path(BASE_PATH, "runs/prediction_maps/modis_biomass_fia.tif"),
  pred_nefin_tif = file.path(BASE_PATH, "runs/prediction_maps/modis_biomass_nefin.tif"),
  pred_diff_tif = file.path(BASE_PATH, "runs/prediction_maps/modis_prediction_difference.tif"),
  
  # NDVI rasters
  modis_ndvi = file.path(BASE_PATH, "data/processed/ndvi/modis"),
  s2_ndvi = file.path(BASE_PATH, "data/processed/ndvi/s2/S2_NDVI_10m_2020_2025.tif"),
  
  # PRISM climate
  prism = file.path(BASE_PATH, "data/processed/prism"),
  
  # MC jitter library
  mc_replicates = file.path(BASE_PATH, "data/processed/mc_jitter_library/replicates")
)

# =============================================================================
# Data Loading Functions
# =============================================================================

load_app_data <- function(config) {
  data <- list()
  
  # Load FIA data
  if (file.exists(config$fia_complete)) {
    message("Loading FIA data...")
    data$fia <- read_csv(config$fia_complete, show_col_types = FALSE)
    message(sprintf("  Loaded %d FIA plots", nrow(data$fia)))
  }
  
  # Load NEFIN data
  if (file.exists(config$nefin_processed)) {
    message("Loading NEFIN data...")
    data$nefin <- read_csv(config$nefin_processed, show_col_types = FALSE)
    message(sprintf("  Loaded %d NEFIN plots", nrow(data$nefin)))
  }
  
  # Load state boundaries
  if (file.exists(config$states)) {
    message("Loading state boundaries...")
    data$states <- st_read(config$states, quiet = TRUE)
    # Transform to WGS84 for leaflet
    if (st_crs(data$states)$epsg != 4326) {
      data$states <- st_transform(data$states, 4326)
    }
  }
  
  # Load model comparison results
  if (file.exists(config$model_comparison)) {
    message("Loading model comparison results...")
    data$model_comparison <- read_csv(config$model_comparison, show_col_types = FALSE)
  }
  
  # Load holdout results
  if (file.exists(config$holdout_results)) {
    message("Loading holdout prediction results...")
    data$holdout_results <- read_csv(config$holdout_results, show_col_types = FALSE)
  }
  
  # Load covariate uncertainty
  if (file.exists(config$covariate_uncertainty)) {
    message("Loading covariate uncertainty...")
    data$covariate_uncertainty <- read_csv(config$covariate_uncertainty, show_col_types = FALSE)
  }
  
  # Load prediction uncertainty
  if (file.exists(config$prediction_uncertainty)) {
    message("Loading prediction uncertainty...")
    data$prediction_uncertainty <- read_csv(config$prediction_uncertainty, show_col_types = FALSE)
  }
  
  # Load FIA-NEFIN comparison
  if (file.exists(config$fia_nefin_comparison)) {
    message("Loading FIA-NEFIN comparison...")
    data$fia_nefin_comparison <- read_csv(config$fia_nefin_comparison, show_col_types = FALSE)
  }
  
  # Load FIA summary by scale
  if (file.exists(config$fia_summary_by_scale)) {
    message("Loading FIA summary by scale...")
    data$fia_summary_by_scale <- read_csv(config$fia_summary_by_scale, show_col_types = FALSE)
  }
  
  # Store config for later use
  data$config <- config
  
  message("Data loading complete!")
  return(data)
}

# =============================================================================
# UI
# =============================================================================

ui <- page_navbar(
  title = div(
    img(src = "logo.svg", height = "30px", class = "me-2"),
    "NEFIN vs FIA Analysis"
  ),
  id = "main_nav",
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#27ae60",
    secondary = "#95a5a6",
    success = "#2ecc71",
    info = "#3498db",
    warning = "#f39c12",
    danger = "#e74c3c",
    base_font = font_google("Open Sans"),
    heading_font = font_google("Montserrat")
  ),
  
  header = tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),
  
  nav_panel(
    title = "Overview",
    icon = icon("home"),
    overview_ui("overview")
  ),
  
  nav_panel(
    title = "Data",
    icon = icon("database"),
    data_ui("data")
  ),
  
  nav_panel(
    title = "Fuzzing Effects",
    icon = icon("random"),
    fuzzing_ui("fuzzing")
  ),
  
  nav_panel(
    title = "Models",
    icon = icon("chart-bar"),
    models_ui("models")
  ),
  
  nav_panel(
    title = "Maps",
    icon = icon("map"),
    maps_ui("maps")
  ),
  
  nav_spacer(),
  
  nav_item(
    tags$a(
      icon("github"),
      href = "https://github.com/",
      target = "_blank",
      class = "nav-link"
    )
  )
)

# =============================================================================
# Server
# =============================================================================

server <- function(input, output, session) {
  
  # Load data once at startup
  app_data <- reactiveVal()
  
  # Show loading modal
  showModal(modalDialog(
    title = "Loading Data",
    div(
      class = "text-center",
      div(class = "spinner-border text-primary", role = "status"),
      p(class = "mt-3", "Loading datasets and results...")
    ),
    footer = NULL,
    easyClose = FALSE
  ))
  
  # Load data
  observe({
    data <- load_app_data(config)
    app_data(data)
    removeModal()
  })
  
  # Initialize modules
  overview_server("overview", app_data, config)
  data_server("data", app_data, config)
  fuzzing_server("fuzzing", app_data, config)
  models_server("models", app_data, config)
  maps_server("maps", app_data, config)
}

# =============================================================================
# Run App
# =============================================================================

shinyApp(ui = ui, server = server)
