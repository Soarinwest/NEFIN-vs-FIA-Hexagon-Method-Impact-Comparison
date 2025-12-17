# =============================================================================
# R/utils.R
# Helper functions for the NEFIN vs FIA Dashboard
# =============================================================================

# Color palettes
biomass_palette <- function(n = 100) {
  colorRampPalette(c("#f7fcf5", "#c7e9c0", "#74c476", "#238b45", "#00441b"))(n)
}

difference_palette <- function(n = 100) {
  colorRampPalette(c("#2166ac", "#67a9cf", "#f7f7f7", "#ef8a62", "#b2182b"))(n)
}

uncertainty_palette <- function(n = 100) {
  colorRampPalette(c("#fff5f0", "#fcbba1", "#fb6a4a", "#cb181d", "#67000d"))(n)
}

ndvi_palette <- function(n = 100) {
  colorRampPalette(c("#ffffe5", "#d9f0a3", "#78c679", "#238443", "#004529"))(n)
}

# Create demo FIA data for testing without real data
create_demo_fia_data <- function(n = 5000) {
  set.seed(42)
  
  # Generate random points in northeastern US
  lon_range <- c(-79, -67)
  lat_range <- c(41, 47)
  
  data.frame(
    CN = sprintf("FIA_%06d", 1:n),
    STATECD = sample(c(9, 23, 25, 33, 36, 44, 50), n, replace = TRUE),
    lon_original = runif(n, lon_range[1], lon_range[2]),
    lat_original = runif(n, lat_range[1], lat_range[2]),
    MEASYEAR = sample(2015:2024, n, replace = TRUE),
    aglb_Mg_per_ha = rlnorm(n, meanlog = 4.5, sdlog = 0.6),
    ndvi_modis = runif(n, 0.4, 0.9),
    tmean = runif(n, 4, 12),
    ppt = runif(n, 800, 1400)
  )
}

# Create demo NEFIN data for testing
create_demo_nefin_data <- function(n = 3000) {
  set.seed(123)
  
  lon_range <- c(-79, -67)
  lat_range <- c(41, 47)
  
  data.frame(
    CN = sprintf("NEFIN_%06d", 1:n),
    PLOT_ID = sprintf("P%05d", 1:n),
    lon_public = runif(n, lon_range[1], lon_range[2]),
    lat_public = runif(n, lat_range[1], lat_range[2]),
    MEASYEAR = sample(2015:2024, n, replace = TRUE),
    aglb_Mg_per_ha = rlnorm(n, meanlog = 4.6, sdlog = 0.55),
    ndvi_modis = runif(n, 0.45, 0.88),
    tmean = runif(n, 4.5, 11.5),
    ppt = runif(n, 820, 1380)
  )
}

# Format large numbers with commas
fmt_num <- function(x, digits = 1) {
  format(round(x, digits), big.mark = ",", scientific = FALSE)
}

# Format percentages
fmt_pct <- function(x, digits = 1) {
  paste0(round(x, digits), "%")
}

# Standardize scale names (fia -> 64kha)
standardize_scale_name <- function(x) {
  x <- as.character(x)
  x[tolower(x) == "fia"] <- "64kha"
  x
}

# Get ordered scale levels
get_scale_order <- function() {
  c("100ha", "500ha", "1kha", "5kha", "10kha", "50kha", "64kha", "100kha")
}

# Order scales in a dataframe
order_scales <- function(df, scale_col = "grid_scale") {
  if (scale_col %in% names(df)) {
    df[[scale_col]] <- standardize_scale_name(df[[scale_col]])
    scale_order <- get_scale_order()
    existing <- unique(df[[scale_col]])
    ordered_levels <- scale_order[scale_order %in% existing]
    extra <- setdiff(existing, scale_order)
    if (length(extra) > 0) ordered_levels <- c(ordered_levels, sort(extra))
    df[[scale_col]] <- factor(df[[scale_col]], levels = ordered_levels, ordered = TRUE)
  }
  df
}

# Create value box card
value_box_custom <- function(value, subtitle, icon_name = "chart-line", color = "primary") {
  card(
    card_body(
      class = paste0("bg-", color, " text-white text-center p-3"),
      div(
        style = "font-size: 2.5rem; font-weight: bold;",
        value
      ),
      div(
        style = "font-size: 0.9rem; opacity: 0.9;",
        icon(icon_name, style = "margin-right: 5px;"),
        subtitle
      )
    )
  )
}

# Create a stat card with comparison
stat_card <- function(title, value, change = NULL, change_positive = TRUE) {
  change_html <- if (!is.null(change)) {
    color <- if (change_positive) "success" else "danger"
    icon_name <- if (change_positive) "arrow-up" else "arrow-down"
    tags$span(
      class = paste0("text-", color),
      icon(icon_name),
      change
    )
  }
  
  card(
    card_body(
      h6(class = "text-muted mb-1", title),
      div(
        class = "d-flex align-items-baseline",
        tags$span(class = "h3 mb-0", value),
        if (!is.null(change_html)) {
          tags$span(class = "ms-2 small", change_html)
        }
      )
    )
  )
}

# Safe file.exists that handles NULLs
safe_file_exists <- function(path) {
  if (is.null(path)) return(FALSE)
  file.exists(path)
}

# Calculate summary statistics for a numeric vector
calc_summary_stats <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NULL)
  
  list(
    n = length(x),
    mean = mean(x),
    median = median(x),
    sd = sd(x),
    min = min(x),
    max = max(x),
    q25 = quantile(x, 0.25),
    q75 = quantile(x, 0.75)
  )
}

# Create a simple leaflet base map
create_base_map <- function() {
  leaflet() %>%
    addProviderTiles(
      providers$CartoDB.Positron,
      group = "Light"
    ) %>%
    addProviderTiles(
      providers$Esri.WorldImagery,
      group = "Satellite"
    ) %>%
    addProviderTiles(
      providers$OpenTopoMap,
      group = "Terrain"
    ) %>%
    addLayersControl(
      baseGroups = c("Light", "Satellite", "Terrain"),
      options = layersControlOptions(collapsed = TRUE)
    ) %>%
    setView(lng = -73, lat = 44, zoom = 6)
}

# Add uncertainty circles to a leaflet map
add_uncertainty_circles <- function(map, data, lon_col, lat_col, radius_m = 1609) {
  if (nrow(data) == 0) return(map)
  
  # Sample if too many points
  if (nrow(data) > 500) {
    data <- data[sample(1:nrow(data), 500), ]
  }
  
  map %>%
    addCircles(
      data = data,
      lng = ~get(lon_col),
      lat = ~get(lat_col),
      radius = radius_m,
      stroke = TRUE,
      color = "#e74c3c",
      weight = 1,
      fillOpacity = 0.1,
      popup = ~paste0(
        "<strong>Plot: </strong>", CN, "<br>",
        "<strong>Uncertainty radius: </strong>±1.6km"
      )
    )
}

# ggplot2 theme for dashboard
theme_dashboard <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2)),
      plot.subtitle = element_text(color = "gray50"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

# Convert hex grid to leaflet-compatible format
prepare_hex_for_leaflet <- function(hex_sf, data_df, join_col = "hex_id", value_col = "mean") {
  if (is.null(hex_sf) || is.null(data_df)) return(NULL)
  
  # Ensure hex_id is character
  hex_sf$hex_id <- as.character(hex_sf[[join_col]])
  data_df[[join_col]] <- as.character(data_df[[join_col]])
  
  # Join
  result <- hex_sf %>%
    left_join(data_df, by = join_col)
  
  # Transform to WGS84 if needed
  if (st_crs(result)$epsg != 4326) {
    result <- st_transform(result, 4326)
  }
  
  result
}
