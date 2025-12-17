# =============================================================================
# R/mod_maps.R
# Maps tab module - Hex biomass explorer + prediction rasters
# =============================================================================

maps_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    div(
      class = "container-fluid px-0",
      style = "height: calc(100vh - 70px);",
      
      layout_sidebar(
        sidebar = sidebar(
          width = 300,
          title = "Layer Controls",
          open = TRUE,
          
          # Hex Biomass Section
          h6(class = "text-primary mt-2", icon("th"), " Hex Biomass Results"),
          
          selectInput(
            ns("hex_scale"),
            "Hex Scale:",
            choices = c(
              "100 ha" = "100ha",
              "500 ha" = "500ha",
              "1 kha" = "1kha", 
              "5 kha" = "5kha",
              "10 kha" = "10kha",
              "50 kha" = "50kha",
              "64 kha (FIA)" = "64kha",
              "100 kha" = "100kha"
            ),
            selected = "10kha"
          ),
          
          selectInput(
            ns("hex_variable"),
            "Display Variable:",
            choices = c(
              "Mean Biomass" = "mean_biomass",
              "Positional SD" = "positional_sd",
              "Sampling SE" = "sampling_se",
              "Total SD" = "total_sd",
              "Sample Size (n)" = "n_plots",
              "Positional Fraction" = "positional_fraction"
            ),
            selected = "mean_biomass"
          ),
          
          selectInput(
            ns("hex_year"),
            "Year:",
            choices = c("All Years (Mean)" = "all", 
                        "2020" = "2020", "2021" = "2021", "2022" = "2022", 
                        "2023" = "2023", "2024" = "2024"),
            selected = "all"
          ),
          
          hr(),
          
          # Prediction Rasters Section
          h6(class = "text-info", icon("layer-group"), " Prediction Surfaces"),
          
          checkboxGroupInput(
            ns("raster_layers"),
            NULL,
            choices = c(
              "Biomass (FIA-trained)" = "pred_fia",
              "Biomass (NEFIN-trained)" = "pred_nefin",
              "Prediction Difference" = "pred_diff"
            ),
            selected = NULL
          ),
          
          hr(),
          
          # Plot Points Section
          h6(class = "text-success", icon("map-marker-alt"), " Plot Locations"),
          
          checkboxGroupInput(
            ns("plot_layers"),
            NULL,
            choices = c(
              "FIA Plots" = "fia",
              "NEFIN Plots" = "nefin"
            ),
            selected = NULL
          ),
          
          hr(),
          
          # Display Options
          h6(class = "text-muted", icon("sliders-h"), " Display Options"),
          
          sliderInput(
            ns("opacity"),
            "Layer Opacity:",
            min = 0.1, max = 1, value = 0.7, step = 0.1
          ),
          
          checkboxInput(ns("show_states"), "Show state boundaries", value = TRUE),
          
          hr(),
          
          # Actions
          actionButton(ns("reset_view"), "Reset View", 
                       class = "btn-secondary btn-sm w-100 mb-2",
                       icon = icon("sync")),
          
          # Info panel
          div(
            class = "alert alert-light mt-3 small",
            id = ns("layer_info"),
            uiOutput(ns("layer_description"))
          )
        ),
        
        # Main map area
        div(
          style = "position: relative; height: 100%;",
          
          leafletOutput(ns("main_map"), height = "100%") %>%
            withSpinner(type = 6, color = "#27ae60"),
          
          # Click info panel (floating)
          absolutePanel(
            id = ns("click_panel"),
            class = "card shadow",
            style = "background: white; padding: 0; width: 320px; max-height: 300px; overflow-y: auto;",
            bottom = 20,
            left = 20,
            draggable = TRUE,
            
            div(
              class = "card-header bg-light py-2",
              h6(class = "mb-0", icon("info-circle"), " Hex Details")
            ),
            div(
              class = "card-body p-2",
              uiOutput(ns("click_info"))
            )
          )
        )
      )
    )
  )
}

maps_server <- function(id, app_data, config) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive for current hex data
    current_hex_data <- reactive({
      req(input$hex_scale)
      
      scale <- input$hex_scale
      
      # Load hex grid geometry
      hex_path <- config$hex_grids[[scale]]
      results_path <- config$hex_results[[scale]]
      
      if (is.null(hex_path) || !file.exists(hex_path)) {
        message(sprintf("Hex grid not found: %s", hex_path))
        return(NULL)
      }
      
      # Read hex geometry
      hex_sf <- tryCatch({
        sf::st_read(hex_path, quiet = TRUE)
      }, error = function(e) {
        message(sprintf("Error reading hex: %s", e$message))
        return(NULL)
      })
      
      if (is.null(hex_sf)) return(NULL)
      
      # Transform to WGS84 for leaflet
      if (sf::st_crs(hex_sf)$epsg != 4326) {
        hex_sf <- sf::st_transform(hex_sf, 4326)
      }
      
      # Find the hex ID column in geometry
      hex_id_col <- intersect(names(hex_sf), c("hex_id", "HEX_ID", "id", "ID", "hexid"))[1]
      if (!is.na(hex_id_col) && hex_id_col != "hex_id") {
        names(hex_sf)[names(hex_sf) == hex_id_col] <- "hex_id"
      }
      
      # Load results CSV if available
      if (!is.null(results_path) && file.exists(results_path)) {
        message(sprintf("Loading hex results from: %s", results_path))
        results <- read_csv(results_path, show_col_types = FALSE)
        
        message(sprintf("Results columns: %s", paste(names(results), collapse = ", ")))
        
        # Find the hex ID column in results
        results_id_col <- intersect(names(results), c("hex_id", "HEX_ID", "id", "ID", "hexid"))[1]
        if (!is.na(results_id_col) && results_id_col != "hex_id") {
          names(results)[names(results) == results_id_col] <- "hex_id"
        }
        
        # Standardize column names - map from results CSV to expected names
        # Results CSV has: hex_id, year_label, window, source, mean, se, n_plots, positional_sd, n_reps, total_sd
        col_mapping <- c(
          "mean" = "mean_biomass",
          "se" = "sampling_se", 
          "positional_sd" = "positional_sd",
          "total_sd" = "total_sd",
          "n_plots" = "n_plots"
        )
        
        for (old_name in names(col_mapping)) {
          if (old_name %in% names(results)) {
            names(results)[names(results) == old_name] <- col_mapping[old_name]
          }
        }
        
        # Filter by year if specified
        if (input$hex_year != "all" && "year_label" %in% names(results)) {
          results <- results %>% filter(year_label == as.integer(input$hex_year))
        }
        
        # Aggregate if multiple years/windows and showing "all"
        if (input$hex_year == "all" && nrow(results) > 0) {
          # Check if we have multiple rows per hex
          if (any(c("year_label", "window") %in% names(results))) {
            message("Aggregating across years/windows...")
            results <- results %>%
              group_by(hex_id) %>%
              summarize(
                mean_biomass = mean(mean_biomass, na.rm = TRUE),
                positional_sd = mean(positional_sd, na.rm = TRUE),
                sampling_se = mean(sampling_se, na.rm = TRUE),
                total_sd = mean(total_sd, na.rm = TRUE),
                n_plots = mean(n_plots, na.rm = TRUE),
                .groups = "drop"
              )
          }
        }
        
        # Calculate positional fraction
        if (all(c("positional_sd", "sampling_se") %in% names(results))) {
          results <- results %>%
            mutate(
              positional_fraction = positional_sd^2 / (positional_sd^2 + sampling_se^2 + 1e-10)
            )
        }
        
        message(sprintf("Joining %d hex geometries with %d result rows", nrow(hex_sf), nrow(results)))
        
        # Join geometry with results
        hex_sf <- hex_sf %>%
          left_join(results, by = "hex_id")
        
        message(sprintf("After join: %d rows with mean_biomass data: %d", 
                        nrow(hex_sf), sum(!is.na(hex_sf$mean_biomass))))
      } else {
        message(sprintf("Results file not found: %s", results_path))
      }
      
      hex_sf
    })
    
    # Layer description
    output$layer_description <- renderUI({
      var <- input$hex_variable
      
      desc <- switch(var,
        "mean_biomass" = "Mean aboveground live biomass (Mg/ha) within each hex",
        "positional_sd" = "Standard deviation from coordinate fuzzing uncertainty",
        "sampling_se" = "Standard error from plot sampling variability", 
        "total_sd" = "Combined uncertainty (positional + sampling)",
        "n_plots" = "Number of FIA plots within the hex",
        "positional_fraction" = "Fraction of variance due to positional uncertainty (0-1)",
        "Select a variable"
      )
      
      tags$p(class = "mb-0", desc)
    })
    
    # Click data storage
    click_data <- reactiveVal(NULL)
    
    # Main map
    output$main_map <- renderLeaflet({
      
      map <- leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
        addProviderTiles(providers$CartoDB.Positron, group = "Light") %>%
        addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
        addProviderTiles(providers$OpenTopoMap, group = "Terrain") %>%
        setView(lng = -72.5, lat = 44, zoom = 7) %>%
        addScaleBar(position = "bottomleft") %>%
        addLayersControl(
          baseGroups = c("Light", "Satellite", "Terrain"),
          options = layersControlOptions(collapsed = TRUE)
        )
      
      # Add state boundaries if available
      data <- app_data()
      if (input$show_states && !is.null(data$states)) {
        map <- map %>%
          addPolygons(
            data = data$states,
            fillColor = "transparent",
            color = "#333",
            weight = 1.5,
            opacity = 0.6,
            group = "States"
          )
      }
      
      map
    })
    
    # Update hex layer when scale/variable/year changes
    observe({
      hex_sf <- current_hex_data()
      
      proxy <- leafletProxy("main_map", session)
      proxy %>% clearGroup("Hex")
      
      if (is.null(hex_sf)) return()
      
      var <- input$hex_variable
      
      # Check if variable exists
      if (!var %in% names(hex_sf)) {
        message(sprintf("Variable %s not in hex data", var))
        return()
      }
      
      values <- hex_sf[[var]]
      
      if (all(is.na(values))) return()
      
      # Choose palette based on variable
      pal <- if (var == "mean_biomass") {
        colorNumeric("YlGn", domain = values, na.color = "transparent")
      } else if (var %in% c("positional_sd", "sampling_se", "total_sd")) {
        colorNumeric("OrRd", domain = values, na.color = "transparent")
      } else if (var == "positional_fraction") {
        colorNumeric("RdYlBu", domain = c(0, 1), na.color = "transparent", reverse = TRUE)
      } else {
        colorNumeric("viridis", domain = values, na.color = "transparent")
      }
      
      # Build popup content
      popup_content <- sprintf(
        "<strong>Hex:</strong> %s<br>
         <strong>%s:</strong> %.2f<br>
         <strong>N plots:</strong> %s",
        hex_sf$hex_id,
        var,
        values,
        if ("n_plots" %in% names(hex_sf)) round(hex_sf$n_plots, 0) else "N/A"
      )
      
      proxy %>%
        addPolygons(
          data = hex_sf,
          fillColor = ~pal(values),
          fillOpacity = input$opacity,
          color = "white",
          weight = 0.5,
          opacity = 0.8,
          group = "Hex",
          layerId = ~hex_id,
          popup = popup_content,
          highlightOptions = highlightOptions(
            weight = 2,
            color = "#333",
            fillOpacity = 0.9,
            bringToFront = TRUE
          )
        ) %>%
        addLegend(
          layerId = "hex_legend",
          position = "bottomright",
          pal = pal,
          values = values,
          title = switch(var,
            "mean_biomass" = "Biomass (Mg/ha)",
            "positional_sd" = "Positional SD",
            "sampling_se" = "Sampling SE",
            "total_sd" = "Total SD",
            "n_plots" = "N Plots",
            "positional_fraction" = "Pos. Fraction",
            var
          ),
          opacity = 1
        )
    })
    
    # Update raster layers
    observe({
      proxy <- leafletProxy("main_map", session)
      proxy %>% 
        clearGroup("Raster FIA") %>%
        clearGroup("Raster NEFIN") %>%
        clearGroup("Raster Diff")
      
      # Add prediction rasters if selected
      if ("pred_fia" %in% input$raster_layers && file.exists(config$pred_fia_tif)) {
        tryCatch({
          r <- terra::rast(config$pred_fia_tif)
          # Project to WGS84 if needed
          if (!terra::same.crs(r, "EPSG:4326")) {
            r <- terra::project(r, "EPSG:4326", method = "bilinear")
          }
          
          # Sample for web display
          if (terra::ncell(r) > 500000) {
            r <- terra::aggregate(r, fact = ceiling(sqrt(terra::ncell(r) / 500000)))
          }
          
          vals <- terra::values(r, na.rm = TRUE)
          pal <- colorNumeric("YlGn", domain = range(vals), na.color = "transparent")
          
          proxy %>%
            addRasterImage(
              raster::raster(r),
              colors = pal,
              opacity = input$opacity * 0.8,
              group = "Raster FIA"
            )
        }, error = function(e) message(sprintf("Error loading FIA raster: %s", e$message)))
      }
      
      if ("pred_nefin" %in% input$raster_layers && file.exists(config$pred_nefin_tif)) {
        tryCatch({
          r <- terra::rast(config$pred_nefin_tif)
          if (!terra::same.crs(r, "EPSG:4326")) {
            r <- terra::project(r, "EPSG:4326", method = "bilinear")
          }
          if (terra::ncell(r) > 500000) {
            r <- terra::aggregate(r, fact = ceiling(sqrt(terra::ncell(r) / 500000)))
          }
          
          vals <- terra::values(r, na.rm = TRUE)
          pal <- colorNumeric("YlGn", domain = range(vals), na.color = "transparent")
          
          proxy %>%
            addRasterImage(
              raster::raster(r),
              colors = pal,
              opacity = input$opacity * 0.8,
              group = "Raster NEFIN"
            )
        }, error = function(e) message(sprintf("Error loading NEFIN raster: %s", e$message)))
      }
      
      if ("pred_diff" %in% input$raster_layers && file.exists(config$pred_diff_tif)) {
        tryCatch({
          r <- terra::rast(config$pred_diff_tif)
          if (!terra::same.crs(r, "EPSG:4326")) {
            r <- terra::project(r, "EPSG:4326", method = "bilinear")
          }
          if (terra::ncell(r) > 500000) {
            r <- terra::aggregate(r, fact = ceiling(sqrt(terra::ncell(r) / 500000)))
          }
          
          vals <- terra::values(r, na.rm = TRUE)
          max_abs <- max(abs(range(vals, na.rm = TRUE)))
          pal <- colorNumeric("RdBu", domain = c(-max_abs, max_abs), na.color = "transparent")
          
          proxy %>%
            addRasterImage(
              raster::raster(r),
              colors = pal,
              opacity = input$opacity * 0.8,
              group = "Raster Diff"
            )
        }, error = function(e) message(sprintf("Error loading diff raster: %s", e$message)))
      }
    })
    
    # Update plot points
    observe({
      data <- app_data()
      proxy <- leafletProxy("main_map", session)
      
      proxy %>%
        clearGroup("FIA") %>%
        clearGroup("NEFIN")
      
      if ("fia" %in% input$plot_layers && !is.null(data$fia)) {
        fia <- data$fia
        
        # Find coordinate columns
        lon_col <- names(fia)[grepl("^lon|longitude", names(fia), ignore.case = TRUE)][1]
        lat_col <- names(fia)[grepl("^lat|latitude", names(fia), ignore.case = TRUE)][1]
        
        if (!is.na(lon_col) && !is.na(lat_col)) {
          fia_sample <- fia %>%
            filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]])) %>%
            sample_n(min(2000, n()))
          
          proxy %>%
            addCircleMarkers(
              data = fia_sample,
              lng = ~get(lon_col),
              lat = ~get(lat_col),
              radius = 3,
              color = "#e74c3c",
              fillOpacity = 0.6,
              stroke = FALSE,
              group = "FIA",
              popup = ~paste0("FIA Plot<br>Biomass: ", 
                              if("aglb_Mg_per_ha" %in% names(.)) round(aglb_Mg_per_ha, 1) else "N/A",
                              " Mg/ha")
            )
        }
      }
      
      if ("nefin" %in% input$plot_layers && !is.null(data$nefin)) {
        nefin <- data$nefin
        
        lon_col <- names(nefin)[grepl("^lon|longitude", names(nefin), ignore.case = TRUE)][1]
        lat_col <- names(nefin)[grepl("^lat|latitude", names(nefin), ignore.case = TRUE)][1]
        
        if (!is.na(lon_col) && !is.na(lat_col)) {
          nefin_sample <- nefin %>%
            filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]])) %>%
            sample_n(min(2000, n()))
          
          proxy %>%
            addCircleMarkers(
              data = nefin_sample,
              lng = ~get(lon_col),
              lat = ~get(lat_col),
              radius = 3,
              color = "#27ae60",
              fillOpacity = 0.6,
              stroke = FALSE,
              group = "NEFIN",
              popup = ~paste0("NEFIN Plot<br>Biomass: ",
                              if("aglb_Mg_per_ha" %in% names(.)) round(aglb_Mg_per_ha, 1) else "N/A",
                              " Mg/ha")
            )
        }
      }
    })
    
    # Handle hex clicks
    observeEvent(input$main_map_shape_click, {
      click <- input$main_map_shape_click
      click_data(click)
    })
    
    # Click info display
    output$click_info <- renderUI({
      click <- click_data()
      hex_sf <- current_hex_data()
      
      if (is.null(click) || is.null(hex_sf)) {
        return(tags$p(class = "text-muted", "Click on a hex to see details"))
      }
      
      hex_id <- click$id
      row <- hex_sf %>% filter(hex_id == !!hex_id)
      
      if (nrow(row) == 0) {
        return(tags$p(class = "text-muted", "Hex not found"))
      }
      
      # Build info display
      items <- list(
        tags$h6(class = "text-primary", paste("Hex:", hex_id)),
        tags$hr(class = "my-2")
      )
      
      # Add available metrics
      metrics <- c(
        "mean_biomass" = "Mean Biomass (Mg/ha)",
        "positional_sd" = "Positional SD",
        "sampling_se" = "Sampling SE",
        "total_sd" = "Total SD",
        "n_plots" = "N Plots",
        "positional_fraction" = "Positional Fraction"
      )
      
      for (col in names(metrics)) {
        if (col %in% names(row) && !is.na(row[[col]][1])) {
          val <- row[[col]][1]
          formatted <- if (col == "n_plots") round(val, 0) else round(val, 2)
          items <- c(items, list(
            tags$p(
              class = "mb-1",
              tags$strong(paste0(metrics[col], ": ")),
              formatted
            )
          ))
        }
      }
      
      do.call(tags$div, items)
    })
    
    # Reset view
    observeEvent(input$reset_view, {
      leafletProxy("main_map", session) %>%
        setView(lng = -72.5, lat = 44, zoom = 7)
    })
    
  })
}
