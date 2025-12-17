# =============================================================================
# R/mod_data.R
# Data tab module - Explore FIA vs NEFIN datasets
# =============================================================================

data_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    div(
      class = "container-fluid px-4 py-3",
      
      h2(icon("database"), " The Data: FIA vs NEFIN"),
      p(class = "lead text-muted",
        "Understanding the two datasets and the fuzzing problem"),
      
      # Section 2.1: FIA vs NEFIN Comparison
      h4(class = "mt-4", "2.1 Side-by-Side Comparison"),
      
      layout_columns(
        col_widths = c(6, 6),
        
        # FIA map panel
        card(
          card_header(
            class = "bg-danger bg-opacity-25",
            icon("map-marker-alt"),
            " FIA Plots (Fuzzed ±1.6km)"
          ),
          card_body(
            leafletOutput(ns("fia_map"), height = "350px") %>%
              withSpinner(type = 6)
          ),
          card_footer(
            class = "text-center",
            uiOutput(ns("fia_stats"))
          )
        ),
        
        # NEFIN map panel
        card(
          card_header(
            class = "bg-success bg-opacity-25",
            icon("crosshairs"),
            " NEFIN Plots (True Coordinates)"
          ),
          card_body(
            leafletOutput(ns("nefin_map"), height = "350px") %>%
              withSpinner(type = 6)
          ),
          card_footer(
            class = "text-center",
            uiOutput(ns("nefin_stats"))
          )
        )
      ),
      
      # Uncertainty visualization
      card(
        class = "mt-3",
        card_header(
          class = "bg-warning bg-opacity-25",
          icon("circle-notch"),
          " Coordinate Uncertainty Visualization"
        ),
        card_body(
          p("Click on an FIA plot to see its 1.6km uncertainty radius"),
          leafletOutput(ns("uncertainty_map"), height = "300px") %>%
            withSpinner(type = 6)
        )
      ),
      
      # Section 2.1b: FIA vs NEFIN Displacement
      h4(class = "mt-4", "2.1b Actual Coordinate Displacement"),
      p(class = "text-muted", 
        "For plots in both FIA and NEFIN, we can compare the fuzzed vs true coordinates"),
      
      card(
        card_body(
          layout_sidebar(
            sidebar = sidebar(
              width = 220,
              title = "Displacement View",
              
              selectInput(
                ns("displacement_region"),
                "Focus Region:",
                choices = c(
                  "All Matched Plots" = "all",
                  "Vermont" = "vt",
                  "New Hampshire" = "nh",
                  "Maine" = "me",
                  "New York" = "ny"
                ),
                selected = "all"
              ),
              
              checkboxInput(ns("show_displacement_lines"), 
                            "Show displacement lines", value = TRUE),
              
              checkboxInput(ns("color_by_displacement"),
                            "Color by displacement distance", value = TRUE),
              
              hr(),
              
              div(
                class = "text-center p-2 bg-info bg-opacity-10 rounded mb-2",
                h6(class = "text-muted mb-1", "Mean Displacement"),
                h4(class = "text-info mb-0", textOutput(ns("mean_displacement_text")))
              ),
              div(
                class = "text-center p-2 bg-danger bg-opacity-10 rounded",
                h6(class = "text-muted mb-1", "Max Displacement"),
                h4(class = "text-danger mb-0", textOutput(ns("max_displacement_text")))
              )
            ),
            
            leafletOutput(ns("displacement_map"), height = "400px") %>%
              withSpinner(type = 6)
          )
        ),
        card_footer(
          class = "text-muted small",
          tags$strong("Red"), " = FIA published (fuzzed) | ",
          tags$strong("Green"), " = NEFIN (true) | ",
          tags$strong("Purple lines"), " = displacement vector"
        )
      ),
      
      # Displacement histogram
      layout_columns(
        col_widths = c(6, 6),
        
        card(
          card_header("Distribution of Actual Displacements"),
          card_body(
            plotlyOutput(ns("displacement_hist"), height = "280px") %>%
              withSpinner(type = 6)
          )
        ),
        
        card(
          card_header("Displacement by State"),
          card_body(
            plotlyOutput(ns("displacement_by_state"), height = "280px") %>%
              withSpinner(type = 6)
          )
        )
      ),
      
      # Section 2.2: Temporal Coverage
      h4(class = "mt-4", "2.2 Temporal Coverage"),
      
      layout_columns(
        col_widths = c(8, 4),
        
        card(
          card_body(
            plotlyOutput(ns("temporal_plot"), height = "300px") %>%
              withSpinner(type = 6)
          )
        ),
        
        card(
          card_header("Measurement Years"),
          card_body(
            tableOutput(ns("year_table"))
          )
        )
      ),
      
      # Section 2.3: Biomass Distributions
      h4(class = "mt-4", "2.3 Biomass Distributions"),
      
      layout_columns(
        col_widths = c(6, 6),
        
        card(
          card_body(
            plotlyOutput(ns("biomass_violin"), height = "350px") %>%
              withSpinner(type = 6)
          )
        ),
        
        card(
          card_body(
            plotlyOutput(ns("biomass_density"), height = "350px") %>%
              withSpinner(type = 6)
          )
        )
      ),
      
      # Section 2.4: Covariate Layers
      h4(class = "mt-4", "2.4 Environmental Covariates"),
      
      card(
        card_body(
          layout_sidebar(
            sidebar = sidebar(
              width = 250,
              title = "Layer Controls",
              
              checkboxGroupInput(
                ns("raster_layers"),
                "Raster Layers:",
                choices = c(
                  "MODIS NDVI (250m)" = "modis_ndvi",
                  "Sentinel-2 NDVI (10m)" = "s2_ndvi",
                  "Mean Temperature" = "tmean",
                  "Precipitation" = "ppt"
                ),
                selected = "modis_ndvi"
              ),
              
              hr(),
              
              checkboxGroupInput(
                ns("plot_layers"),
                "Plot Layers:",
                choices = c(
                  "FIA Plots" = "fia",
                  "NEFIN Plots" = "nefin"
                ),
                selected = c("fia", "nefin")
              ),
              
              hr(),
              
              selectInput(
                ns("color_by"),
                "Color Plots By:",
                choices = c(
                  "Biomass" = "aglb_Mg_per_ha",
                  "NDVI" = "ndvi_modis",
                  "Temperature" = "tmean",
                  "Precipitation" = "ppt"
                ),
                selected = "aglb_Mg_per_ha"
              ),
              
              hr(),
              
              sliderInput(
                ns("opacity"),
                "Layer Opacity:",
                min = 0, max = 1, value = 0.7, step = 0.1
              )
            ),
            
            leafletOutput(ns("covariate_map"), height = "500px") %>%
              withSpinner(type = 6)
          )
        )
      ),
      
      # Data summary table
      h4(class = "mt-4", "2.5 Dataset Summary"),
      
      card(
        card_body(
          DTOutput(ns("summary_table")) %>%
            withSpinner(type = 6)
        )
      )
    )
  )
}

data_server <- function(id, app_data, config) {
  moduleServer(id, function(input, output, session) {
    
    # FIA map
    output$fia_map <- renderLeaflet({
      data <- app_data()
      if (is.null(data$fia)) return(create_base_map())
      
      fia <- data$fia
      lon_col <- names(fia)[grepl("lon", names(fia), ignore.case = TRUE)][1]
      lat_col <- names(fia)[grepl("lat", names(fia), ignore.case = TRUE)][1]
      
      if (is.null(lon_col) || is.null(lat_col)) return(create_base_map())
      
      fia_sample <- fia %>%
        filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]])) %>%
        sample_n(min(2000, n()))
      
      # Color by biomass
      pal <- colorNumeric(
        palette = "YlGn",
        domain = fia_sample$aglb_Mg_per_ha,
        na.color = "gray"
      )
      
      create_base_map() %>%
        addCircleMarkers(
          data = fia_sample,
          lng = ~get(lon_col),
          lat = ~get(lat_col),
          radius = 4,
          color = ~pal(aglb_Mg_per_ha),
          fillOpacity = 0.7,
          stroke = FALSE,
          popup = ~paste0(
            "<strong>Plot:</strong> ", CN, "<br>",
            "<strong>Biomass:</strong> ", round(aglb_Mg_per_ha, 1), " Mg/ha<br>",
            "<strong>Year:</strong> ", MEASYEAR
          )
        ) %>%
        addLegend(
          position = "bottomright",
          pal = pal,
          values = fia_sample$aglb_Mg_per_ha,
          title = "Biomass (Mg/ha)",
          opacity = 0.7
        )
    })
    
    # NEFIN map
    output$nefin_map <- renderLeaflet({
      data <- app_data()
      if (is.null(data$nefin)) return(create_base_map())
      
      nefin <- data$nefin
      lon_col <- names(nefin)[grepl("lon", names(nefin), ignore.case = TRUE)][1]
      lat_col <- names(nefin)[grepl("lat", names(nefin), ignore.case = TRUE)][1]
      
      if (is.null(lon_col) || is.null(lat_col)) return(create_base_map())
      
      nefin_sample <- nefin %>%
        filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]])) %>%
        sample_n(min(2000, n()))
      
      pal <- colorNumeric(
        palette = "YlGn",
        domain = nefin_sample$aglb_Mg_per_ha,
        na.color = "gray"
      )
      
      create_base_map() %>%
        addCircleMarkers(
          data = nefin_sample,
          lng = ~get(lon_col),
          lat = ~get(lat_col),
          radius = 4,
          color = ~pal(aglb_Mg_per_ha),
          fillOpacity = 0.7,
          stroke = FALSE,
          popup = ~paste0(
            "<strong>Plot:</strong> ", CN, "<br>",
            "<strong>Biomass:</strong> ", round(aglb_Mg_per_ha, 1), " Mg/ha"
          )
        ) %>%
        addLegend(
          position = "bottomright",
          pal = pal,
          values = nefin_sample$aglb_Mg_per_ha,
          title = "Biomass (Mg/ha)",
          opacity = 0.7
        )
    })
    
    # FIA stats
    output$fia_stats <- renderUI({
      data <- app_data()
      if (is.null(data$fia)) return(NULL)
      
      fia <- data$fia
      tags$div(
        tags$strong(format(nrow(fia), big.mark = ",")), " plots | ",
        tags$strong(round(mean(fia$aglb_Mg_per_ha, na.rm = TRUE), 1)), " Mg/ha mean | ",
        "Fuzzed ±1.6km"
      )
    })
    
    # NEFIN stats
    output$nefin_stats <- renderUI({
      data <- app_data()
      if (is.null(data$nefin)) return(NULL)
      
      nefin <- data$nefin
      tags$div(
        tags$strong(format(nrow(nefin), big.mark = ",")), " plots | ",
        tags$strong(round(mean(nefin$aglb_Mg_per_ha, na.rm = TRUE), 1)), " Mg/ha mean | ",
        "True coordinates"
      )
    })
    
    # Uncertainty map
    output$uncertainty_map <- renderLeaflet({
      data <- app_data()
      if (is.null(data$fia)) return(create_base_map())
      
      fia <- data$fia
      lon_col <- names(fia)[grepl("lon", names(fia), ignore.case = TRUE)][1]
      lat_col <- names(fia)[grepl("lat", names(fia), ignore.case = TRUE)][1]
      
      if (is.null(lon_col) || is.null(lat_col)) return(create_base_map())
      
      # Sample a few plots to show uncertainty circles
      fia_sample <- fia %>%
        filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]])) %>%
        sample_n(min(100, n()))
      
      map <- create_base_map() %>%
        setView(lng = -73, lat = 44, zoom = 8)
      
      # Add points
      map <- map %>%
        addCircleMarkers(
          data = fia_sample,
          lng = ~get(lon_col),
          lat = ~get(lat_col),
          radius = 5,
          color = "#e74c3c",
          fillOpacity = 0.8,
          stroke = TRUE,
          weight = 1,
          layerId = ~CN,
          popup = ~paste0(
            "<strong>Plot:</strong> ", CN, "<br>",
            "<strong>Biomass:</strong> ", round(aglb_Mg_per_ha, 1), " Mg/ha<br>",
            "<em>Click to show uncertainty radius</em>"
          )
        )
      
      # Add uncertainty circles
      map <- map %>%
        addCircles(
          data = fia_sample,
          lng = ~get(lon_col),
          lat = ~get(lat_col),
          radius = 1609,  # 1 mile in meters
          stroke = TRUE,
          color = "#e74c3c",
          weight = 1,
          dashArray = "5,5",
          fillOpacity = 0.05,
          group = "Uncertainty"
        )
      
      map
    })
    
    # Temporal plot
    output$temporal_plot <- renderPlotly({
      data <- app_data()
      
      plot_data <- data.frame(year = integer(), count = integer(), source = character())
      
      if (!is.null(data$fia) && "MEASYEAR" %in% names(data$fia)) {
        fia_years <- data$fia %>%
          filter(!is.na(MEASYEAR)) %>%
          count(MEASYEAR) %>%
          rename(year = MEASYEAR) %>%
          mutate(source = "FIA")
        plot_data <- bind_rows(plot_data, fia_years)
      }
      
      if (!is.null(data$nefin) && "MEASYEAR" %in% names(data$nefin)) {
        nefin_years <- data$nefin %>%
          filter(!is.na(MEASYEAR)) %>%
          count(MEASYEAR) %>%
          rename(year = MEASYEAR) %>%
          mutate(source = "NEFIN")
        plot_data <- bind_rows(plot_data, nefin_years)
      }
      
      if (nrow(plot_data) == 0) {
        return(plotly_empty() %>% layout(title = "No temporal data available"))
      }
      
      p <- ggplot(plot_data, aes(x = year, y = n, fill = source)) +
        geom_area(alpha = 0.7, position = "identity") +
        scale_fill_manual(values = c("FIA" = "#e74c3c", "NEFIN" = "#27ae60")) +
        labs(
          title = "Measurement Years",
          x = "Year",
          y = "Number of Plots",
          fill = "Dataset"
        ) +
        theme_dashboard()
      
      ggplotly(p) %>%
        layout(legend = list(orientation = "h", y = -0.1))
    })
    
    # Year summary table
    output$year_table <- renderTable({
      data <- app_data()
      
      result <- data.frame(
        Metric = c("Year Range", "Peak Year", "Total Plots")
      )
      
      if (!is.null(data$fia) && "MEASYEAR" %in% names(data$fia)) {
        fia_years <- data$fia$MEASYEAR[!is.na(data$fia$MEASYEAR)]
        peak_fia <- as.numeric(names(sort(table(fia_years), decreasing = TRUE)[1]))
        result$FIA <- c(
          paste(min(fia_years), "-", max(fia_years)),
          as.character(peak_fia),
          format(length(fia_years), big.mark = ",")
        )
      }
      
      if (!is.null(data$nefin) && "MEASYEAR" %in% names(data$nefin)) {
        nefin_years <- data$nefin$MEASYEAR[!is.na(data$nefin$MEASYEAR)]
        if (length(nefin_years) > 0) {
          peak_nefin <- as.numeric(names(sort(table(nefin_years), decreasing = TRUE)[1]))
          result$NEFIN <- c(
            paste(min(nefin_years), "-", max(nefin_years)),
            as.character(peak_nefin),
            format(length(nefin_years), big.mark = ",")
          )
        }
      }
      
      result
    }, striped = TRUE, hover = TRUE, width = "100%")
    
    # Biomass violin plot
    output$biomass_violin <- renderPlotly({
      data <- app_data()
      
      plot_data <- data.frame(
        biomass = numeric(),
        source = character()
      )
      
      if (!is.null(data$fia) && "aglb_Mg_per_ha" %in% names(data$fia)) {
        fia_bio <- data.frame(
          biomass = data$fia$aglb_Mg_per_ha,
          source = "FIA"
        ) %>% filter(!is.na(biomass))
        plot_data <- bind_rows(plot_data, fia_bio)
      }
      
      if (!is.null(data$nefin) && "aglb_Mg_per_ha" %in% names(data$nefin)) {
        nefin_bio <- data.frame(
          biomass = data$nefin$aglb_Mg_per_ha,
          source = "NEFIN"
        ) %>% filter(!is.na(biomass))
        plot_data <- bind_rows(plot_data, nefin_bio)
      }
      
      if (nrow(plot_data) == 0) {
        return(plotly_empty() %>% layout(title = "No biomass data"))
      }
      
      p <- ggplot(plot_data, aes(x = source, y = biomass, fill = source)) +
        geom_violin(alpha = 0.7) +
        geom_boxplot(width = 0.2, alpha = 0.8) +
        scale_fill_manual(values = c("FIA" = "#e74c3c", "NEFIN" = "#27ae60")) +
        labs(
          title = "Biomass Distribution",
          x = "",
          y = "Biomass (Mg/ha)"
        ) +
        theme_dashboard() +
        theme(legend.position = "none")
      
      ggplotly(p)
    })
    
    # Biomass density plot
    output$biomass_density <- renderPlotly({
      data <- app_data()
      
      plot_data <- data.frame(
        biomass = numeric(),
        source = character()
      )
      
      if (!is.null(data$fia) && "aglb_Mg_per_ha" %in% names(data$fia)) {
        fia_bio <- data.frame(
          biomass = data$fia$aglb_Mg_per_ha,
          source = "FIA"
        ) %>% filter(!is.na(biomass), biomass > 0)
        plot_data <- bind_rows(plot_data, fia_bio)
      }
      
      if (!is.null(data$nefin) && "aglb_Mg_per_ha" %in% names(data$nefin)) {
        nefin_bio <- data.frame(
          biomass = data$nefin$aglb_Mg_per_ha,
          source = "NEFIN"
        ) %>% filter(!is.na(biomass), biomass > 0)
        plot_data <- bind_rows(plot_data, nefin_bio)
      }
      
      if (nrow(plot_data) == 0) {
        return(plotly_empty() %>% layout(title = "No biomass data"))
      }
      
      p <- ggplot(plot_data, aes(x = biomass, fill = source)) +
        geom_density(alpha = 0.5) +
        scale_fill_manual(values = c("FIA" = "#e74c3c", "NEFIN" = "#27ae60")) +
        labs(
          title = "Biomass Density Comparison",
          x = "Biomass (Mg/ha)",
          y = "Density",
          fill = "Dataset"
        ) +
        theme_dashboard()
      
      ggplotly(p) %>%
        layout(legend = list(orientation = "h", y = -0.15))
    })
    
    # Covariate map
    output$covariate_map <- renderLeaflet({
      data <- app_data()
      
      map <- create_base_map() %>%
        setView(lng = -73, lat = 44, zoom = 6)
      
      # Add plot layers based on selection
      if ("fia" %in% input$plot_layers && !is.null(data$fia)) {
        fia <- data$fia
        lon_col <- names(fia)[grepl("lon", names(fia), ignore.case = TRUE)][1]
        lat_col <- names(fia)[grepl("lat", names(fia), ignore.case = TRUE)][1]
        color_col <- input$color_by
        
        if (!is.null(lon_col) && !is.null(lat_col) && color_col %in% names(fia)) {
          fia_sample <- fia %>%
            filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]]),
                   !is.na(.data[[color_col]])) %>%
            sample_n(min(1500, n()))
          
          pal <- colorNumeric("YlGn", domain = fia_sample[[color_col]])
          
          map <- map %>%
            addCircleMarkers(
              data = fia_sample,
              lng = ~get(lon_col),
              lat = ~get(lat_col),
              radius = 4,
              color = ~pal(get(color_col)),
              fillOpacity = input$opacity,
              stroke = FALSE,
              group = "FIA",
              popup = ~paste0(
                "<strong>FIA Plot:</strong> ", CN, "<br>",
                "<strong>", color_col, ":</strong> ", round(get(color_col), 2)
              )
            )
        }
      }
      
      if ("nefin" %in% input$plot_layers && !is.null(data$nefin)) {
        nefin <- data$nefin
        lon_col <- names(nefin)[grepl("lon", names(nefin), ignore.case = TRUE)][1]
        lat_col <- names(nefin)[grepl("lat", names(nefin), ignore.case = TRUE)][1]
        color_col <- input$color_by
        
        if (!is.null(lon_col) && !is.null(lat_col) && color_col %in% names(nefin)) {
          nefin_sample <- nefin %>%
            filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]]),
                   !is.na(.data[[color_col]])) %>%
            sample_n(min(1500, n()))
          
          pal <- colorNumeric("YlOrRd", domain = nefin_sample[[color_col]])
          
          map <- map %>%
            addCircleMarkers(
              data = nefin_sample,
              lng = ~get(lon_col),
              lat = ~get(lat_col),
              radius = 4,
              color = ~pal(get(color_col)),
              fillOpacity = input$opacity,
              stroke = FALSE,
              group = "NEFIN",
              popup = ~paste0(
                "<strong>NEFIN Plot:</strong> ", CN, "<br>",
                "<strong>", color_col, ":</strong> ", round(get(color_col), 2)
              )
            )
        }
      }
      
      map %>%
        addLayersControl(
          overlayGroups = c("FIA", "NEFIN"),
          options = layersControlOptions(collapsed = FALSE)
        )
    })
    
    # Summary table
    output$summary_table <- renderDT({
      data <- app_data()
      
      # Build summary
      summary_list <- list()
      
      if (!is.null(data$fia)) {
        fia <- data$fia
        summary_list$FIA <- data.frame(
          Metric = c("Total Plots", "Mean Biomass (Mg/ha)", "Median Biomass", 
                     "SD Biomass", "Year Range"),
          Value = c(
            format(nrow(fia), big.mark = ","),
            round(mean(fia$aglb_Mg_per_ha, na.rm = TRUE), 1),
            round(median(fia$aglb_Mg_per_ha, na.rm = TRUE), 1),
            round(sd(fia$aglb_Mg_per_ha, na.rm = TRUE), 1),
            if ("MEASYEAR" %in% names(fia)) 
              paste(range(fia$MEASYEAR, na.rm = TRUE), collapse = "-") 
            else "N/A"
          )
        )
      }
      
      if (!is.null(data$nefin)) {
        nefin <- data$nefin
        summary_list$NEFIN <- data.frame(
          Metric = c("Total Plots", "Mean Biomass (Mg/ha)", "Median Biomass",
                     "SD Biomass", "Year Range"),
          Value = c(
            format(nrow(nefin), big.mark = ","),
            round(mean(nefin$aglb_Mg_per_ha, na.rm = TRUE), 1),
            round(median(nefin$aglb_Mg_per_ha, na.rm = TRUE), 1),
            round(sd(nefin$aglb_Mg_per_ha, na.rm = TRUE), 1),
            if ("MEASYEAR" %in% names(nefin))
              paste(range(nefin$MEASYEAR, na.rm = TRUE), collapse = "-")
            else "N/A"
          )
        )
      }
      
      if (length(summary_list) == 0) {
        return(data.frame(Message = "No data loaded"))
      }
      
      # Combine
      result <- summary_list[[1]]
      names(result)[2] <- names(summary_list)[1]
      
      if (length(summary_list) > 1) {
        for (i in 2:length(summary_list)) {
          result[[names(summary_list)[i]]] <- summary_list[[i]]$Value
        }
      }
      
      datatable(
        result,
        options = list(
          dom = 't',
          pageLength = 10,
          ordering = FALSE
        ),
        rownames = FALSE
      )
    })
    
    # =========================================================================
    # Displacement Analysis (FIA fuzzed vs NEFIN true)
    # =========================================================================
    
    # Compute matched plots with displacement
    matched_displacement <- reactive({
      data <- app_data()
      
      if (is.null(data$fia) || is.null(data$nefin)) {
        # Return demo data
        set.seed(42)
        n <- 200
        
        true_lons <- runif(n, -73.5, -70.5)
        true_lats <- runif(n, 42.5, 45.5)
        
        angles <- runif(n, 0, 2 * pi)
        distances <- runif(n, 50, 1600)
        
        fuzz_lons <- true_lons + (distances * cos(angles)) / (111320 * cos(true_lats * pi / 180))
        fuzz_lats <- true_lats + (distances * sin(angles)) / 111320
        
        return(data.frame(
          CN = paste0("DEMO_", 1:n),
          nefin_lon = true_lons,
          nefin_lat = true_lats,
          fia_lon = fuzz_lons,
          fia_lat = fuzz_lats,
          displacement_m = distances,
          state = sample(c("VT", "NH", "ME", "NY", "MA"), n, replace = TRUE)
        ))
      }
      
      fia <- data$fia
      nefin <- data$nefin
      
      # Find coordinate columns
      fia_lon <- names(fia)[grepl("lon", names(fia), ignore.case = TRUE)][1]
      fia_lat <- names(fia)[grepl("lat", names(fia), ignore.case = TRUE)][1]
      nefin_lon <- names(nefin)[grepl("lon", names(nefin), ignore.case = TRUE)][1]
      nefin_lat <- names(nefin)[grepl("lat", names(nefin), ignore.case = TRUE)][1]
      
      if (is.null(fia_lon) || is.null(nefin_lon)) return(NULL)
      
      # Join on CN
      matched <- inner_join(
        fia %>% 
          select(CN, fia_lon = all_of(fia_lon), fia_lat = all_of(fia_lat),
                 any_of(c("STATECD", "STATENM"))),
        nefin %>% 
          select(CN, nefin_lon = all_of(nefin_lon), nefin_lat = all_of(nefin_lat)),
        by = "CN"
      ) %>%
        filter(!is.na(fia_lon), !is.na(nefin_lon), !is.na(fia_lat), !is.na(nefin_lat))
      
      if (nrow(matched) == 0) return(NULL)
      
      # Calculate displacement
      matched <- matched %>%
        mutate(
          displacement_m = sqrt(
            ((fia_lon - nefin_lon) * 111320 * cos(nefin_lat * pi / 180))^2 +
            ((fia_lat - nefin_lat) * 111320)^2
          ),
          state = if ("STATENM" %in% names(.)) STATENM else 
                  if ("STATECD" %in% names(.)) as.character(STATECD) else "Unknown"
        )
      
      matched
    })
    
    # Mean displacement text
    output$mean_displacement_text <- renderText({
      matched <- matched_displacement()
      if (is.null(matched) || nrow(matched) == 0) {
        return("N/A")
      }
      paste0(round(mean(matched$displacement_m, na.rm = TRUE), 0), " m")
    })
    
    # Max displacement text
    output$max_displacement_text <- renderText({
      matched <- matched_displacement()
      if (is.null(matched) || nrow(matched) == 0) {
        return("N/A")
      }
      paste0(round(max(matched$displacement_m, na.rm = TRUE), 0), " m")
    })
    
    # Displacement map
    output$displacement_map <- renderLeaflet({
      matched <- matched_displacement()
      
      map <- create_base_map() %>%
        setView(lng = -72.5, lat = 44, zoom = 7)
      
      if (is.null(matched) || nrow(matched) == 0) {
        return(map %>% 
          addControl(
            html = "<div class='alert alert-warning'>No matched plots found</div>",
            position = "topright"
          ))
      }
      
      # Sample for performance
      if (nrow(matched) > 500) {
        matched <- matched %>% sample_n(500)
      }
      
      # Color palette for displacement
      if (input$color_by_displacement) {
        pal <- colorNumeric("YlOrRd", domain = matched$displacement_m)
        fia_colors <- pal(matched$displacement_m)
      } else {
        fia_colors <- rep("#e74c3c", nrow(matched))
      }
      
      # Add displacement lines if enabled
      if (input$show_displacement_lines) {
        for (i in 1:nrow(matched)) {
          map <- map %>%
            addPolylines(
              lng = c(matched$nefin_lon[i], matched$fia_lon[i]),
              lat = c(matched$nefin_lat[i], matched$fia_lat[i]),
              weight = 1.5,
              color = "#9b59b6",
              opacity = 0.5,
              group = "Displacement"
            )
        }
      }
      
      # Add FIA points (fuzzed)
      map <- map %>%
        addCircleMarkers(
          data = matched,
          lng = ~fia_lon,
          lat = ~fia_lat,
          radius = 4,
          color = fia_colors,
          fillOpacity = 0.8,
          stroke = TRUE,
          weight = 1,
          group = "FIA (fuzzed)",
          popup = ~paste0(
            "<strong>FIA Plot (Fuzzed)</strong><br>",
            "CN: ", CN, "<br>",
            "Displacement: ", round(displacement_m, 0), " m"
          )
        ) %>%
        # Add NEFIN points (true)
        addCircleMarkers(
          data = matched,
          lng = ~nefin_lon,
          lat = ~nefin_lat,
          radius = 4,
          color = "#27ae60",
          fillOpacity = 0.8,
          stroke = TRUE,
          weight = 1,
          group = "NEFIN (true)",
          popup = ~paste0(
            "<strong>NEFIN Plot (True)</strong><br>",
            "CN: ", CN
          )
        ) %>%
        addLayersControl(
          overlayGroups = c("FIA (fuzzed)", "NEFIN (true)", "Displacement"),
          options = layersControlOptions(collapsed = FALSE)
        )
      
      # Add color legend if coloring by displacement
      if (input$color_by_displacement) {
        map <- map %>%
          addLegend(
            position = "bottomright",
            pal = pal,
            values = matched$displacement_m,
            title = "Displacement (m)",
            opacity = 1
          )
      }
      
      map
    })
    
    # Displacement histogram
    output$displacement_hist <- renderPlotly({
      matched <- matched_displacement()
      
      if (is.null(matched) || nrow(matched) == 0) {
        return(plotly_empty() %>% layout(title = "No displacement data"))
      }
      
      p <- ggplot(matched, aes(x = displacement_m)) +
        geom_histogram(bins = 40, fill = "#9b59b6", alpha = 0.7, color = "white") +
        geom_vline(aes(xintercept = mean(displacement_m)), 
                   color = "#e74c3c", linetype = "dashed", linewidth = 1) +
        geom_vline(xintercept = 1609, color = "#333", linetype = "dotted", linewidth = 1) +
        annotate("text", x = 1609, y = Inf, label = "Max fuzzing (1.6km)", 
                 hjust = -0.1, vjust = 2, size = 3) +
        labs(
          title = "Distribution of FIA-NEFIN Coordinate Displacement",
          x = "Displacement (meters)",
          y = "Number of Plots"
        ) +
        theme_dashboard()
      
      ggplotly(p)
    })
    
    # Displacement by state
    output$displacement_by_state <- renderPlotly({
      matched <- matched_displacement()
      
      if (is.null(matched) || nrow(matched) == 0 || !"state" %in% names(matched)) {
        return(plotly_empty() %>% layout(title = "No state data"))
      }
      
      state_summary <- matched %>%
        group_by(state) %>%
        summarize(
          mean_displacement = mean(displacement_m, na.rm = TRUE),
          sd_displacement = sd(displacement_m, na.rm = TRUE),
          n = n(),
          .groups = "drop"
        ) %>%
        filter(n >= 5) %>%
        arrange(desc(mean_displacement))
      
      p <- ggplot(state_summary, aes(x = reorder(state, mean_displacement), 
                                     y = mean_displacement, fill = n)) +
        geom_col(alpha = 0.8) +
        geom_errorbar(aes(ymin = mean_displacement - sd_displacement,
                          ymax = mean_displacement + sd_displacement),
                      width = 0.2) +
        coord_flip() +
        scale_fill_viridis_c(option = "D", name = "N plots") +
        labs(
          title = "Mean Displacement by State",
          x = "",
          y = "Mean Displacement (m)"
        ) +
        theme_dashboard()
      
      ggplotly(p)
    })
    
  })
}
