# =============================================================================
# R/mod_overview.R
# Overview tab module - Landing page with key findings
# =============================================================================

overview_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Hero section
    div(
      class = "bg-dark text-white py-5 mb-4",
      style = "margin-top: -1rem; margin-left: -1rem; margin-right: -1rem;",
      div(
        class = "container text-center",
        h1(
          class = "display-4 fw-bold",
          "Does ±1.6km Coordinate Fuzzing Degrade Forest Biomass Predictions?"
        ),
        p(
          class = "lead mt-3",
          "A Monte Carlo analysis of 14,656 FIA plots across the Northeastern US"
        ),
        hr(class = "my-4 border-light"),
        div(
          class = "row justify-content-center mt-4",
          div(
            class = "col-md-3 col-sm-6 mb-3",
            div(
              class = "bg-success bg-opacity-25 rounded p-3",
              h2(class = "display-5 fw-bold mb-0", textOutput(ns("stat_r2_degrad"))),
              p(class = "mb-0 text-light", "Model R² Degradation")
            )
          ),
          div(
            class = "col-md-3 col-sm-6 mb-3",
            div(
              class = "bg-warning bg-opacity-25 rounded p-3",
              h2(class = "display-5 fw-bold mb-0", textOutput(ns("stat_pred_unc"))),
              p(class = "mb-0 text-light", "Individual Prediction Uncertainty")
            )
          ),
          div(
            class = "col-md-3 col-sm-6 mb-3",
            div(
              class = "bg-info bg-opacity-25 rounded p-3",
              h2(class = "display-5 fw-bold mb-0", textOutput(ns("stat_mc_reps"))),
              p(class = "mb-0 text-light", "Monte Carlo Replicates")
            )
          )
        )
      )
    ),
    
    # Main content
    div(
      class = "container-fluid px-4",
      
      # Study area map
      layout_columns(
        col_widths = c(7, 5),
        card(
          card_header(
            class = "bg-primary text-white",
            icon("map-marker-alt"),
            " Study Area"
          ),
          card_body(
            leafletOutput(ns("study_area_map"), height = "450px") %>%
              withSpinner(type = 6, color = "#27ae60")
          ),
          card_footer(
            class = "text-muted small",
            textOutput(ns("map_caption"))
          )
        ),
        
        # Key findings sidebar
        card(
          card_header(
            class = "bg-secondary text-white",
            icon("lightbulb"),
            " Key Findings"
          ),
          card_body(
            h5(icon("check-circle", class = "text-success"), " Model Training is Robust"),
            p(
              class = "text-muted ms-4",
              "Cross-validation R² remains stable across fuzzing replicates. 
               The ensemble of trees learns to average out positional noise."
            ),
            hr(),
            
            h5(icon("exclamation-triangle", class = "text-warning"), " Individual Predictions Vary"),
            p(
              class = "text-muted ms-4",
              "Single-plot predictions show 3-5 Mg/ha uncertainty due to fuzzing.
               This matters for plot-level applications."
            ),
            hr(),
            
            h5(icon("map", class = "text-info"), " Uncertainty is Spatially Clustered"),
            p(
              class = "text-muted ms-4",
              "High-biomass forests and heterogeneous landscapes show the largest
               fuzzing-induced uncertainty."
            ),
            hr(),
            
            h5(icon("tree", class = "text-success"), " Biomass Drives Uncertainty"),
            p(
              class = "text-muted ms-4",
              "Correlation between biomass and prediction uncertainty: r = 0.26.
               Denser forests = more variation when fuzzing shifts covariates."
            ),
            hr(),
            
            h5(icon("robot", class = "text-primary"), " Tree Models More Sensitive"),
            p(
              class = "text-muted ms-4",
              "Random Forest and Gradient Boosting show higher MC variance than 
               linear/ridge regression, though all models remain stable overall."
            )
          )
        )
      ),
      
      # Methods overview
      h3(class = "mt-4 mb-3", icon("flask"), " Methodology"),
      
      layout_columns(
        col_widths = c(4, 4, 4),
        
        card(
          card_header(class = "bg-light", "1. Data Sources"),
          card_body(
            tags$ul(
              tags$li(tags$strong("FIA:"), " 14,656 plots with fuzzed coordinates (±1.6km)"),
              tags$li(tags$strong("NEFIN:"), " 15,306 plots with true coordinates"),
              tags$li(tags$strong("MODIS NDVI:"), " 250m 5-year composites"),
              tags$li(tags$strong("PRISM:"), " 4km temperature & precipitation")
            )
          )
        ),
        
        card(
          card_header(class = "bg-light", "2. Monte Carlo Simulation"),
          card_body(
            tags$ul(
              tags$li("Generate 50-100 jittered positions per plot"),
              tags$li("Extract covariates at each jittered location"),
              tags$li("Train/predict with each covariate set"),
              tags$li("Quantify prediction variance from fuzzing")
            )
          )
        ),
        
        card(
          card_header(class = "bg-light", "3. Model Comparison"),
          card_body(
            tags$ul(
              tags$li("Linear regression (baseline)"),
              tags$li("Ridge regression (regularized)"),
              tags$li("Gradient Boosting (XGBoost)"),
              tags$li("Random Forest"),
              tags$li("5-fold cross-validation + holdout test")
            )
          )
        )
      ),
      
      # Call to action
      div(
        class = "text-center mt-4 mb-4",
        p(class = "lead text-muted",
          "Explore the data, methods, and findings using the tabs above."),
        actionButton(
          ns("go_to_data"),
          "Start Exploring",
          class = "btn btn-primary btn-lg",
          icon = icon("arrow-right")
        )
      )
    )
  )
}

overview_server <- function(id, app_data, config = NULL) {
  moduleServer(id, function(input, output, session) {
    
    # Key statistics
    output$stat_r2_degrad <- renderText({
      "<0.2%"
    })
    
    output$stat_pred_unc <- renderText({
      "3-5 Mg/ha"
    })
    
    output$stat_mc_reps <- renderText({
      "50-100"
    })
    
    # Study area map
    output$study_area_map <- renderLeaflet({
      data <- app_data()
      
      map <- create_base_map()
      
      # Add state boundaries if available
      if (!is.null(data$states)) {
        map <- map %>%
          addPolygons(
            data = data$states,
            fillColor = "transparent",
            color = "#333",
            weight = 2,
            opacity = 1,
            label = ~STUSPS
          )
      }
      
      # Add FIA plot density as heatmap
      if (!is.null(data$fia)) {
        fia <- data$fia
        
        # Find coordinate columns
        lon_col <- names(fia)[grepl("lon", names(fia), ignore.case = TRUE)][1]
        lat_col <- names(fia)[grepl("lat", names(fia), ignore.case = TRUE)][1]
        
        if (!is.null(lon_col) && !is.null(lat_col)) {
          fia_coords <- fia %>%
            filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]]))
          
          map <- map %>%
            addHeatmap(
              data = fia_coords,
              lng = ~get(lon_col),
              lat = ~get(lat_col),
              intensity = 1,
              blur = 15,
              max = 0.5,
              radius = 8,
              gradient = c("0" = "blue", "0.5" = "yellow", "1" = "red")
            )
        }
      }
      
      # Add NEFIN plots as points
      if (!is.null(data$nefin)) {
        nefin <- data$nefin
        
        lon_col <- names(nefin)[grepl("lon", names(nefin), ignore.case = TRUE)][1]
        lat_col <- names(nefin)[grepl("lat", names(nefin), ignore.case = TRUE)][1]
        
        if (!is.null(lon_col) && !is.null(lat_col)) {
          nefin_coords <- nefin %>%
            filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]])) %>%
            sample_n(min(1000, n()))  # Limit points for performance
          
          map <- map %>%
            addCircleMarkers(
              data = nefin_coords,
              lng = ~get(lon_col),
              lat = ~get(lat_col),
              radius = 3,
              color = "#e67e22",
              fillOpacity = 0.5,
              stroke = FALSE,
              group = "NEFIN Plots"
            ) %>%
            addLayersControl(
              overlayGroups = "NEFIN Plots",
              options = layersControlOptions(collapsed = TRUE)
            )
        }
      }
      
      map
    })
    
    # Map caption
    output$map_caption <- renderText({
      data <- app_data()
      
      fia_count <- if (!is.null(data$fia)) nrow(data$fia) else 0
      nefin_count <- if (!is.null(data$nefin)) nrow(data$nefin) else 0
      
      paste0(
        "FIA plot density heatmap (", format(fia_count, big.mark = ","), " plots) | ",
        "Orange points: NEFIN plots (", format(nefin_count, big.mark = ","), " plots)"
      )
    })
    
    # Navigation to Data tab
    observeEvent(input$go_to_data, {
      updateNavbarPage(session$parentSession, "main_nav", selected = "The Data")
    })
    
  })
}
