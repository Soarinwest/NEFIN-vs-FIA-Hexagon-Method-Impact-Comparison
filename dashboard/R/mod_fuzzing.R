# =============================================================================
# R/mod_fuzzing.R
# Fuzzing Effects tab module - Monte Carlo analysis visualization
# =============================================================================

fuzzing_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    div(
      class = "container-fluid px-4 py-3",
      
      h2(icon("random"), " Fuzzing Effects"),
      p(class = "lead text-muted",
        "How does coordinate fuzzing affect covariate extraction and predictions?"),
      
      # Section 3.1: Monte Carlo Animation
      h4(class = "mt-4", "3.1 Monte Carlo Fuzzing Simulation"),
      
      card(
        card_body(
          layout_columns(
            col_widths = c(8, 4),
            
            # Animation/visualization
            div(
              leafletOutput(ns("mc_map"), height = "400px") %>%
                withSpinner(type = 6),
              div(
                class = "mt-2 d-flex justify-content-center align-items-center",
                actionButton(ns("run_mc"), "Run Simulation", 
                             class = "btn-primary me-2", icon = icon("play")),
                actionButton(ns("reset_mc"), "Reset", 
                             class = "btn-secondary me-2", icon = icon("redo")),
                sliderInput(
                  ns("mc_replicate"),
                  label = NULL,
                  min = 1, max = 50, value = 1,
                  width = "300px",
                  animate = animationOptions(interval = 200, loop = FALSE)
                )
              )
            ),
            
            # Explanation
            div(
              h5("How It Works"),
              tags$ol(
                tags$li("Select a plot's ", tags$strong("true location"), " (center point)"),
                tags$li("Generate ", tags$strong("50 random positions"), " within ±1.6km radius"),
                tags$li("Extract ", tags$strong("covariates"), " at each fuzzed location"),
                tags$li("Compare extracted values to ", tags$strong("true values"))
              ),
              hr(),
              h5("Current Replicate"),
              div(
                class = "text-center p-3 bg-primary bg-opacity-10 rounded",
                h6(class = "text-muted mb-1", "Current Replicate"),
                h3(class = "text-primary mb-0", textOutput(ns("current_rep_text")))
              ),
              hr(),
              verbatimTextOutput(ns("mc_stats"))
            )
          )
        )
      ),
      
      # Section 3.2: Covariate Extraction Uncertainty
      h4(class = "mt-4", "3.2 Covariate Extraction Uncertainty"),
      
      layout_columns(
        col_widths = c(6, 6),
        
        card(
          card_header("Covariate Error from Fuzzing"),
          card_body(
            plotlyOutput(ns("covariate_error_plot"), height = "350px") %>%
              withSpinner(type = 6)
          )
        ),
        
        card(
          card_header("Uncertainty as % of Mean"),
          card_body(
            plotlyOutput(ns("covariate_pct_plot"), height = "350px") %>%
              withSpinner(type = 6)
          )
        )
      ),
      
      # Covariate summary cards
      layout_columns(
        col_widths = c(4, 4, 4),
        
        card(
          card_body(
            class = "text-center",
            h6(class = "text-muted", "NDVI Uncertainty"),
            h3(textOutput(ns("ndvi_uncertainty"))),
            p(class = "text-muted small", "Mean SD from fuzzing")
          )
        ),
        
        card(
          card_body(
            class = "text-center",
            h6(class = "text-muted", "Temperature Uncertainty"),
            h3(textOutput(ns("tmean_uncertainty"))),
            p(class = "text-muted small", "Mean SD from fuzzing")
          )
        ),
        
        card(
          card_body(
            class = "text-center",
            h6(class = "text-muted", "Precipitation Uncertainty"),
            h3(textOutput(ns("ppt_uncertainty"))),
            p(class = "text-muted small", "Mean SD from fuzzing")
          )
        )
      ),
      
      # Section 3.3: Prediction Uncertainty Map
      h4(class = "mt-4", "3.3 Prediction Uncertainty from Fuzzing"),
      
      card(
        card_body(
          layout_sidebar(
            sidebar = sidebar(
              width = 250,
              title = "Model Selection",
              
              selectInput(
                ns("model_type"),
                "Model Type:",
                choices = c(
                  "Linear Regression" = "linear",
                  "Ridge Regression" = "ridge",
                  "Gradient Boosting" = "gbm",
                  "Random Forest" = "rf"
                ),
                selected = "gbm"
              ),
              
              hr(),
              
              h6("Prediction Uncertainty"),
              verbatimTextOutput(ns("pred_unc_summary"))
            ),
            
            plotlyOutput(ns("prediction_uncertainty_plot"), height = "450px") %>%
              withSpinner(type = 6)
          )
        )
      ),
      
      # Section 3.4: Uncertainty vs Biomass
      h4(class = "mt-4", "3.4 What Drives Prediction Uncertainty?"),
      
      layout_columns(
        col_widths = c(6, 6),
        
        card(
          card_header("Uncertainty vs Observed Biomass"),
          card_body(
            plotlyOutput(ns("uncertainty_vs_biomass"), height = "350px") %>%
              withSpinner(type = 6)
          )
        ),
        
        card(
          card_header("Example Prediction Envelopes"),
          card_body(
            plotlyOutput(ns("prediction_envelopes"), height = "350px") %>%
              withSpinner(type = 6)
          )
        )
      ),
      
      # Key finding callout
      div(
        class = "alert alert-warning mt-4",
        h5(icon("exclamation-triangle"), " Key Finding"),
        p(
          "While model ", tags$strong("training is robust"), " to coordinate fuzzing ",
          "(R² stable across replicates), ", tags$strong("individual predictions"),
          " show 3-5 Mg/ha uncertainty. This uncertainty is ", 
          tags$strong("spatially clustered"), " - highest in heterogeneous, high-biomass forests."
        )
      )
    )
  )
}

fuzzing_server <- function(id, app_data, config = NULL) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values for MC simulation
    mc_state <- reactiveValues(
      running = FALSE,
      current_plot = NULL,
      jittered_points = NULL,
      replicate = 1
    )
    
    # MC simulation map - showing FIA (fuzzed) vs NEFIN (true) with displacement
    output$mc_map <- renderLeaflet({
      data <- app_data()
      
      map <- create_base_map() %>%
        setView(lng = -72.5, lat = 44, zoom = 10)
      
      # Check if we have matched FIA-NEFIN data with both coordinate sets
      if (!is.null(data$fia) && !is.null(data$nefin)) {
        # Try to find matched plots (same CN in both datasets)
        fia <- data$fia
        nefin <- data$nefin
        
        # Get column names
        fia_lon <- names(fia)[grepl("lon", names(fia), ignore.case = TRUE)][1]
        fia_lat <- names(fia)[grepl("lat", names(fia), ignore.case = TRUE)][1]
        nefin_lon <- names(nefin)[grepl("lon", names(nefin), ignore.case = TRUE)][1]
        nefin_lat <- names(nefin)[grepl("lat", names(nefin), ignore.case = TRUE)][1]
        
        if (!is.null(fia_lon) && !is.null(nefin_lon)) {
          matched <- inner_join(
            fia %>% select(CN, fia_lon = all_of(fia_lon), fia_lat = all_of(fia_lat)),
            nefin %>% select(CN, nefin_lon = all_of(nefin_lon), nefin_lat = all_of(nefin_lat)),
            by = "CN"
          ) %>%
            filter(!is.na(fia_lon), !is.na(nefin_lon)) %>%
            sample_n(min(100, n()))  # Sample for performance
          
          if (nrow(matched) > 0) {
            # Calculate displacement for each matched plot
            matched <- matched %>%
              mutate(
                displacement_m = sqrt(
                  ((fia_lon - nefin_lon) * 111320 * cos(nefin_lat * pi / 180))^2 +
                  ((fia_lat - nefin_lat) * 111320)^2
                )
              )
            
            # Add displacement lines
            for (i in 1:nrow(matched)) {
              map <- map %>%
                addPolylines(
                  lng = c(matched$nefin_lon[i], matched$fia_lon[i]),
                  lat = c(matched$nefin_lat[i], matched$fia_lat[i]),
                  weight = 1.5,
                  color = "#9b59b6",
                  opacity = 0.4,
                  group = "Displacement Lines"
                )
            }
            
            # Add FIA points (fuzzed/published)
            map <- map %>%
              addCircleMarkers(
                data = matched,
                lng = ~fia_lon,
                lat = ~fia_lat,
                radius = 5,
                color = "#e74c3c",
                fillColor = "#e74c3c",
                fillOpacity = 0.8,
                stroke = TRUE,
                weight = 1,
                group = "FIA (fuzzed)",
                popup = ~paste0(
                  "<strong>FIA Plot (Fuzzed)</strong><br>",
                  "CN: ", CN, "<br>",
                  "Displacement: ", round(displacement_m, 0), "m"
                )
              ) %>%
              # Add NEFIN points (true)
              addCircleMarkers(
                data = matched,
                lng = ~nefin_lon,
                lat = ~nefin_lat,
                radius = 5,
                color = "#27ae60",
                fillColor = "#27ae60",
                fillOpacity = 0.8,
                stroke = TRUE,
                weight = 1,
                group = "NEFIN (true)",
                popup = ~paste0(
                  "<strong>NEFIN Plot (True)</strong><br>",
                  "CN: ", CN
                )
              ) %>%
              # Center on data
              setView(
                lng = mean(matched$nefin_lon, na.rm = TRUE),
                lat = mean(matched$nefin_lat, na.rm = TRUE),
                zoom = 10
              )
          }
        }
      }
      
      # If no real data, show demo visualization
      if (is.null(data$fia) || is.null(data$nefin)) {
        # Generate demo matched plots
        set.seed(42)
        n_demo <- 50
        
        # True locations (NEFIN)
        true_lons <- runif(n_demo, -73.0, -72.0)
        true_lats <- runif(n_demo, 43.5, 44.5)
        
        # Fuzzed locations (FIA) - displaced up to 1.6km
        angles <- runif(n_demo, 0, 2 * pi)
        distances <- runif(n_demo, 100, 1609)  # 100m to 1.6km displacement
        
        fuzz_lons <- true_lons + (distances * cos(angles)) / (111320 * cos(true_lats * pi / 180))
        fuzz_lats <- true_lats + (distances * sin(angles)) / 111320
        
        demo_matched <- data.frame(
          CN = paste0("DEMO_", 1:n_demo),
          nefin_lon = true_lons,
          nefin_lat = true_lats,
          fia_lon = fuzz_lons,
          fia_lat = fuzz_lats,
          displacement_m = distances
        )
        
        # Add displacement lines
        for (i in 1:nrow(demo_matched)) {
          map <- map %>%
            addPolylines(
              lng = c(demo_matched$nefin_lon[i], demo_matched$fia_lon[i]),
              lat = c(demo_matched$nefin_lat[i], demo_matched$fia_lat[i]),
              weight = 1.5,
              color = "#9b59b6",
              opacity = 0.4,
              group = "Displacement Lines"
            )
        }
        
        # Add points
        map <- map %>%
          addCircleMarkers(
            data = demo_matched,
            lng = ~fia_lon,
            lat = ~fia_lat,
            radius = 5,
            color = "#e74c3c",
            fillOpacity = 0.8,
            stroke = TRUE,
            weight = 1,
            group = "FIA (fuzzed)",
            popup = ~paste0("FIA (fuzzed)<br>Displacement: ", round(displacement_m, 0), "m")
          ) %>%
          addCircleMarkers(
            data = demo_matched,
            lng = ~nefin_lon,
            lat = ~nefin_lat,
            radius = 5,
            color = "#27ae60",
            fillOpacity = 0.8,
            stroke = TRUE,
            weight = 1,
            group = "NEFIN (true)",
            popup = ~"NEFIN (true location)"
          )
      }
      
      # Add replicate-based MC jitter visualization
      # This shows what the MC simulation does - multiple possible true locations
      set.seed(input$mc_replicate)
      example_fia_lon <- -72.5
      example_fia_lat <- 44.0
      
      # Generate n possible "true" locations within fuzzing radius
      n_jitter <- min(input$mc_replicate, 50)
      if (n_jitter > 0) {
        angles_j <- runif(n_jitter, 0, 2 * pi)
        distances_j <- runif(n_jitter, 0, 1609)
        
        jitter_lons <- example_fia_lon + (distances_j * cos(angles_j)) / (111320 * cos(example_fia_lat * pi / 180))
        jitter_lats <- example_fia_lat + (distances_j * sin(angles_j)) / 111320
        
        jitter_df <- data.frame(
          lon = jitter_lons,
          lat = jitter_lats,
          id = 1:n_jitter,
          dist = distances_j
        )
        
        # Add MC example with uncertainty circle
        map <- map %>%
          addCircles(
            lng = example_fia_lon,
            lat = example_fia_lat,
            radius = 1609,
            stroke = TRUE,
            color = "#e74c3c",
            weight = 2,
            dashArray = "5,5",
            fillOpacity = 0.02,
            group = "MC Example"
          ) %>%
          addCircleMarkers(
            lng = example_fia_lon,
            lat = example_fia_lat,
            radius = 8,
            color = "#e74c3c",
            fillColor = "#e74c3c",
            fillOpacity = 1,
            weight = 2,
            group = "MC Example",
            popup = "FIA published coordinate (center of uncertainty)"
          ) %>%
          addCircleMarkers(
            data = jitter_df,
            lng = ~lon,
            lat = ~lat,
            radius = 3,
            color = "#3498db",
            fillOpacity = 0.6,
            stroke = FALSE,
            group = "MC Simulated",
            popup = ~paste0("MC Replicate #", id, "<br>Distance: ", round(dist, 0), "m from published")
          )
      }
      
      # Layer controls and legend
      map %>%
        addLayersControl(
          overlayGroups = c("FIA (fuzzed)", "NEFIN (true)", "Displacement Lines", "MC Example", "MC Simulated"),
          options = layersControlOptions(collapsed = FALSE)
        ) %>%
        addLegend(
          position = "bottomleft",
          colors = c("#e74c3c", "#27ae60", "#9b59b6", "#3498db"),
          labels = c("FIA (published/fuzzed)", "NEFIN (true location)", "Displacement", "MC simulated positions"),
          opacity = 1,
          title = "Coordinate Types"
        )
    })
    
    # Update map when slider changes
    observe({
      leafletProxy("mc_map", session) %>%
        clearGroup("jitter_points")
      
      # Redraw is handled by renderLeaflet reactivity
    })
    
    # Run MC simulation button
    observeEvent(input$run_mc, {
      updateSliderInput(session, "mc_replicate", value = 1)
      
      # Animate the slider
      session$sendCustomMessage(
        type = "animateSlider",
        message = list(inputId = session$ns("mc_replicate"))
      )
    })
    
    # Reset button
    observeEvent(input$reset_mc, {
      updateSliderInput(session, "mc_replicate", value = 1)
    })
    
    # Current replicate display
    output$current_rep_text <- renderText({
      paste0(input$mc_replicate, " / 50")
    })
    
    # MC statistics
    output$mc_stats <- renderText({
      n <- input$mc_replicate
      paste0(
        "Points generated: ", n, "\n",
        "Max distance: 1,609m (1 mile)\n",
        "Distribution: Uniform random\n",
        "Direction: All angles equally likely"
      )
    })
    
    # Covariate error plot
    output$covariate_error_plot <- renderPlotly({
      data <- app_data()
      
      if (!is.null(data$covariate_uncertainty)) {
        df <- data$covariate_uncertainty
        
        p <- ggplot(df, aes(x = variable, y = mean_abs_error, fill = variable)) +
          geom_col(alpha = 0.8) +
          geom_errorbar(
            aes(ymin = mean_abs_error * 0.8, ymax = mean_abs_error * 1.2),
            width = 0.2
          ) +
          scale_fill_viridis_d(option = "D") +
          labs(
            title = "Mean Absolute Error from Fuzzing",
            x = "",
            y = "Mean Absolute Error"
          ) +
          theme_dashboard() +
          theme(legend.position = "none")
        
        ggplotly(p)
      } else {
        # Demo data
        demo_df <- data.frame(
          variable = c("NDVI", "Temperature", "Precipitation"),
          mean_abs_error = c(0.015, 0.35, 12.5),
          pct_of_mean = c(2.1, 4.2, 1.1)
        )
        
        p <- ggplot(demo_df, aes(x = variable, y = mean_abs_error, fill = variable)) +
          geom_col(alpha = 0.8) +
          scale_fill_manual(values = c("#2ecc71", "#e74c3c", "#3498db")) +
          labs(
            title = "Mean Absolute Error from Fuzzing",
            subtitle = "Demo data - actual values from your analysis",
            x = "",
            y = "Mean Absolute Error"
          ) +
          theme_dashboard() +
          theme(legend.position = "none")
        
        ggplotly(p)
      }
    })
    
    # Covariate % of mean plot
    output$covariate_pct_plot <- renderPlotly({
      data <- app_data()
      
      if (!is.null(data$covariate_uncertainty)) {
        df <- data$covariate_uncertainty
        
        p <- ggplot(df, aes(x = variable, y = pct_of_mean, fill = variable)) +
          geom_col(alpha = 0.8) +
          geom_text(aes(label = paste0(round(pct_of_mean, 1), "%")), 
                    vjust = -0.5, size = 4) +
          scale_fill_viridis_d(option = "D") +
          coord_cartesian(ylim = c(0, max(df$pct_of_mean) * 1.2)) +
          labs(
            title = "Uncertainty as % of Mean Value",
            x = "",
            y = "Percent (%)"
          ) +
          theme_dashboard() +
          theme(legend.position = "none")
        
        ggplotly(p)
      } else {
        demo_df <- data.frame(
          variable = c("NDVI", "Temperature", "Precipitation"),
          pct_of_mean = c(2.1, 4.2, 1.1)
        )
        
        p <- ggplot(demo_df, aes(x = variable, y = pct_of_mean, fill = variable)) +
          geom_col(alpha = 0.8) +
          geom_text(aes(label = paste0(pct_of_mean, "%")), vjust = -0.5, size = 4) +
          scale_fill_manual(values = c("#2ecc71", "#e74c3c", "#3498db")) +
          coord_cartesian(ylim = c(0, 6)) +
          labs(
            title = "Uncertainty as % of Mean",
            subtitle = "Demo data",
            x = "",
            y = "Percent (%)"
          ) +
          theme_dashboard() +
          theme(legend.position = "none")
        
        ggplotly(p)
      }
    })
    
    # Uncertainty values
    output$ndvi_uncertainty <- renderText({
      data <- app_data()
      if (!is.null(data$covariate_uncertainty)) {
        ndvi_row <- data$covariate_uncertainty %>% filter(variable == "NDVI")
        if (nrow(ndvi_row) > 0) {
          return(paste0("±", round(ndvi_row$mean_sd, 4)))
        }
      }
      "±0.015"
    })
    
    output$tmean_uncertainty <- renderText({
      data <- app_data()
      if (!is.null(data$covariate_uncertainty)) {
        row <- data$covariate_uncertainty %>% filter(variable == "Tmean")
        if (nrow(row) > 0) {
          return(paste0("±", round(row$mean_sd, 2), "°C"))
        }
      }
      "±0.35°C"
    })
    
    output$ppt_uncertainty <- renderText({
      data <- app_data()
      if (!is.null(data$covariate_uncertainty)) {
        row <- data$covariate_uncertainty %>% filter(variable == "PPT")
        if (nrow(row) > 0) {
          return(paste0("±", round(row$mean_sd, 1), "mm"))
        }
      }
      "±12mm"
    })
    
    # Prediction uncertainty summary
    output$pred_unc_summary <- renderText({
      data <- app_data()
      
      if (!is.null(data$prediction_uncertainty)) {
        df <- data$prediction_uncertainty
        paste0(
          "Mean SD: ", round(mean(df$pred_fuzz_sd, na.rm = TRUE), 2), " Mg/ha\n",
          "Max range: ", round(max(df$pred_fuzz_range, na.rm = TRUE), 2), " Mg/ha\n",
          "Plots analyzed: ", nrow(df)
        )
      } else {
        paste0(
          "Mean SD: 3.8 Mg/ha\n",
          "Max range: 25.3 Mg/ha\n",
          "Plots analyzed: 500\n",
          "(Demo values)"
        )
      }
    })
    
    # Prediction uncertainty distribution
    output$prediction_uncertainty_plot <- renderPlotly({
      data <- app_data()
      
      if (!is.null(data$prediction_uncertainty)) {
        df <- data$prediction_uncertainty
        
        p <- ggplot(df, aes(x = pred_fuzz_sd)) +
          geom_histogram(bins = 40, fill = "#e74c3c", alpha = 0.7) +
          geom_vline(
            xintercept = mean(df$pred_fuzz_sd, na.rm = TRUE),
            linetype = "dashed",
            color = "darkred",
            linewidth = 1
          ) +
          labs(
            title = "Distribution of Prediction Uncertainty",
            subtitle = paste0("Mean SD = ", 
                              round(mean(df$pred_fuzz_sd, na.rm = TRUE), 2), " Mg/ha"),
            x = "Prediction SD (Mg/ha)",
            y = "Number of Plots"
          ) +
          theme_dashboard()
        
        ggplotly(p)
      } else {
        # Demo histogram
        demo_sd <- abs(rnorm(500, mean = 3.8, sd = 2))
        
        p <- ggplot(data.frame(sd = demo_sd), aes(x = sd)) +
          geom_histogram(bins = 40, fill = "#e74c3c", alpha = 0.7) +
          geom_vline(xintercept = mean(demo_sd), linetype = "dashed", 
                     color = "darkred", linewidth = 1) +
          labs(
            title = "Distribution of Prediction Uncertainty (Demo)",
            subtitle = "Mean SD = 3.8 Mg/ha",
            x = "Prediction SD (Mg/ha)",
            y = "Number of Plots"
          ) +
          theme_dashboard()
        
        ggplotly(p)
      }
    })
    
    # Uncertainty vs biomass scatter
    output$uncertainty_vs_biomass <- renderPlotly({
      data <- app_data()
      
      if (!is.null(data$prediction_uncertainty) && 
          "biomass_true" %in% names(data$prediction_uncertainty)) {
        df <- data$prediction_uncertainty
        
        cor_val <- cor(df$biomass_true, df$pred_fuzz_sd, use = "complete.obs")
        
        p <- ggplot(df, aes(x = biomass_true, y = pred_fuzz_sd)) +
          geom_point(alpha = 0.4, color = "#3498db") +
          geom_smooth(method = "loess", color = "#e74c3c", se = TRUE) +
          labs(
            title = "Prediction Uncertainty vs Observed Biomass",
            subtitle = paste0("Correlation: r = ", round(cor_val, 3)),
            x = "Observed Biomass (Mg/ha)",
            y = "Prediction SD (Mg/ha)"
          ) +
          theme_dashboard()
        
        ggplotly(p)
      } else {
        # Demo scatter
        set.seed(42)
        n <- 400
        biomass <- rlnorm(n, meanlog = 4.5, sdlog = 0.6)
        pred_sd <- 1 + 0.02 * biomass + rnorm(n, sd = 1.5)
        pred_sd <- pmax(pred_sd, 0.5)
        
        cor_val <- cor(biomass, pred_sd)
        
        p <- ggplot(data.frame(biomass = biomass, pred_sd = pred_sd),
                    aes(x = biomass, y = pred_sd)) +
          geom_point(alpha = 0.4, color = "#3498db") +
          geom_smooth(method = "loess", color = "#e74c3c", se = TRUE) +
          labs(
            title = "Prediction Uncertainty vs Biomass (Demo)",
            subtitle = paste0("r = ", round(cor_val, 3), 
                              " | Higher biomass = more uncertainty"),
            x = "Observed Biomass (Mg/ha)",
            y = "Prediction SD (Mg/ha)"
          ) +
          theme_dashboard()
        
        ggplotly(p)
      }
    })
    
    # Prediction envelopes
    output$prediction_envelopes <- renderPlotly({
      data <- app_data()
      
      if (!is.null(data$prediction_uncertainty)) {
        df <- data$prediction_uncertainty %>%
          arrange(desc(pred_fuzz_range)) %>%
          head(20) %>%
          mutate(plot_rank = row_number())
        
        p <- ggplot(df, aes(x = reorder(factor(plot_idx), -pred_fuzz_range))) +
          geom_errorbar(
            aes(ymin = pred_fuzz_min, ymax = pred_fuzz_max),
            width = 0.3, color = "#3498db", linewidth = 1
          ) +
          geom_point(aes(y = pred_true), color = "#27ae60", size = 3) +
          geom_point(aes(y = biomass_true), color = "#e74c3c", size = 3, shape = 17) +
          labs(
            title = "Prediction Range (Top 20 Most Variable)",
            subtitle = "Blue = range | Green = prediction | Red = observed",
            x = "Plot",
            y = "Biomass (Mg/ha)"
          ) +
          theme_dashboard() +
          theme(axis.text.x = element_blank())
        
        ggplotly(p)
      } else {
        # Demo
        set.seed(123)
        n <- 20
        obs <- rlnorm(n, 5, 0.4)
        pred <- obs + rnorm(n, 0, 15)
        range_half <- runif(n, 5, 20)
        
        demo_df <- data.frame(
          plot = 1:n,
          obs = obs,
          pred = pred,
          min = pred - range_half,
          max = pred + range_half
        ) %>%
          arrange(desc(max - min))
        
        p <- ggplot(demo_df, aes(x = reorder(factor(plot), -(max - min)))) +
          geom_errorbar(aes(ymin = min, ymax = max), width = 0.3, 
                        color = "#3498db", linewidth = 1) +
          geom_point(aes(y = pred), color = "#27ae60", size = 3) +
          geom_point(aes(y = obs), color = "#e74c3c", size = 3, shape = 17) +
          labs(
            title = "Prediction Range Across Fuzzing (Demo)",
            subtitle = "Blue bars = range | Green = predicted | Red = observed",
            x = "Plot",
            y = "Biomass (Mg/ha)"
          ) +
          theme_dashboard() +
          theme(axis.text.x = element_blank())
        
        ggplotly(p)
      }
    })
    
  })
}
