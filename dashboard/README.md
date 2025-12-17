# NEFIN vs FIA: Interactive Dashboard

## Does Coordinate Precision Matter for Forest Biomass Predictions?

An interactive, map-centric dashboard that tells the story of how FIA's ±1.6km coordinate fuzzing affects forest biomass predictions and whether NEFIN's true coordinates provide measurable improvement.

![Dashboard Preview](www/preview.png)

## Features

### 5 Interactive Tabs

1. **Overview** - Landing page with key findings and study area map
2. **The Data** - Explore FIA vs NEFIN datasets side-by-side
3. **Fuzzing Effects** - Monte Carlo analysis visualization
4. **Model Comparison** - Compare model performance (Linear, Ridge, GBM, RF)
5. **Maps** - Full exploration mode with layer controls

### Key Visualizations

- Interactive Leaflet maps with plot locations and hex grids
- Plotly charts with linked brushing
- Real-time Monte Carlo simulation animation
- Feature importance with uncertainty bars
- Observed vs predicted scatter plots
- Spatial residual maps

## Quick Start

### 1. Install Dependencies

```r
# Install required packages
install.packages(c(
  "shiny", "bslib", "leaflet", "leaflet.extras", "plotly",
  "sf", "terra", "dplyr", "tidyr", "readr", "ggplot2",
  "viridis", "scales", "DT", "htmltools", "shinycssloaders"
))
```

### 2. Configure Data Paths

Edit `app.R` and update the `config` list to point to your data files:

```r
config <- list(
  fia_complete = "path/to/your/fia_complete.csv",
  nefin_processed = "path/to/your/nefin_processed.csv",
  # ... etc
)
```

### 3. Run the Dashboard

```r
shiny::runApp("dashboard")
```

Or from the command line:

```bash
cd dashboard
Rscript -e "shiny::runApp('.')"
```

## Data Requirements

The dashboard expects these files (paths configurable in `app.R`):

### Required
- `fia_complete.csv` - FIA plot data with biomass and covariates
- `nefin_processed.csv` - NEFIN plot data with biomass and covariates

### Optional (for full functionality)
- `states_5070.geojson` - State boundaries
- Hex grid GeoJSONs at various scales
- Model comparison results from your analysis
- Fuzzing analysis results

### Expected Columns

**FIA data:**
- `CN` - Plot identifier
- `aglb_Mg_per_ha` - Aboveground live biomass (Mg/ha)
- `lon_original` / `lat_original` - Coordinates
- `MEASYEAR` - Measurement year
- `ndvi_modis`, `tmean`, `ppt` - Covariates

**NEFIN data:**
- `CN` - Plot identifier
- `aglb_Mg_per_ha` - Biomass
- `lon_public` / `lat_public` - Coordinates
- Covariates as above

## Demo Mode

If data files are not found, the dashboard automatically loads demo data for testing the interface. Demo data includes:

- 5,000 synthetic FIA plots
- 3,000 synthetic NEFIN plots
- Realistic biomass distributions
- Random coordinates in the northeastern US

## Project Structure

```
dashboard/
├── app.R                 # Main Shiny app
├── R/
│   ├── mod_overview.R    # Overview tab module
│   ├── mod_data.R        # Data tab module
│   ├── mod_fuzzing.R     # Fuzzing effects module
│   ├── mod_models.R      # Model comparison module
│   ├── mod_maps.R        # Full maps module
│   └── utils.R           # Helper functions
├── www/
│   ├── styles.css        # Custom CSS
│   ├── logo.svg          # Project logo
│   └── logo.png          # Logo (PNG version)
├── data/
│   └── dashboard/        # Dashboard-specific data (optional)
└── README.md             # This file
```

## Deployment Options

### Local Development
```r
shiny::runApp("dashboard", launch.browser = TRUE)
```

### shinyapps.io
```r
library(rsconnect)
deployApp("dashboard")
```

### Docker
```dockerfile
FROM rocker/shiny-verse:4.3
RUN install2.r --error \
    leaflet leaflet.extras plotly terra sf bslib thematic
COPY dashboard /srv/shiny-server/
```

### Posit Connect
Push via RStudio IDE or rsconnect package.

## Customization

### Theming
The dashboard uses `bslib` with the Flatly bootswatch theme. Modify in `app.R`:

```r
theme = bs_theme(
  version = 5,
  bootswatch = "flatly",  # Change theme here
  primary = "#27ae60",    # Customize colors
  ...
)
```

### Adding Layers
Add new raster or vector layers in `mod_maps.R` by extending the observer that handles layer toggles.

### Color Palettes
Customize palettes in `utils.R`:
- `biomass_palette()` - Sequential green
- `difference_palette()` - Diverging blue-red
- `uncertainty_palette()` - Sequential orange-red

## Performance Tips

1. **Large datasets**: Pre-filter or sample data before adding to maps
2. **Rasters**: Convert to Cloud-Optimized GeoTIFFs (COGs) and use tile server
3. **Hex grids**: Simplify geometry for web display
4. **Caching**: Use `bindCache()` for expensive computations

## Credits

- **Author**: Soren Donisvitch
- **Affiliation**: Forest Ecosystem Monitoring Cooperative (FEMC) / UVM
- **Data Sources**: 
  - USDA Forest Service Forest Inventory and Analysis (FIA)
  - Northeast Forest Inventory Network (NEFIN)
  - MODIS/Sentinel-2 NDVI
  - PRISM Climate Data

## License

This dashboard is part of the NEFIN vs FIA comparison project. See main project repository for license terms.

## Support

For issues or questions:
1. Check the main project documentation
2. Open an issue on the project repository
3. Contact the author

---

*Built with R Shiny + Leaflet + Plotly*
