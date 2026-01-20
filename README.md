# FIA-NEFIN Hexagon Analysis Pipeline (Consolidated)

## Overview

This is the reorganized pipeline for comparing FIA (fuzzed coordinates) vs NEFIN (true coordinates) forest inventory data aggregated to hexagonal spatial units across multiple scales.

**Active Scripts:** 24 (down from 47)  
**Archived Scripts:** 22 (Phase 2, legacy, one-time fixes)

---

## Directory Structure

```
R/
├── run_pipeline.R              # Master pipeline controller
├── validate_dashboard_data.R   # Data validation checks
│
├── 00_utils/                   # Shared utility functions
│   ├── utils_spatial.R         # Coordinate transforms, hex assignment
│   ├── utils_metrics.R         # SE calculation, weighted means  
│   └── utils_scale_names.R     # Scale name standardization (fia→64kha)
│
├── 01_data_prep/               # Data acquisition & preparation
│   ├── 01_init_project.R       # Create directories
│   ├── 02_create_hex_grid.R    # Generate hex grids at multiple scales
│   ├── 03_filter_hex_grid.R    # Clip grids to study area
│   ├── 04_fia_pull.R           # Download FIA data from FIADB
│   ├── 05_process_nefin.R      # Standardize NEFIN to FIA format
│   ├── 06_assign_plots.R       # Assign FIA plots to hexes
│   ├── 07_assign_nefin.R       # Assign NEFIN plots to hexes
│   └── 08_extract_covariates.R # Extract PRISM climate + NDVI
│
├── 02_uncertainty/             # Monte Carlo positional uncertainty
│   ├── 01_build_jitter_library.R  # Generate N jittered coordinate sets
│   └── 02_compute_metrics.R       # Compute hex metrics with MC integration
│
├── 03_comparison/              # FIA vs FIA+NEFIN comparison
│   ├── 01_compare_fia_nefin.R  # Per-hex comparison (FIXED augmented_se bug)
│   ├── 02_consolidate_results.R # Merge all scale results
│   └── 03_process_all_scales.R  # Multi-scale processing loop
│
├── 04_analysis/                # Statistical analysis
│   ├── 01_error_analysis.R     # Error decomposition (sampling vs positional)
│   ├── 02_advanced_analysis.R  # Summary statistics
│   ├── 03_phase1_hypothesis_tests.R  # Wilcoxon tests, CIs, effect sizes
│   └── 04_nefin_dominance_analysis.R # Bias risk from NEFIN dominance
│
├── 05_visualization/           # Publication figures
│   ├── 01_nefin_impact_figures.R   # NEFIN impact publication figures
│   ├── 02_visualize_results.R      # Consolidated results dashboard
│   └── 03_spatial_visualizations.R # Maps with hex grids
│
└── _archive/                   # Archived scripts (not part of core pipeline)
    ├── phase2/                 # Spatial modeling scripts (future work)
    ├── legacy/                 # Superseded scripts
    └── one_time_fixes/         # One-time data fixes
```

---

## Quick Start

### Full Pipeline
```bash
cd NEFIN-vs-FIA-Hexagon-Method-Impact-Comparison
Rscript run_pipeline.R --all
```

### Run Specific Stages
```bash
# Data preparation only
Rscript run_pipeline.R --data

# Just rerun comparison after fixing something
Rscript run_pipeline.R --compare --analyze --visualize

# Generate figures only
Rscript run_pipeline.R --visualize
```

### Validate Data
```bash
Rscript run_pipeline.R --validate
```

---

## Pipeline Stages

### Stage 1: Data Preparation (`01_data_prep/`)
- **Time:** ~30-60 min (FIA download is slowest)
- **Outputs:** `data/processed/`, hex grid GeoJSONs

### Stage 2: Monte Carlo Uncertainty (`02_uncertainty/`)
- **Time:** ~2-4 hours (100 MC replicates × 8 scales)
- **Outputs:** `data/processed/mc_jitter_library/`

### Stage 3: Comparison (`03_comparison/`)
- **Time:** ~10-15 min
- **Outputs:** `runs/consolidated_*/fia_nefin_comparison_all_scales.csv`

### Stage 4: Analysis (`04_analysis/`)
- **Time:** ~5 min
- **Outputs:** `runs/consolidated_*/advanced_analysis/`, `phase1_statistics/`, `dominance_analysis/`

### Stage 5: Visualization (`05_visualization/`)
- **Time:** ~5 min
- **Outputs:** `runs/consolidated_*/nefin_impact_figures/`, `visualizations/`

---

## Key Bug Fixes in This Version

### `compare_fia_nefin.R` - Fixed augmented_se calculation

**Old (buggy):**
```r
# When FIA had ≥3 plots, NEFIN was ignored:
w_fia = min(1, n_fia / 3)  # = 1 if n_fia >= 3
augmented_se = w_fia * fia_se + (1-w_fia) * nefin_se  # = fia_se when w_fia=1
```

**New (correct):**
```r
# Proper inverse-variance weighting:
augmented_se = 1 / sqrt(1/fia_se^2 + 1/nefin_se^2)  # Always < either SE
```

---

## Configuration

Edit `process.yml` to configure:
- Study area bounds
- Hex grid scales (default: 100ha to 100kha)
- MC replicates (default: 100)
- FIA states to include
- Output directories

---

## Archived Scripts

Scripts in `_archive/` are not needed for the core Phase 1 analysis:

### `_archive/phase2/` - Spatial Modeling (Future Work)
- `fia_nefin_spatial_modeling.R` - RF/GAM biomass models
- `fuzzing_effect_analysis.R` - Fuzzing impact visualization
- `sensor_resolution_comparison.R` - MODIS vs Sentinel-2
- `create_prediction_maps.R` - Hex-level prediction maps

### `_archive/legacy/` - Superseded
- `master_process_all.R` - Old pipeline controller
- `03_fia_pull_states.R` - State-by-state FIA pull (merged into 04_fia_pull.R)
- `process_multi_period.R` - Temporal analysis

### `_archive/one_time_fixes/`
- `batch_rename_fia_to_64kha.R` - One-time scale rename
- `fix_fia_to_64kha_columns.R` - One-time column fix

---

## Dependencies

```r
# Core
library(dplyr)
library(tidyr)
library(readr)
library(sf)
library(ggplot2)

# FIA
library(DBI)
library(RSQLite)

# Spatial
library(terra)
library(h3)

# Optional (for some visualizations)
library(patchwork)
library(viridis)
library(scales)
```

---

## Contact

Soren Donisvitch  
UVM FEMC / Forest Inventory Analysis
