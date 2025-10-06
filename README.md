# NEFIN-vs-FIA-Hexagon-Method-Impact-Comparison

FIA-NEFIN Multi-Stage Processing Pipeline

Efficient, cacheable pipeline for computing forest metrics from FIA data with Monte Carlo positional uncertainty.

Architecture

The pipeline is split into 5 discrete stages, allowing you to:



Run expensive operations (jittering) once and cache results

Compute multiple metrics quickly using the same jitter library

Resume from any stage without re-running previous work



Pipeline Stages

┌─────────────────────────────────────────────────────────────┐

│ STAGE 1: Hex Grid Filtering                                │

│ Filter national hex grid to your study area                │

│ Output: data/hex/hex\_grid.geojson                          │

│ Run: Once per study area                                   │

└─────────────────────────────────────────────────────────────┘

&nbsp;                           ↓

┌─────────────────────────────────────────────────────────────┐

│ STAGE 2: FIA Data Pull                                     │

│ Download and filter FIA database                           │

│ Output: data/interim/fia\_region/ OR data/interim/fia/     │

│ Run: Once per data update                                  │

└─────────────────────────────────────────────────────────────┘

&nbsp;                           ↓

┌─────────────────────────────────────────────────────────────┐

│ STAGE 3: Plot-to-Hex Assignment (NEW)                     │

│ Assign FIA plots to hexes, preserve original coordinates  │

│ Output: data/processed/plot\_hex\_assignments.csv           │

│ Run: Once per hex grid + FIA data combination             │

└─────────────────────────────────────────────────────────────┘

&nbsp;                           ↓

┌─────────────────────────────────────────────────────────────┐

│ STAGE 4: Monte Carlo Jitter Library (NEW)                 │

│ Pre-generate N jittered coordinate sets for all plots     │

│ Output: data/processed/mc\_jitter\_library/                 │

│ Run: Once (or when changing jitter parameters)            │

│ ⏱️  EXPENSIVE: ~10-30 min for 100 reps                     │

└─────────────────────────────────────────────────────────────┘

&nbsp;                           ↓

┌─────────────────────────────────────────────────────────────┐

│ STAGE 5: Metric Computation (NEW)                         │

│ Compute metrics using cached jitter library               │

│ Output: runs/{run\_id}/hex\_{metric}\_results.csv           │

│ Run: As many times as you want, for any metric            │

│ ⚡ FAST: ~30 sec - 2 min (no re-jittering!)               │

└─────────────────────────────────────────────────────────────┘

Quick Start

First-Time Setup

Run the full pipeline once:

bash# Install all stages and build cache

Rscript run\_pipeline.R --all

This will:



Initialize project structure

Filter hex grid to your states

Pull FIA data

Assign plots to hexes

Build jitter library (100 replicates, ~15-20 min)

Compute initial metric (AGLB)



Compute Additional Metrics (Fast!)

After setup, computing new metrics is fast because it reuses the jitter library:

bash# Edit configs/process.yml, change:

\# metric: "carbon"



\# Then run (takes ~1 min):

Rscript run\_pipeline.R --compute

You can compute as many metrics as you want without re-jittering!

Available Metrics



aglb - Aboveground live biomass (Mg/ha)

carbon - Carbon stock (Mg C/ha)

mortality - Annual mortality rate (% or Mg/ha/yr)

growth - Annual growth increment (Mg/ha/yr) requires previous visit data

regeneration - Sapling→tree transition ratio



Configuration

configs/process.yml

yamlproject\_dir: "."

hex\_path: "data/hex/hex\_grid.geojson"



\# Choose metric: aglb, carbon, mortality, growth, regeneration

metric: "aglb"



\# Metric-specific parameters

metric\_params:

&nbsp; carbon\_fraction: 0.5

&nbsp; dbh\_sapling: 2.5

&nbsp; dbh\_tree: 5.0



\# Processing years

years: \[2018, 2019, 2020]

level\_window: 3



\# Monte Carlo settings (used in Stage 4)

mc\_reps: 100              # More reps = better precision, longer Stage 4

jitter\_radius\_m: 1609.34  # ~1 mile



\# Jitter constraints

mask:

&nbsp; use\_hex\_union: true

&nbsp; use\_state\_constraint: true

&nbsp; state\_geo\_path: "data/boundaries/states\_5070.geojson"

&nbsp; state\_field: "STATEFP"

&nbsp; max\_reroll: 10

Running Individual Stages

Stage 3: Re-assign Plots (if hex grid changes)

bashRscript run\_pipeline.R --assign --overwrite

Stage 4: Rebuild Jitter Library (if parameters change)

bash# After editing mc\_reps, jitter\_radius\_m, or mask settings:

Rscript run\_pipeline.R --jitter --overwrite

Stage 5: Compute Metrics Only (default)

bash# Just run the computation (fastest):

Rscript run\_pipeline.R --compute



\# Or equivalently:

Rscript run\_pipeline.R

File Structure

project/

├── R/

│   ├── stage2\_assign\_plots.R       # Stage 3: Plot→Hex assignment

│   ├── stage3\_build\_jitter\_library.R # Stage 4: MC jitter generation

│   ├── stage4\_compute\_metrics.R    # Stage 5: Fast metric computation

│   ├── compute\_metrics.R           # Metric computation functions

│   ├── process\_to\_hex.R            # Legacy + shared utilities

│   ├── fiaDataPull.R              # FIA data acquisition

│   └── ...

├── data/

│   ├── interim/

│   │   └── fia\_region/            # FIA data (Stage 2)

│   └── processed/

│       ├── plot\_hex\_assignments.csv      # Stage 3 output

│       └── mc\_jitter\_library/            # Stage 4 output

│           ├── jitter\_library.parquet   # Jittered coordinates

│           └── manifest.yml             # Metadata

├── runs/

│   └── {run\_id}/                  # Stage 5 outputs

│       ├── hex\_aglb\_results.csv   # Metric results

│       └── viz/                   # Optional visualizations

├── configs/

│   ├── process.yml                # Main configuration

│   ├── fia\_pull.yml              # FIA download settings

│   └── hex\_filter.yml            # Hex grid filtering

└── run\_pipeline.R                 # Master pipeline controller

Performance

Stage timings (40k plots, 100 MC reps, 7 states):

StageTimeCacheableNotes1. Hex Filter~5 sec✓Run once per study area2. FIA Pull~10-30 min✓Run once per data update3. Plot Assignment~30 sec✓Run once per hex×FIA combo4. Jitter Library~15-20 min✓Expensive but cached5. Metric Computation~1-2 min✗Fast, run repeatedly

Key insight: Stage 4 is expensive but runs once. Stage 5 is fast and runs many times for different metrics.

Advanced Usage

Batch Compute Multiple Metrics

bash# Create a simple loop script

for METRIC in aglb carbon mortality; do

&nbsp; sed -i "s/^metric:.\*/metric: \\"$METRIC\\"/" configs/process.yml

&nbsp; Rscript run\_pipeline.R --compute

done

Different Jitter Libraries

You can maintain multiple jitter libraries with different parameters:

bash# Build high-precision library (200 reps)

sed -i 's/mc\_reps:.\*/mc\_reps: 200/' configs/process.yml

Rscript run\_pipeline.R --jitter --overwrite



\# This creates a new library; metrics computed will use it

Inspect Cached Data

R# Check plot assignments

assignments <- readr::read\_csv("data/processed/plot\_hex\_assignments.csv")

summary(assignments)



\# Check jitter library

library(arrow)

jitters <- read\_parquet("data/processed/mc\_jitter\_library/jitter\_library.parquet")

head(jitters)



\# See how many jittered versions per plot

table(jitters$replicate\_id)

Troubleshooting

"Jitter library not found"

Run Stage 4:

bashRscript run\_pipeline.R --jitter

"Plot assignments not found"

Run Stage 3:

bashRscript run\_pipeline.R --assign

Changed hex grid or FIA data

Rebuild Stages 3-4:

bashRscript run\_pipeline.R --assign --jitter --overwrite

Metrics seem wrong

Check that the jitter library matches your current hex grid and FIA data:

Rlibrary(yaml)

meta <- read\_yaml("data/processed/mc\_jitter\_library/manifest.yml")

print(meta)

Benefits of This Architecture



Separation of Concerns: Expensive jittering separated from fast metric computation

Reproducibility: Same jitter library used across all metrics

Efficiency: Compute 10 metrics in the time old pipeline took for 1

Flexibility: Change metrics without re-running expensive operations

Scalability: Add more metrics without performance penalty

Transparency: Original FIA coordinates preserved for validation



Migration from Old Pipeline

If you have existing runs from process\_to\_hex.R:

bash# One-time: Build the cache

Rscript run\_pipeline.R --assign --jitter



\# Then use new fast pipeline

Rscript run\_pipeline.R --compute

The new pipeline produces the same results but runs much faster for subsequent metrics.

