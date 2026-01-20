#!/usr/bin/env Rscript
# =============================================================================
# diagnose_se_issue.R
# Check why augmented_se == fia_only_se (the fix isn't working)
# =============================================================================

library(readr)
library(dplyr)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║  Diagnosing SE Issue                                                     ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

# Find consolidated directory
runs_dir <- "runs"
consolidated_dirs <- list.dirs(runs_dir, recursive = FALSE)
consolidated_dirs <- consolidated_dirs[grepl("consolidated", consolidated_dirs)]
if (length(consolidated_dirs) == 0) {
  stop("No consolidated directory found in runs/")
}
latest_dir <- sort(consolidated_dirs, decreasing = TRUE)[1]
cat("Using:", latest_dir, "\n\n")

# Load comparison data
comp_file <- file.path(latest_dir, "fia_nefin_comparison_all_scales.csv")
if (!file.exists(comp_file)) {
  stop("Comparison file not found: ", comp_file)
}

df <- read_csv(comp_file, show_col_types = FALSE)
cat("Loaded", nrow(df), "rows\n\n")

# Check columns
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("COLUMN CHECK\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
se_cols <- names(df)[grepl("se|SE", names(df), ignore.case = TRUE)]
cat("SE-related columns:", paste(se_cols, collapse = ", "), "\n\n")

# Check if augmented_se == fia_only_se
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("SE COMPARISON\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

if ("augmented_se" %in% names(df) && "fia_only_se" %in% names(df)) {
  same <- sum(df$augmented_se == df$fia_only_se, na.rm = TRUE)
  diff <- sum(df$augmented_se != df$fia_only_se, na.rm = TRUE)
  cat("Rows where augmented_se == fia_only_se:", same, "\n")
  cat("Rows where augmented_se != fia_only_se:", diff, "\n\n")
  
  if (diff > 0) {
    cat("✓ FIX IS WORKING - some hexes have different SE values\n")
    cat("\nSample of hexes where SE was reduced:\n")
    df %>%
      filter(augmented_se < fia_only_se) %>%
      select(hex_id, grid_scale, fia_only_se, augmented_se, 
             matches("nefin_se|nefin_n")) %>%
      head(10) %>%
      print()
  } else {
    cat("✗ FIX NOT WORKING - all SE values are identical\n")
  }
}

# Check NEFIN SE values
cat("\n═══════════════════════════════════════════════════════════════════════════\n")
cat("NEFIN SE CHECK\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

nefin_se_col <- names(df)[grepl("nefin.*se", names(df), ignore.case = TRUE)][1]
nefin_n_col <- names(df)[grepl("nefin.*n|n_plots_nefin", names(df), ignore.case = TRUE)][1]

if (!is.na(nefin_se_col)) {
  cat("NEFIN SE column:", nefin_se_col, "\n")
  cat("  NA values:", sum(is.na(df[[nefin_se_col]])), "/", nrow(df), "\n")
  cat("  Zero values:", sum(df[[nefin_se_col]] == 0, na.rm = TRUE), "\n")
  cat("  Valid (>0) values:", sum(df[[nefin_se_col]] > 0, na.rm = TRUE), "\n")
  
  if (sum(!is.na(df[[nefin_se_col]])) > 0) {
    cat("  Range:", range(df[[nefin_se_col]], na.rm = TRUE), "\n")
  }
} else {
  cat("✗ No NEFIN SE column found!\n")
  cat("  This is the problem - NEFIN SE is not being computed/saved\n")
}

if (!is.na(nefin_n_col)) {
  cat("\nNEFIN plot count column:", nefin_n_col, "\n")
  cat("  Distribution:\n")
  print(table(df[[nefin_n_col]], useNA = "ifany"))
  
  # Key diagnostic: how many hexes have >= 2 NEFIN plots?
  n_with_2plus <- sum(df[[nefin_n_col]] >= 2, na.rm = TRUE)
  cat("\n  Hexes with >=2 NEFIN plots:", n_with_2plus, "/", nrow(df), "\n")
  cat("  (Need >=2 plots to compute SD/SE)\n")
}

# Check the script version
cat("\n═══════════════════════════════════════════════════════════════════════════\n")
cat("SCRIPT VERSION CHECK\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

script_paths <- c(
  "R/03_comparison/01_compare_fia_nefin.R",
  "R/compare_fia_nefin.R"
)

for (path in script_paths) {
  if (file.exists(path)) {
    cat("Found:", path, "\n")
    # Check for the fix
    content <- readLines(path)
    has_fix <- any(grepl("1/sqrt\\(1/fia_only_se", content))
    has_old <- any(grepl("w_fia.*min\\(1.*n_fia", content))
    
    if (has_fix) {
      cat("  ✓ Contains inverse-variance fix\n")
    } else {
      cat("  ✗ Does NOT contain inverse-variance fix\n")
    }
    
    if (has_old) {
      cat("  ⚠ Still contains old buggy weighting code\n")
    }
  }
}

# Summary and recommendation
cat("\n═══════════════════════════════════════════════════════════════════════════\n")
cat("DIAGNOSIS SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

if (!is.na(nefin_se_col) && sum(is.na(df[[nefin_se_col]])) == nrow(df)) {
  cat("
PROBLEM: nefin_se is NA for ALL hexes

This happens because:
1. Most hexes have only 1 NEFIN plot
2. sd() returns NA when n=1 (can't compute SD from single value)
3. Therefore nefin_se = sd()/sqrt(n) = NA/1 = NA
4. When nefin_se is NA, the inverse-variance formula can't be used
5. Falls back to: augmented_se = fia_only_se

SOLUTION OPTIONS:
A) Use a different SE estimator for single-plot hexes (e.g., pooled variance)
B) Only claim SE improvement for hexes with >=2 NEFIN plots
C) Use NEFIN data to improve the MEAN estimate only, not SE
")
} else if (!is.na(nefin_n_col) && sum(df[[nefin_n_col]] >= 2, na.rm = TRUE) == 0) {
  cat("
PROBLEM: No hexes have >=2 NEFIN plots

The inverse-variance SE combination requires valid SE estimates from both
FIA and NEFIN. With only 1 plot per hex, we can't estimate NEFIN's SE.
")
} else {
  cat("
The data structure looks okay. The issue may be elsewhere.
Please share the output above for further diagnosis.
")
}
