# R/clip_existing_hex_grids.R
# Clip waterbodies from existing hex grids

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(fs)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

clip_water_from_existing <- function(hex_path,
                                     water_path,
                                     out_path = NULL,
                                     backup = TRUE,
                                     min_land_fraction = 0.5) {
  
  message("\n═══════════════════════════════════════════════════════════")
  message("CLIPPING WATERBODIES FROM EXISTING HEX GRID")
  message("═══════════════════════════════════════════════════════════\n")
  
  # Check inputs exist
  if (!fs::file_exists(hex_path)) stop("Hex grid not found: ", hex_path)
  if (!fs::file_exists(water_path)) stop("Waterbodies not found: ", water_path)
  
  # Auto-generate output path if not provided
  if (is.null(out_path)) {
    # Insert "_nowater" before file extension
    base <- fs::path_ext_remove(hex_path)
    ext <- fs::path_ext(hex_path)
    out_path <- paste0(base, "_nowater.", ext)
  }
  
  # Backup original if requested
  if (backup && !fs::file_exists(paste0(hex_path, ".backup"))) {
    message("→ Creating backup: ", hex_path, ".backup")
    fs::file_copy(hex_path, paste0(hex_path, ".backup"))
  }
  
  # Turn off spherical geometry
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  # Load hex grid
  message("→ Loading hex grid: ", hex_path)
  hexes <- sf::st_read(hex_path, quiet = TRUE)
  hexes <- sf::st_make_valid(hexes)
  message("  Initial hexes: ", nrow(hexes))
  
  original_crs <- sf::st_crs(hexes)
  
  # Load waterbodies
  message("→ Loading waterbodies: ", water_path)
  water <- sf::st_read(water_path, quiet = TRUE)
  water <- sf::st_make_valid(water)
  message("  Waterbody features: ", nrow(water))
  
  # Transform both to Albers (5070) for accurate area calculations
  crs_ref <- sf::st_crs(5070)
  hexes_5070 <- sf::st_transform(hexes, crs_ref)
  water_5070 <- sf::st_transform(water, crs_ref)
  
  # Calculate original total hex area
  original_area_m2 <- sum(as.numeric(sf::st_area(hexes_5070)))
  original_area_acres <- original_area_m2 / 4046.86
  
  message("\n→ Creating waterbodies mask...")
  water_union <- sf::st_union(water_5070) |> sf::st_make_valid()
  
  water_area_m2 <- as.numeric(sf::st_area(water_union))
  water_area_acres <- water_area_m2 / 4046.86
  
  message("  Total water area: ", format(round(water_area_acres), big.mark = ","), " acres")
  
  # Find hexes that intersect water
  message("\n→ Identifying hexes that intersect water...")
  intersects_water <- sf::st_intersects(hexes_5070, water_union, sparse = FALSE)[,1]
  n_intersect <- sum(intersects_water)
  
  message("  Hexes intersecting water: ", n_intersect, " (", 
          round(100 * n_intersect / nrow(hexes_5070), 1), "%)")
  
  if (n_intersect == 0) {
    message("\n✓ No hexes intersect water - no clipping needed!")
    message("  Original grid is already water-free.")
    return(invisible(hex_path))
  }
  
  # Process hexes
  message("\n→ Clipping water from hexes...")
  
  # Keep non-intersecting hexes as-is
  hexes_no_water <- hexes_5070[!intersects_water, ]
  
  # Clip water from intersecting hexes
  hexes_with_water <- hexes_5070[intersects_water, ]
  
  message("  Processing ", n_intersect, " hexes with water...")
  
  # Subtract water (this may create multipolygons or remove hexes entirely)
  hexes_clipped <- sf::st_difference(hexes_with_water, water_union) |>
    sf::st_make_valid()
  
  # Some hexes might be empty or nearly empty after clipping
  hexes_clipped$area_after <- as.numeric(sf::st_area(hexes_clipped))
  
  # Calculate original area for comparison
  hexes_clipped <- hexes_clipped |>
    dplyr::mutate(
      hex_id_original = hex_id,
      land_fraction = area_after / (6 * (2590 * sqrt(3/4)))  # Approximate, will recalculate properly
    )
  
  # Better approach: get original hex area
  hexes_with_water$area_before <- as.numeric(sf::st_area(hexes_with_water))
  hexes_clipped$area_before <- hexes_with_water$area_before[match(hexes_clipped$hex_id, hexes_with_water$hex_id)]
  hexes_clipped$land_fraction <- hexes_clipped$area_after / hexes_clipped$area_before
  
  # Filter out hexes with insufficient land
  message("  Filtering hexes with <", round(100 * min_land_fraction), "% land coverage...")
  
  hexes_keep <- hexes_clipped |>
    dplyr::filter(land_fraction >= min_land_fraction) |>
    dplyr::select(-area_after, -area_before, -land_fraction, -hex_id_original)
  
  n_removed <- nrow(hexes_clipped) - nrow(hexes_keep)
  
  message("    Kept: ", nrow(hexes_keep), " hexes")
  message("    Removed: ", n_removed, " hexes (mostly water)")
  
  # Combine with non-water hexes
  hexes_final <- rbind(
    hexes_no_water,
    hexes_keep
  )
  
  # Re-assign sequential hex IDs
  message("\n→ Re-assigning hex IDs...")
  old_ids <- hexes_final$hex_id
  hexes_final$hex_id <- sprintf("H%06d", seq_len(nrow(hexes_final)))
  
  # Calculate final statistics
  final_area_m2 <- sum(as.numeric(sf::st_area(hexes_final)))
  final_area_acres <- final_area_m2 / 4046.86
  
  water_removed_acres <- original_area_acres - final_area_acres
  
  # Transform back to original CRS
  message("→ Transforming back to original CRS...")
  hexes_out <- sf::st_transform(hexes_final, original_crs)
  
  # Write output
  message("→ Writing to: ", out_path)
  fs::dir_create(fs::path_dir(out_path), recurse = TRUE)
  sf::st_write(hexes_out, out_path, delete_dsn = TRUE, quiet = TRUE)
  
  # Summary
  message("\n═══════════════════════════════════════════════════════════")
  message("SUMMARY")
  message("═══════════════════════════════════════════════════════════")
  message("  Input:  ", hex_path)
  message("  Output: ", out_path)
  if (backup) {
    message("  Backup: ", hex_path, ".backup")
  }
  message("")
  message("  Original hexes:      ", nrow(hexes))
  message("  Final hexes:         ", nrow(hexes_out))
  message("  Hexes removed:       ", nrow(hexes) - nrow(hexes_out), 
          " (", round(100 * (nrow(hexes) - nrow(hexes_out)) / nrow(hexes), 1), "%)")
  message("")
  message("  Original area:       ", format(round(original_area_acres), big.mark = ","), " acres")
  message("  Water removed:       ", format(round(water_removed_acres), big.mark = ","), " acres")
  message("  Final area:          ", format(round(final_area_acres), big.mark = ","), " acres")
  message("  Water coverage:      ", round(100 * water_removed_acres / original_area_acres, 1), "%")
  message("═══════════════════════════════════════════════════════════\n")
  
  invisible(list(
    input = hex_path,
    output = out_path,
    hexes_removed = nrow(hexes) - nrow(hexes_out),
    water_acres = water_removed_acres,
    water_pct = 100 * water_removed_acres / original_area_acres
  ))
}

# Batch process multiple grids
clip_all_grids <- function(hex_pattern = "data/hex/hex_grid_*.geojson",
                           water_path = "data/hex/Waterbodies_FeaturesToJSON.geojson",
                           suffix = "_nowater",
                           in_place = FALSE,
                           backup = TRUE) {
  
  message("\n╔══════════════════════════════════════════════════════════╗")
  message("║  Batch Water Clipping for All Hex Grids                  ║")
  message("╚══════════════════════════════════════════════════════════╝\n")
  
  # Find all matching hex grids
  hex_files <- Sys.glob(hex_pattern)
  
  # Exclude files that already have the suffix
  hex_files <- hex_files[!grepl(paste0(suffix, "\\."), hex_files)]
  
  if (!length(hex_files)) {
    message("✗ No hex grids found matching: ", hex_pattern)
    return(invisible(NULL))
  }
  
  message("Found ", length(hex_files), " hex grids to process:")
  for (f in hex_files) {
    message("  - ", f)
  }
  
  message("\nWaterbodies: ", water_path)
  
  if (in_place) {
    message("\n⚠ IN-PLACE MODE: Will overwrite original files")
    if (backup) {
      message("  (Backups will be created as .backup)")
    }
  } else {
    message("\nNew files will be created with '", suffix, "' suffix")
  }
  
  message("")
  readline(prompt = "Press [Enter] to continue or [Ctrl+C] to cancel...")
  
  results <- list()
  
  for (i in seq_along(hex_files)) {
    hex_file <- hex_files[i]
    
    message("\n\n")
    message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    message("Processing ", i, "/", length(hex_files), ": ", basename(hex_file))
    message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    
    if (in_place) {
      # Create temp file, then replace original
      temp_file <- paste0(hex_file, ".temp")
      
      result <- tryCatch({
        clip_water_from_existing(
          hex_path = hex_file,
          water_path = water_path,
          out_path = temp_file,
          backup = backup
        )
        
        # Replace original with clipped version
        fs::file_move(temp_file, hex_file)
        message("\n✓ Updated: ", hex_file)
        
        list(file = hex_file, status = "success")
      }, error = function(e) {
        message("\n✗ Error processing ", hex_file, ": ", e$message)
        if (fs::file_exists(temp_file)) fs::file_delete(temp_file)
        list(file = hex_file, status = "error", message = e$message)
      })
    } else {
      # Create new file with suffix
      base <- fs::path_ext_remove(hex_file)
      ext <- fs::path_ext(hex_file)
      out_file <- paste0(base, suffix, ".", ext)
      
      result <- tryCatch({
        clip_water_from_existing(
          hex_path = hex_file,
          water_path = water_path,
          out_path = out_file,
          backup = FALSE  # No backup needed when creating new file
        )
        list(file = out_file, status = "success")
      }, error = function(e) {
        message("\n✗ Error processing ", hex_file, ": ", e$message)
        list(file = hex_file, status = "error", message = e$message)
      })
    }
    
    results[[i]] <- result
  }
  
  # Summary
  message("\n\n")
  message("╔══════════════════════════════════════════════════════════╗")
  message("║  Batch Processing Complete                               ║")
  message("╚══════════════════════════════════════════════════════════╝\n")
  
  n_success <- sum(sapply(results, function(x) x$status == "success"))
  n_error <- sum(sapply(results, function(x) x$status == "error"))
  
  message("Successfully processed: ", n_success, "/", length(hex_files))
  if (n_error > 0) {
    message("Errors: ", n_error)
    message("\nFailed files:")
    for (r in results) {
      if (r$status == "error") {
        message("  ✗ ", r$file, " - ", r$message)
      }
    }
  }
  
  if (n_success > 0) {
    message("\n✓ Water-clipped grids ready to use!")
  }
  
  invisible(results)
}

# CLI interface
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  
  get_arg <- function(flag, default = NULL) {
    hit <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
  }
  
  # Check for batch mode
  if ("--all" %in% args) {
    clip_all_grids(
      hex_pattern = get_arg("--pattern", "data/hex/hex_grid_*.geojson"),
      water_path = get_arg("--water", "data/hex/Waterbodies_FeaturesToJSON.geojson"),
      suffix = get_arg("--suffix", "_nowater"),
      in_place = "--in-place" %in% args,
      backup = !("--no-backup" %in% args)
    )
  } else {
    # Single file mode
    hex_path <- get_arg("--hex")
    
    if (is.null(hex_path)) {
      cat("\nUsage:\n")
      cat("  Single file:\n")
      cat("    Rscript R/clip_existing_hex_grids.R --hex=path/to/grid.geojson --water=path/to/water.geojson\n")
      cat("\n")
      cat("  Batch mode:\n")
      cat("    Rscript R/clip_existing_hex_grids.R --all\n")
      cat("    Rscript R/clip_existing_hex_grids.R --all --in-place\n")
      cat("\n")
      cat("Options:\n")
      cat("  --hex=PATH          Input hex grid\n")
      cat("  --water=PATH        Waterbodies (default: data/hex/Waterbodies_FeaturesToJSON.geojson)\n")
      cat("  --out=PATH          Output path (default: input_nowater.geojson)\n")
      cat("  --all               Process all grids matching pattern\n")
      cat("  --pattern=GLOB      Pattern to match (default: data/hex/hex_grid_*.geojson)\n")
      cat("  --suffix=TEXT       Suffix for new files (default: _nowater)\n")
      cat("  --in-place          Overwrite original files\n")
      cat("  --no-backup         Skip backup creation\n")
      cat("\n")
      quit(status = 1)
    }
    
    clip_water_from_existing(
      hex_path = hex_path,
      water_path = get_arg("--water", "data/hex/Waterbodies_FeaturesToJSON.geojson"),
      out_path = get_arg("--out", NULL),
      backup = !("--no-backup" %in% args)
    )
  }
}