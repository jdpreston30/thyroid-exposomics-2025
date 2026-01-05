#' Convert Raw Files to mzML Format
#'
#' This function handles the computationally intensive workflow of converting
#' .raw files to mzML format. It checks if conversion has already been completed
#' and skips redundant processing to save time in pipeline runs.
#'
#' @param file_list Tibble with ID, study, type, and files columns
#' @param tumor_raw_dir Path to tumor raw data directory
#' @param cadaver_raw_dir Path to cadaver raw data directory
#'
#' @return Tibble with file inventory including paths to raw and mzML files
#'
#' @export
convert_raw_to_mzml <- function(file_list,
                                tumor_raw_dir,
                                cadaver_raw_dir) {
  
  # Create output directories
  tumor_mzml_dir <- file.path(tumor_raw_dir, "mzML_validation")
  cadaver_mzml_dir <- file.path(cadaver_raw_dir, "mzML_validation")
  dir.create(tumor_mzml_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(cadaver_mzml_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Build file inventory
  file_inventory <- file_list |>
    separate_rows(files, sep = ", ") |>
    mutate(
      raw_file = paste0(files, ".raw"),
      raw_dir = case_when(
        study == "DTC" ~ tumor_raw_dir,
        study == "cadaver" ~ cadaver_raw_dir,
        TRUE ~ NA_character_
      ),
      mzml_dir = case_when(
        study == "DTC" ~ tumor_mzml_dir,
        study == "cadaver" ~ cadaver_mzml_dir,
        TRUE ~ NA_character_
      ),
      raw_path = file.path(raw_dir, raw_file),
      mzml_file = gsub("\\.raw$", ".mzML", raw_file, ignore.case = TRUE),
      mzml_path = file.path(mzml_dir, mzml_file)
    )
  
  # Check if all conversions are already complete
  all_mzml_exist <- all(file.exists(file_inventory$mzml_path))
  
  if (all_mzml_exist) {
    cat("\n✓ All mzML files already exist.\n")
    cat("  Skipping conversion to save time.\n")
    cat(sprintf("  Found %d mzML files.\n", 
                sum(file.exists(file_inventory$mzml_path))))
    return(file_inventory)
  }
  
  # Validate all raw files exist
  missing_files <- file_inventory |>
    filter(!file.exists(raw_path)) |>
    select(ID, study, type, raw_file, raw_path)
  
  if (nrow(missing_files) > 0) {
    cat("WARNING: The following files are missing:\n")
    print(missing_files)
    stop(sprintf("%d raw files not found. Please check file paths.", nrow(missing_files)))
  } else {
    cat(sprintf("✓ All %d raw files found and validated\n", nrow(file_inventory)))
  }
  
  # Setup conversion tools
  mono_path <- "/Library/Frameworks/Mono.framework/Versions/Current/Commands/mono"
  trfp_exe <- path.expand("~/bin/ThermoRawFileParser/ThermoRawFileParser.exe")
  
  if (!file.exists(mono_path)) {
    stop(sprintf("Mono not found at: %s\nPlease ensure Mono is installed.", mono_path))
  }
  if (!file.exists(trfp_exe)) {
    stop(sprintf("ThermoRawFileParser not found at: %s\nPlease ensure it's installed correctly.", trfp_exe))
  }
  
  # Convert files that don't exist yet
  files_to_convert <- file_inventory |>
    filter(!file.exists(mzml_path))
  
  if (nrow(files_to_convert) > 0) {
    cat(sprintf("\n=== Converting %d files to mzML ===\n", nrow(files_to_convert)))
    
    for (i in 1:nrow(files_to_convert)) {
      row <- files_to_convert[i, ]
      cat(sprintf("\n[%d/%d] Converting %s...\n", i, nrow(files_to_convert), row$raw_file))
      
      conversion_cmd <- sprintf('"%s" "%s" -i="%s" -o="%s" -f=1', 
                                mono_path, trfp_exe, row$raw_path, row$mzml_dir)
      
      result <- system(conversion_cmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE)
      
      if (file.exists(row$mzml_path)) {
        cat(sprintf("✓ Success: %s\n", row$mzml_file))
      } else {
        warning(sprintf("✗ Failed: %s (exit code: %d)\n", row$mzml_file, result))
      }
    }
    cat(sprintf("\n=== Conversion Complete ===\n"))
  } else {
    cat("\n✓ All files already converted, skipping conversion step\n")
  }
  
  # Return file inventory for downstream processing
  return(file_inventory)
}
