#' Copy Raw Validation Plot RDS Files to Local Repository
#'
#' Copies individual validation plot RDS files from config-based directories
#' to a local output folder for repository storage.
#'
#' @param validation_curated Tibble with columns: source, plot
#' @param config List: configuration object with paths$validation_plot_directory
#' @param output_dir Character string: destination directory for copied files
#' @param overwrite Logical: if FALSE, skip files that already exist in output_dir
#'
#' @return Invisible NULL
#'
#' @export
copy_raw_validation_plots <- function(validation_curated,
                                      config,
                                      output_dir = "Outputs/Validation/ggplot_objects_raw",
                                      overwrite = TRUE,
                                      skip_if_disabled = TRUE) {
  
  # Check if function should skip based on config
  if (skip_if_disabled && !is.null(config$analysis$copy_rds_val_plots) && 
      !config$analysis$copy_rds_val_plots) {
    cat("\n=== Skipping copy of raw validation plots (copy_rds_val_plots = FALSE) ===\n")
    return(invisible(NULL))
  }
  
  # Map source types to their RDS folders
  source_to_folder <- list(
    "VD" = "variant_rtx",
    "IARC" = "iarc_tumor_rtx",
    "IARC and VD" = "variant_rtx"
  )
  
  # Get base directory
  rds_base_dir <- config$paths$validation_plot_directory
  
  # Create output directory if needed
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Process each row in validation_curated
  cat(sprintf("\nCopying %d RDS files to local repository...\n", nrow(validation_curated)))
  
  files_copied <- 0
  
  for (i in 1:nrow(validation_curated)) {
    row <- validation_curated[i, ]
    source_type <- row$source
    plot_tag <- row$plot
    
    # Determine source folder
    folder_name <- source_to_folder[[source_type]]
    
    if (is.null(folder_name)) {
      warning(sprintf("Unknown source type '%s' for %s - skipping", source_type, plot_tag))
      next
    }
    
    # Construct source path
    source_file <- file.path(rds_base_dir, folder_name, paste0(plot_tag, ".rds"))
    
    # Check if source file exists
    if (!file.exists(source_file)) {
      warning(sprintf("Source file not found: %s - skipping", source_file))
      next
    }
    
    # Construct destination path
    dest_file <- file.path(output_dir, paste0(plot_tag, ".rds"))
    
    # Skip if file exists and overwrite is FALSE
    if (!overwrite && file.exists(dest_file)) {
      cat(sprintf("[%d/%d] Skipped (exists): %s\n", i, nrow(validation_curated), plot_tag))
      next
    }
    
    # Copy file
    file.copy(source_file, dest_file, overwrite = TRUE)
    files_copied <- files_copied + 1
    
    cat(sprintf("[%d/%d] Copied: %s\n", i, nrow(validation_curated), plot_tag))
  }
  
  cat(sprintf("\nCompleted copying %d RDS files to: %s\n", files_copied, output_dir))
  
  invisible(NULL)
}
