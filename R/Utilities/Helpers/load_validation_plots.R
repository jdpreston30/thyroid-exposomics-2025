#' Load Compiled Validation Plots
#'
#' Loads validation plots from compiled RDS file on OneDrive. Checks if
#' validation_plots already exists in environment and prompts user before reloading.
#'
#' @param onedrive_base_path Path to OneDrive validation plot directory
#' @param force_reload Logical, if TRUE skips prompt and reloads regardless
#'
#' @return Named list of ggplot objects
#' @export
load_validation_plots <- function(onedrive_base_path, force_reload = FALSE) {
  
  compiled_rds_path <- file.path(onedrive_base_path, "curated", "validation_plots_compiled.rds")
  
  #- Check if already exists
  if (exists("validation_plots", envir = .GlobalEnv) && !force_reload) {
    cat("⚠️  validation_plots already exists in environment.\n")
    response <- readline(prompt = "Reload from RDS? (y/n): ")
    
    if (tolower(trimws(response)) != "y") {
      cat("⏭️  Skipping reload, using existing validation_plots\n")
      return(invisible(get("validation_plots", envir = .GlobalEnv)))
    }
  }
  
  #- Load from RDS
  if (file.exists(compiled_rds_path)) {
    validation_plots <- readRDS(compiled_rds_path)
    cat(sprintf("✓ Loaded %d validation plots from compiled RDS\n", length(validation_plots)))
    return(validation_plots)
  } else {
    stop("⚠️  Compiled validation plots RDS not found. Run section 8.4.3 to create it.")
  }
}
