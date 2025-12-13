#' Quick Load Validation Plot
#'
#' Loads a single validation plot RDS file and assigns to global environment
#' with the plot name. Convenience function for manual inspection.
#'
#' @param plot_name Plot tag (unquoted), e.g., F1_S1_CP2486
#'
#' @return Invisibly returns the plot object
#' @export
#'
#' @examples
#' l(F1_S1_CP2486)  # Loads and assigns to F1_S1_CP2486
l <- function(plot_name) {
  # Capture the unquoted name
  plot_tag <- deparse(substitute(plot_name))
  
  # Construct path to RDS file
  base_path <- config$paths$validation_plot_directory_onedrive
  rds_path <- file.path(base_path, "curated", "original", paste0(plot_tag, ".rds"))
  
  # Check if file exists
  if (!file.exists(rds_path)) {
    stop(sprintf("Plot file not found: %s", rds_path))
  }
  
  # Load the plot
  plot_obj <- readRDS(rds_path)
  
  # Assign to global environment with the plot tag name
  assign(plot_tag, plot_obj, envir = .GlobalEnv)
  
  cat(sprintf("✓ Loaded %s\n", plot_tag))
  
  invisible(plot_obj)
}
