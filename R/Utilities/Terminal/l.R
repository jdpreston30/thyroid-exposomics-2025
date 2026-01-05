#' Quick Load Validation Plot(s)
#'
#' Loads one or more validation plot RDS files and assigns to global environment
#' with their plot names. Convenience function for manual inspection.
#'
#' @param ... Plot tags (unquoted), e.g., F1_S1_CP2486 or F1_S1_CP2486, F6_S1_CP3113
#'
#' @return Invisibly returns the plot object(s)
#' @export
#'
#' @examples
#' l(F1_S1_CP2486)  # Loads single plot
#' l(F1_S1_CP2486, F6_S1_CP3113, F1_S1_CP2107)  # Loads multiple plots
l <- function(...) {
  # Capture all unquoted names
  plot_names <- as.list(substitute(list(...)))[-1]
  plot_tags <- sapply(plot_names, deparse)
  
  if (length(plot_tags) == 0) {
    stop("No plot names provided. Usage: l(plot_name) or l(plot1, plot2, ...)")
  }
  
  base_path <- config$paths$validation_plot_directory_onedrive
  loaded_plots <- list()
  
  for (plot_tag in plot_tags) {
    # Construct path to RDS file
    rds_path <- file.path(base_path, "curated", "original", paste0(plot_tag, ".rds"))
    
    # Check if file exists
    if (!file.exists(rds_path)) {
      warning(sprintf("Plot file not found: %s", rds_path))
      next
    }
    
    # Load the plot
    plot_obj <- readRDS(rds_path)
    
    # Assign to global environment with the plot tag name
    assign(plot_tag, plot_obj, envir = .GlobalEnv)
    loaded_plots[[plot_tag]] <- plot_obj
    
    cat(sprintf("✓ Loaded %s\n", plot_tag))
  }
  
  if (length(loaded_plots) == 1) {
    invisible(loaded_plots[[1]])
  } else {
    invisible(loaded_plots)
  }
}
