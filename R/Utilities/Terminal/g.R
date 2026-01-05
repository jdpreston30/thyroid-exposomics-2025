#' Apply Global Adjustments to Plot(s)
#'
#' Takes one or more plot names, applies global_VP_adjust, and updates them in global environment.
#' Convenience function for manual inspection and adjustment.
#'
#' @param ... Plot tags (unquoted), e.g., F1_S1_CP2486, F6_S1_CP3113
#'
#' @return Invisibly returns the adjusted plot object(s)
#' @export
#'
#' @examples
#' g(F1_S1_CP2486)  # Applies adjustments to single plot
#' g(F1_S1_CP2486, F6_S1_CP3113, F1_S1_CP2107)  # Multiple plots
g <- function(...) {
  # Capture all unquoted names
  plot_names <- as.list(substitute(list(...)))[-1]
  plot_tags <- sapply(plot_names, deparse)
  
  if (length(plot_tags) == 0) {
    stop("No plot names provided. Usage: g(plot_name) or g(plot1, plot2, ...)")
  }
  
  # Create a temporary list with all plots
  temp_list <- list()
  
  for (plot_tag in plot_tags) {
    # Get plot from global environment
    if (!exists(plot_tag, envir = .GlobalEnv)) {
      stop(sprintf("Plot '%s' not found in global environment. Load it first with l(%s)", 
                   plot_tag, plot_tag))
    }
    temp_list[[plot_tag]] <- get(plot_tag, envir = .GlobalEnv)
  }
  
  # Apply adjustments to all plots at once
  adjusted_list <- global_VP_adjust(
    validation_plots = temp_list,
    validation_curated = validation_check_files,
    y_title = "← Standard | Tumor →",
    remove_vertical_grid = TRUE,
    remove_horizontal_grid = FALSE
  )
  
  # Update all in global environment
  for (plot_tag in plot_tags) {
    assign(plot_tag, adjusted_list[[plot_tag]], envir = .GlobalEnv)
  }
  
  cat(sprintf("✓ Applied global adjustments to %d plot(s)\n", length(plot_tags)))
  
  invisible(adjusted_list)
}
