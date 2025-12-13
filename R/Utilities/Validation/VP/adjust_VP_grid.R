#' Remove Grid Lines from Validation Plots
#'
#' Removes thin light gray vertical (and optionally horizontal) grid lines 
#' from plot background.
#'
#' @param validation_plots Named list of plot objects
#' @param remove_vertical Logical, remove vertical grid lines (default: TRUE)
#' @param remove_horizontal Logical, remove horizontal grid lines (default: FALSE)
#'
#' @return Named list of plots with updated theme
#' @export
adjust_VP_grid <- function(validation_plots, 
                           remove_vertical = TRUE, 
                           remove_horizontal = FALSE) {
  
  if (!is.list(validation_plots) || length(validation_plots) == 0) {
    stop("validation_plots must be a non-empty list")
  }
  
  cat(sprintf("\n🔧 Removing grid lines from %d plots...\n", length(validation_plots)))
  
  adjusted_plots <- validation_plots
  
  for (plot_tag in names(validation_plots)) {
    plot_obj <- validation_plots[[plot_tag]]
    
    # Validate plot structure
    if (!is.list(plot_obj) || !"plot" %in% names(plot_obj)) {
      warning(sprintf("Skipping %s: invalid plot structure", plot_tag))
      next
    }
    
    # Build theme modifications
    theme_updates <- list()
    
    if (remove_vertical) {
      theme_updates$panel.grid.major.x <- ggplot2::element_blank()
      theme_updates$panel.grid.minor.x <- ggplot2::element_blank()
    }
    
    if (remove_horizontal) {
      theme_updates$panel.grid.major.y <- ggplot2::element_blank()
      theme_updates$panel.grid.minor.y <- ggplot2::element_blank()
    }
    
    # Apply theme updates
    adjusted_plots[[plot_tag]]$plot <- plot_obj$plot +
      do.call(ggplot2::theme, theme_updates)
  }
  
  grid_msg <- c()
  if (remove_vertical) grid_msg <- c(grid_msg, "vertical")
  if (remove_horizontal) grid_msg <- c(grid_msg, "horizontal")
  
  cat(sprintf("✓ Removed %s grid lines\n", paste(grid_msg, collapse = " and ")))
  
  return(adjusted_plots)
}
