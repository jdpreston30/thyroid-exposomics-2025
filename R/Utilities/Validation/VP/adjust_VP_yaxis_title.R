#' Adjust Y-Axis Title in Validation Plots
#'
#' Universally changes the y-axis title across all validation plots.
#'
#' @param validation_plots Named list of plot objects
#' @param y_title Character string for new y-axis title (default: "Intensity")
#'
#' @return Named list of plots with updated y-axis titles
#' @export
adjust_VP_yaxis_title <- function(validation_plots, y_title = "Intensity") {
  
  if (!is.list(validation_plots) || length(validation_plots) == 0) {
    stop("validation_plots must be a non-empty list")
  }
  
  cat(sprintf("\n🔧 Adjusting y-axis title for %d plots...\n", length(validation_plots)))
  
  adjusted_plots <- validation_plots
  
  for (plot_tag in names(validation_plots)) {
    plot_obj <- validation_plots[[plot_tag]]
    
    # Validate plot structure
    if (!is.list(plot_obj) || !"plot" %in% names(plot_obj)) {
      warning(sprintf("Skipping %s: invalid plot structure", plot_tag))
      next
    }
    
    # Update y-axis title and optimize legend layout
    adjusted_plots[[plot_tag]]$plot <- plot_obj$plot +
      ggplot2::labs(y = y_title) +
      ggplot2::guides(color = ggplot2::guide_legend(
        nrow = 1,
        override.aes = list(size = 0.5),
        keywidth = ggplot2::unit(0.25, "cm"),
        keyheight = ggplot2::unit(0.3, "cm")
      )) +
      ggplot2::theme(
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.margin = ggplot2::margin(0, 0, 0, 0),
        legend.spacing.x = ggplot2::unit(0.1, "cm"),
        legend.key.size = ggplot2::unit(0.25, "cm")
      )
  }
  
  cat(sprintf("✓ Updated y-axis title to '%s' and reduced legend spacing\n", y_title))
  
  return(adjusted_plots)
}
