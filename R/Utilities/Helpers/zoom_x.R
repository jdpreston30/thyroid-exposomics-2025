#' Zoom X-axis on Validation Plot
#'
#' Adjusts the x-axis (RT) limits on a validation plot without removing standard data.
#'
#' @param plot_obj Plot object with 'plot' field
#' @param xl Lower x-axis limit (will be rounded to nearest 0.025)
#' @param xu Upper x-axis limit (will be rounded to nearest 0.025)
#'
#' @return Modified plot object with updated x-axis limits
#' @export
zoom_x <- function(plot_obj, xl, xu, subfolder = "revised") {
  # Validate input
  if (!is.list(plot_obj) || !"plot" %in% names(plot_obj)) {
    stop("Invalid plot object. Expected a plot object with 'plot' field.")
  }
  
  # Create a copy to avoid modifying original
  modified_plot <- plot_obj
  
  # Round x limits to nearest 0.025
  xl <- round(xl / 0.025) * 0.025
  xu <- round(xu / 0.025) * 0.025
  x_limits <- c(xl, xu)
  
  # Adjust x-axis limits on the plot
  modified_plot$plot <- modified_plot$plot +
    ggplot2::scale_x_continuous(
      limits = x_limits,
      expand = ggplot2::expansion(mult = c(0.05, 0.05), add = 0),
      breaks = function(limits) seq(ceiling(limits[1] * 20) / 20, floor(limits[2] * 20) / 20, by = 0.05),
      minor_breaks = function(limits) seq(ceiling(limits[1] * 40) / 40, floor(limits[2] * 40) / 40, by = 0.025)
    )
  
  # Update rt_range metadata
  modified_plot$rt_range <- x_limits
  modified_plot$adjusted_rt_range <- x_limits
  
  # Automatically write small PNG
  write_small(modified_plot, subfolder = subfolder)
  
  return(modified_plot)
}
