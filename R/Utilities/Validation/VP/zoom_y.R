#' Zoom Y-axis on Validation Plot
#'
#' Adjusts the y-axis (intensity) limits on a validation plot without removing data.
#'
#' @param plot_obj Plot object with 'plot' field
#' @param yl Lower y-axis limit (default: 0)
#' @param yu Upper y-axis limit
#' @param subfolder Subfolder for write_small output (default: "revised")
#'
#' @return Modified plot object with updated y-axis limits
#' @export
zoom_y <- function(plot_obj, yl = 0, yu, subfolder = "revised", write_output = TRUE) {
  # Validate input
  if (!is.list(plot_obj) || !"plot" %in% names(plot_obj)) {
    stop("Invalid plot object. Expected a plot object with 'plot' field.")
  }
  
  # Create a copy to avoid modifying original
  modified_plot <- plot_obj
  
  # Set y limits
  y_limits <- c(yl, yu)
  
  # Adjust y-axis limits on the plot
  modified_plot$plot <- modified_plot$plot +
    ggplot2::scale_y_continuous(
      limits = y_limits,
      expand = ggplot2::expansion(mult = c(0, 0.1)),
      n.breaks = 8,
      labels = scales::scientific
    )
  
  # Automatically write small PNG if requested
  if (write_output) {
    write_small(modified_plot, subfolder = subfolder)
  }
  
  return(modified_plot)
}
