remove_standard <- function(plot_obj, xl = NULL, xu = NULL, subfolder = "revised") {
  # Validate input
  if (!is.list(plot_obj) || !"plot" %in% names(plot_obj)) {
    stop("Invalid plot object. Expected a plot object with 'plot' field.")
  }
  
  # Create a copy to avoid modifying original
  modified_plot <- plot_obj
  
  # Round x limits to nearest 0.025 if provided
  if (!is.null(xl)) {
    xl <- round(xl / 0.025) * 0.025
  }
  if (!is.null(xu)) {
    xu <- round(xu / 0.025) * 0.025
  }
  
  # Determine x-axis limits
  if (!is.null(xl) && !is.null(xu)) {
    x_limits <- c(xl, xu)
    new_rt_range <- x_limits
  } else {
    x_limits <- NULL
    new_rt_range <- modified_plot$rt_range
  }
  
  # Modify the plot to remove standard and show only sample (positive Y)
  modified_plot$plot <- modified_plot$plot +
    ggplot2::aes(x = rt, y = plot_intensity, color = mz_label, group = mz_label) +
    ggplot2::scale_x_continuous(
      limits = x_limits,
      expand = ggplot2::expansion(mult = c(0.05, 0.05), add = 0),
      breaks = function(limits) seq(ceiling(limits[1] * 20) / 20, floor(limits[2] * 20) / 20, by = 0.05),
      minor_breaks = function(limits) seq(ceiling(limits[1] * 40) / 40, floor(limits[2] * 40) / 40, by = 0.025)
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.1)),
      limits = function(x) c(0, max(x, na.rm = TRUE) * 1.1),
      n.breaks = 8
    ) +
    ggplot2::labs(
      subtitle = sprintf("Sample: %s  |  RT = %.2f min", 
                         modified_plot$sample_id, 
                         mean(new_rt_range)),
      y = "Intensity (Sample)"
    )
  
  # Filter the plot data to only include positive intensities
  plot_data <- modified_plot$plot$data
  plot_data_filtered <- plot_data[plot_data$plot_intensity > 0, ]
  modified_plot$plot$data <- plot_data_filtered
  
  # Automatically save the plot
  write_small(modified_plot, subfolder = subfolder)
  
  return(modified_plot)
}
