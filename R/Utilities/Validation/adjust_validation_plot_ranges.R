#' Adjust RT Ranges for Validation Plots
#'
#' Takes a list of validation plots and adjusts their x-axis RT ranges based on
#' curated RT adjustment values.
#'
#' @param validation_plots List of plot data objects with plot and rt_range elements
#' @param validation_curated Tibble with columns: plot (plot_tag), rt_range (adjustment value), order, quality
#'
#' @return List of adjusted plot data objects with updated x-axis ranges
#'
#' @export
adjust_validation_plot_ranges <- function(validation_plots, validation_curated) {
  
  cat(sprintf("\nAdjusting RT ranges for %d plots...\n", nrow(validation_curated)))
  
  validation_plots_adjusted <- list()
  
  for (i in 1:nrow(validation_curated)) {
    row <- validation_curated[i, ]
    plot_tag <- row$plot
    rt_adjustment <- row$rt_range
    
    # Check if plot exists
    if (!plot_tag %in% names(validation_plots)) {
      warning(sprintf("Plot %s not found in validation_plots - skipping", plot_tag))
      next
    }
    
    # Get original plot data
    plot_data <- validation_plots[[plot_tag]]
    original_plot <- plot_data$plot
    original_rt_range <- plot_data$rt_range
    
    # Calculate original RT center from rt_range
    original_rt_center <- mean(original_rt_range)
    
    # Calculate new RT range: center +/- rt_adjustment, rounded to nearest 0.025
    new_rt_min <- round((original_rt_center - rt_adjustment) / 0.025) * 0.025
    new_rt_max <- round((original_rt_center + rt_adjustment) / 0.025) * 0.025
    new_rt_range <- c(new_rt_min, new_rt_max)
    
    cat(sprintf("[%d/%d] %s: %.2f +/- %.4f -> [%.3f, %.3f]\n", 
                i, nrow(validation_curated), plot_tag, 
                original_rt_center, rt_adjustment, new_rt_min, new_rt_max))
    
    # Adjust x-axis limits
    adjusted_plot <- original_plot +
      scale_x_continuous(
        limits = new_rt_range,
        expand = expansion(mult = c(0.05, 0.05), add = 0),
        breaks = function(limits) seq(ceiling(limits[1] * 20) / 20, floor(limits[2] * 20) / 20, by = 0.05),
        minor_breaks = function(limits) seq(ceiling(limits[1] * 40) / 40, floor(limits[2] * 40) / 40, by = 0.025)
      )
    
    # Store adjusted plot with updated metadata
    validation_plots_adjusted[[plot_tag]] <- list(
      short_name = plot_data$short_name,
      id = plot_data$id,
      order = row$order,
      plot = adjusted_plot,
      sample_id = plot_data$sample_id,
      standard_file = plot_data$standard_file,
      plot_tag = plot_tag,
      rt_range = plot_data$rt_range,
      original_rt_range = original_rt_range,
      adjusted_rt_range = new_rt_range,
      quality = row$quality
    )
  }
  
  cat(sprintf("\nCompleted adjusting %d plots\n", length(validation_plots_adjusted)))
  
  return(validation_plots_adjusted)
}
