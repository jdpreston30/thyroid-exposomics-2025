#' Deconvolute Combined Fragment Plot into Individual Fragment Plots
#'
#' Takes a validation plot with multiple fragments overlaid and creates
#' individual plots for each fragment, arranged as a patchwork.
#'
#' @param plot_data List containing plot object with multiple fragment layers
#' @param layout Character: "vertical", "horizontal", or "grid" for patchwork arrangement
#'
#' @return Patchwork object with individual fragment plots
#'
#' @export
individual_fragments <- function(plot_data, layout = "vertical") {
  
  library(patchwork)
  library(ggplot2)
  library(dplyr)
  
  # Extract the original plot
  original_plot <- plot_data$plot
  
  # Get the plot data (extract from ggplot object)
  plot_df <- original_plot$data
  
  # Identify column with fragment info (mz_label in rtx plots)
  if (!"mz_label" %in% names(plot_df)) {
    stop("Cannot find 'mz_label' column in plot data. Available columns: ", 
         paste(names(plot_df), collapse = ", "))
  }
  
  # Get unique fragments
  fragments <- unique(plot_df$mz_label)
  cat(sprintf("Creating %d individual fragment plots\n", length(fragments)))
  
  # Get color mapping from original plot
  color_scale <- original_plot$scales$get_scales("colour")
  if (!is.null(color_scale)) {
    mz_colors <- color_scale$palette(length(fragments))
    names(mz_colors) <- fragments
  } else {
    mz_colors <- c("mz0" = "black", "mz1" = "red", "mz2" = "gold", "mz3" = "blue")[fragments]
  }
  
  # Store individual plots
  fragment_plots <- list()
  
  for (frag in fragments) {
    # Filter data for this fragment
    frag_data <- plot_df %>% filter(mz_label == frag)
    
    # Create individual plot
    p <- ggplot(frag_data, aes(x = rt, y = plot_intensity)) +
      geom_line(color = mz_colors[frag], linewidth = 0.4) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
      scale_x_continuous(
        expand = expansion(mult = c(0.05, 0.05), add = 0),
        breaks = function(limits) seq(ceiling(limits[1] * 20) / 20, floor(limits[2] * 20) / 20, by = 0.05),
        minor_breaks = function(limits) seq(ceiling(limits[1] * 40) / 40, floor(limits[2] * 40) / 40, by = 0.025)
      ) +
      labs(
        title = frag,
        subtitle = original_plot$labels$subtitle,
        x = "Retention Time (min)",
        y = original_plot$labels$y
      ) +
      # Copy theme from original plot
      original_plot$theme +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8)
      )
    
    fragment_plots[[frag]] <- p
  }
  
  # Arrange with patchwork
  if (layout == "vertical") {
    combined <- wrap_plots(fragment_plots, ncol = 1)
  } else if (layout == "horizontal") {
    combined <- wrap_plots(fragment_plots, nrow = 1)
  } else if (layout == "grid") {
    ncols <- ceiling(sqrt(length(fragments)))
    combined <- wrap_plots(fragment_plots, ncol = ncols)
  } else {
    stop("layout must be 'vertical', 'horizontal', or 'grid'")
  }
  
  return(combined)
}
