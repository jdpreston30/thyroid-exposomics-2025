#' Adjust Colors for Validation Plots
#'
#' Takes a list of validation plots and assigns consistent fragment-specific colors.
#' mz0 is always black, mz1 is always #E69F00, etc., even if some fragments are missing.
#'
#' @param validation_plots List of plot data objects with plot element
#'
#' @return List of plot data objects with updated color scales
#'
#' @export
adjust_VP_colors <- function(validation_plots) {
  
  #- Define color palette (9 colors for fragments mz0-mz8)
  cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  
  names(cbp2) <- paste0("mz", 0:8)
  
  cat(sprintf("\nAdjusting colors for %d plots...\n", length(validation_plots)))
  
  validation_plots_colored <- list()
  
  for (i in seq_along(validation_plots)) {
    plot_name <- names(validation_plots)[i]
    plot_data <- validation_plots[[i]]
    original_plot <- plot_data$plot
    
    #- Extract unique mz_labels from plot data
    if (!is.null(original_plot$data) && "mz_label" %in% names(original_plot$data)) {
      present_mz <- unique(original_plot$data$mz_label)
      present_mz <- present_mz[!is.na(present_mz)]
      
      #- Reorder data so mz with asterisk is drawn last (frontmost)
      plot_data_reordered <- original_plot$data
      mz_with_asterisk <- present_mz[grepl("\\*", present_mz)]
      
      if (length(mz_with_asterisk) > 0) {
        #- Sort data: non-asterisk first, then asterisk (so asterisk drawn on top)
        plot_data_reordered$draw_order <- ifelse(
          plot_data_reordered$mz_label %in% mz_with_asterisk, 
          2,  # Draw last (on top)
          1   # Draw first (behind)
        )
        plot_data_reordered <- plot_data_reordered[order(plot_data_reordered$draw_order), ]
        plot_data_reordered$draw_order <- NULL
        
        #- Update plot data
        original_plot$data <- plot_data_reordered
      }
      
      #- Get colors for present fragments
      colors_to_use <- cbp2[present_mz]
      
      #- Apply color scale
      adjusted_plot <- original_plot +
        scale_colour_manual(
          values = colors_to_use,
          breaks = present_mz,
          drop = FALSE
        )
      
      cat(sprintf("[%d/%d] %s: Applied colors for %s\n", 
                  i, length(validation_plots), plot_name, 
                  paste(present_mz, collapse = ", ")))
      
    } else {
      #- No mz_label found, skip color adjustment
      adjusted_plot <- original_plot
      cat(sprintf("[%d/%d] %s: No mz_label found - skipping\n", 
                  i, length(validation_plots), plot_name))
    }
    
    #- Store adjusted plot
    validation_plots_colored[[plot_name]] <- plot_data
    validation_plots_colored[[plot_name]]$plot <- adjusted_plot
  }
  
  cat(sprintf("\nCompleted color adjustment for %d plots\n\n", length(validation_plots_colored)))
  
  return(validation_plots_colored)
}
