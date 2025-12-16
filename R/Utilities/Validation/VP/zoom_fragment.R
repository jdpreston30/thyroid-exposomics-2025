#' Zoom to or Remove Fragment in Validation Plot
#'
#' Filters validation plot to show only a specific mz fragment (positive number)
#' or to remove a specific fragment (negative number), adjusting Y-axis accordingly.
#'
#' @param plot_obj Plot object with 'plot' field
#' @param mz_fragment Numeric: positive to keep only that fragment, negative to remove it
#'                    (e.g., 1 keeps only mz1, -1 removes mz1)
#' @param subfolder Subfolder within Outputs/Validation/ to save to (default "revised")
#'
#' @return Modified plot object with filtered fragments
#' @export
zoom_fragment <- function(plot_obj, mz_fragment, subfolder = "revised", write_output = TRUE) {
  # Validate input
  if (!is.list(plot_obj) || !"plot" %in% names(plot_obj)) {
    stop("Invalid plot object. Expected a plot object with 'plot' field.")
  }
  
  if (!is.numeric(mz_fragment) || length(mz_fragment) != 1) {
    stop("mz_fragment must be a single numeric value (positive to keep, negative to remove)")
  }
  
  # Create a copy to avoid modifying original
  modified_plot <- plot_obj
  
  # Get original data
  plot_data <- modified_plot$plot$data
  
  # Check if mz_label exists
  if (!"mz_label" %in% names(plot_data)) {
    stop("Plot data does not contain 'mz_label' column")
  }
  
  available_fragments <- unique(plot_data$mz_label)
  
  # Determine target fragment label
  mz_num <- abs(mz_fragment)
  # Pattern handles both regular and bold markdown format: "mz3: ..." or "**mz3: ..."
  target_pattern <- paste0("(^|\\*\\*)mz", mz_num, "[:\\s]")
  
  # Find fragment that matches the pattern (e.g., "mz3: 221.2215" or "**mz3: 221.2215 **\***")
  matching_fragments <- grep(target_pattern, available_fragments, value = TRUE)
  
  if (length(matching_fragments) == 0) {
    stop(sprintf("Fragment 'mz%d' not found. Available: %s", 
                 mz_num, paste(available_fragments, collapse = ", ")))
  }
  
  # Use the first match (there should only be one)
  target_label <- matching_fragments[1]
  
  # Filter data based on sign
  if (mz_fragment > 0) {
    # Positive: keep only specified fragment
    plot_data_filtered <- plot_data[plot_data$mz_label == target_label, ]
    subtitle_text <- sprintf("Fragment: %s", target_label)
  } else {
    # Negative: remove specified fragment
    plot_data_filtered <- plot_data[plot_data$mz_label != target_label, ]
    subtitle_text <- NULL  # Don't modify subtitle for removals
  }
  
  # Update plot data
  modified_plot$plot$data <- plot_data_filtered
  
  # Adjust y-axis based on filtered data range
  y_max <- max(abs(plot_data_filtered$plot_intensity), na.rm = TRUE)
  
  modified_plot$plot <- modified_plot$plot +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.1)),
      limits = c(-y_max * 1.1, y_max * 1.1),
      n.breaks = 8,
      labels = function(x) scales::scientific(abs(x))
    )
  
  # Update plot subtitle to indicate fragment modification (only for positive fragments)
  if (!is.null(subtitle_text)) {
    modified_plot$plot <- modified_plot$plot +
      ggplot2::labs(
        subtitle = sprintf("Sample: %s  |  RT = %.2f min  |  %s", 
                           modified_plot$sample_id, 
                           mean(modified_plot$rt_range),
                           subtitle_text)
      )
  }
  
  # Automatically write small PNG if requested
  if (write_output) {
    write_small(modified_plot, subfolder = subfolder)
  }
  
  return(modified_plot)
}
