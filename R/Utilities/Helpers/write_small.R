write_small <- function(aps, compound = NULL) {
  # Create output directory
  output_dir <- "Outputs/Validation/revised"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Case 1: aps is the full list and compound is specified (e.g., aps, "CP2382")
  if (!is.null(compound) && is.character(compound)) {
    if (!compound %in% names(aps)) {
      stop(sprintf("Compound '%s' not found in aps list", compound))
    }
    plots_to_save <- aps[[compound]]
  } 
  # Case 2: aps is a single plot object (legacy support)
  else if (is.list(aps) && "plot" %in% names(aps) && "short_name" %in% names(aps)) {
    plots_to_save <- list(aps)
    names(plots_to_save) <- aps$plot_tag
  }
  # Case 3: aps is a named list of plots (e.g., aps$CP2382)
  else if (is.list(aps) && length(aps) > 0) {
    # Check if first element is a plot object
    if (is.list(aps[[1]]) && "plot" %in% names(aps[[1]])) {
      plots_to_save <- aps
    } else {
      stop("Invalid input. Expected plot object(s) with 'plot', 'short_name', and 'id' fields.")
    }
  } else {
    stop("Invalid input. Provide either: 1) (aps, 'CP####'), 2) aps$CP####, or 3) single plot object")
  }
  
  # Save all plots
  saved_files <- character()
  for (plot_name in names(plots_to_save)) {
    plot_obj <- plots_to_save[[plot_name]]
    
    # Extract metadata
    plot <- plot_obj$plot
    short_name <- plot_obj$short_name
    id <- plot_obj$id
    
    # Create filename
    filename <- sprintf("%s (%s).png", short_name, id)
    filepath <- file.path(output_dir, filename)
    
    # Save plot
    ggsave(
      filename = filepath,
      plot = plot,
      width = 11 / 2,
      height = 8.5 / 2,
      dpi = 300,
      bg = "white"
    )
    
    cat(sprintf("✓ Saved: %s\n", filename))
    saved_files <- c(saved_files, filepath)
  }
  
  # Open the first saved file
  if (length(saved_files) > 0) {
    system2("open", shQuote(saved_files[1]))
  }
  
  invisible(saved_files)
}
