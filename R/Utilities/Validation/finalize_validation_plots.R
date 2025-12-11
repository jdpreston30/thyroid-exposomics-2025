library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(ggtext)

#' Finalize Validation Plots
#'
#' Takes a curated validation tibble and reads individual RDS plot files from
#' specified directories, adjusting x-axis ranges based on manually curated RT values,
#' then compiles into PDF.
#'
#' @param validation_curated Tibble with columns: order, source, id, short_name, 
#'   tumor_rt_range, quality, plot (format: F#_S#_ID), state, rt_range
#' @param variant_rtx_folder Character string: folder name for variant RDS files (default: "variant_rtx")
#' @param iarc_tumor_rtx_folder Character string: folder name for IARC tumor RDS files (default: "iarc_tumor_rtx")
#' @param config List: configuration object with paths$validation_plot_directory
#' @param output_dir Character string: output directory for final PDF
#' @param pdf_name Character string: name of output PDF file (default: "finalized_validation.pdf")
#' @param plots_per_page Integer: number of plots per page in 2-column grid (default: 6 for 3x2 layout)
#'
#' @return Invisible list of adjusted plots
#'
#' @export
finalize_validation_plots <- function(validation_curated,
                                      variant_rtx_folder = "variant_rtx",
                                      iarc_tumor_rtx_folder = "iarc_tumor_rtx",
                                      config = config,
                                      output_dir,
                                      pdf_name = "finalized_validation.pdf",
                                      plots_per_page = 6) {
  
  # Check required parameter
  if (missing(output_dir)) {
    stop("output_dir is required and must be specified")
  }
  
  # Determine base RDS directory
  if (!is.null(config) && !is.null(config$paths$validation_plot_directory)) {
    base_rds_dir <- config$paths$validation_plot_directory
  } else {
    base_rds_dir <- "validation_plots"
    warning("config$paths$validation_plot_directory not provided - using local 'validation_plots' directory")
  }
  
  # Map source types to folder paths
  source_folder_map <- list(
    VD = file.path(base_rds_dir, variant_rtx_folder),
    IARC = file.path(base_rds_dir, iarc_tumor_rtx_folder)
  )
  
  cat(sprintf("\nRDS base directory: %s\n", base_rds_dir))
  cat(sprintf("  VD plots: %s\n", source_folder_map$VD))
  cat(sprintf("  IARC tumor plots: %s\n", source_folder_map$IARC))
  
  # Sort validation_curated by order
  validation_curated <- validation_curated |>
    arrange(order)
  
  cat(sprintf("\nFinalizing %d validation plots...\n", nrow(validation_curated)))
  
  # Initialize storage for adjusted plots
  adjusted_plots <- list()
  
  # Process each row
  for (row_idx in 1:nrow(validation_curated)) {
    row <- validation_curated[row_idx, ]
    
    id_val <- row$id
    plot_tag <- row$plot
    rt_range_adjustment <- row$rt_range
    source_type <- row$source
    
    cat(sprintf("[%d/%d] Processing %s (%s)...\n", 
                row_idx, nrow(validation_curated), row$short_name, plot_tag))
    
    # Determine which folder to use based on source
    rds_folder <- source_folder_map[[source_type]]
    
    if (is.null(rds_folder)) {
      warning(sprintf("Unknown source type '%s' for %s - skipping", source_type, id_val))
      next
    }
    
    # Check if folder exists
    if (!dir.exists(rds_folder)) {
      warning(sprintf("RDS folder does not exist: %s - skipping %s", rds_folder, id_val))
      next
    }
    
    # Construct RDS file path
    rds_file <- file.path(rds_folder, paste0(plot_tag, ".rds"))
    
    # Check if RDS file exists
    if (!file.exists(rds_file)) {
      warning(sprintf("RDS file not found: %s - skipping", rds_file))
      next
    }
    
    # Read RDS file
    plot_data <- tryCatch({
      readRDS(rds_file)
    }, error = function(e) {
      warning(sprintf("Error reading RDS file %s: %s - skipping", rds_file, e$message))
      return(NULL)
    })
    
    if (is.null(plot_data)) {
      next
    }
    
    # Extract plot and metadata from RDS
    original_plot <- plot_data$plot
    original_rt_range <- plot_data$rt_range
    
    # Calculate new RT range: mean(original) +/- rt_range_adjustment
    original_rt_mean <- mean(original_rt_range)
    new_rt_range <- c(original_rt_mean - rt_range_adjustment, 
                      original_rt_mean + rt_range_adjustment)
    
    cat(sprintf("  Adjusting RT range: [%.2f, %.2f] -> [%.2f, %.2f]\n",
                original_rt_range[1], original_rt_range[2],
                new_rt_range[1], new_rt_range[2]))
    
    # Adjust x-axis limits on the plot
    adjusted_plot <- original_plot +
      scale_x_continuous(
        limits = new_rt_range,
        expand = expansion(mult = c(0.05, 0.05), add = 0),
        breaks = function(limits) seq(ceiling(limits[1] * 20) / 20, floor(limits[2] * 20) / 20, by = 0.05),
        minor_breaks = function(limits) seq(ceiling(limits[1] * 40) / 40, floor(limits[2] * 40) / 40, by = 0.025)
      )
    
    # Store adjusted plot with metadata
    adjusted_plots[[plot_tag]] <- list(
      plot = adjusted_plot,
      id = id_val,
      short_name = row$short_name,
      order = row$order,
      quality = row$quality,
      original_rt_range = original_rt_range,
      adjusted_rt_range = new_rt_range
    )
  }
  
  if (length(adjusted_plots) == 0) {
    stop("No plots were successfully processed")
  }
  
  # Compile into PDF
  cat(sprintf("\nCompiling %d plots into PDF...\n", length(adjusted_plots)))
  
  pdf_path <- file.path(output_dir, pdf_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Ensure no lingering graphics devices before opening PDF
  while (dev.cur() > 1) {
    dev.off()
    cat("Closed lingering graphics device\n")
  }
  
  pdf(pdf_path, width = 8.5, height = 11, family = "Helvetica", onefile = TRUE)
  
  # Get plot keys in order
  plot_keys <- names(adjusted_plots)
  
  # Create pages
  n_pages <- ceiling(length(plot_keys) / plots_per_page)
  n_cols <- 2
  n_rows <- plots_per_page / n_cols
  
  for (page_num in 1:n_pages) {
    start_idx <- (page_num - 1) * plots_per_page + 1
    end_idx <- min(page_num * plots_per_page, length(plot_keys))
    page_plot_keys <- plot_keys[start_idx:end_idx]
    
    # Extract plots for this page
    page_plots <- lapply(page_plot_keys, function(key) {
      adjusted_plots[[key]]$plot
    })
    
    # Create page title
    title_text <- sprintf("Finalized Validation Plots - Page %d/%d", page_num, n_pages)
    title_grob <- grid::textGrob(
      title_text,
      gp = grid::gpar(fontsize = 14, fontface = "bold"),
      x = 0.5, y = 0.95, just = "top"
    )
    
    # Create grid
    grid_plot <- gridExtra::arrangeGrob(
      grobs = page_plots,
      ncol = n_cols,
      top = title_grob
    )
    
    tryCatch({
      if (page_num > 1) {
        grid::grid.newpage()
      }
      grid::grid.draw(grid_plot)
    }, error = function(e) {
      cat(sprintf("  Error drawing page %d: %s\n", page_num, e$message))
    })
  }
  
  # Ensure PDF device is properly closed
  pdf_closed <- tryCatch({
    dev.off()
    cat(sprintf("\nPDF saved to: %s\n", pdf_path))
    TRUE
  }, error = function(e) {
    cat(sprintf("Error closing PDF: %s\n", e$message))
    # Force close any remaining devices
    while (dev.cur() > 1) {
      try(dev.off(), silent = TRUE)
    }
    FALSE
  })
  
  invisible(adjusted_plots)
}
