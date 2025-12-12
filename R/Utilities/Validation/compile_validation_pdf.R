#' Compile Validation Plots into PDF
#'
#' Takes a compound_plots object (output from rtx()) and compiles it into a PDF
#'
#' @param compound_plots List of compound plots from rtx() function
#' @param output_dir Directory to save PDF
#' @param pdf_name Name of PDF file
#' @param add_plot_tags Logical, whether to add plot tags to subtitles (default TRUE)
#'
#' @return Path to created PDF file
#' @export
compile_validation_pdf <- function(compound_plots, 
                                   output_dir, 
                                   pdf_name,
                                   add_plot_tags = TRUE,
                                   skip_if_disabled = TRUE) {
  
  if (is.null(compound_plots) || length(compound_plots) == 0) {
    stop("No plots provided. compound_plots must be a non-empty list.")
  }
  
  cat(sprintf("\n📄 Compiling %d compounds into PDF...\n", length(compound_plots)))
  
  pdf_path <- file.path(output_dir, pdf_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Ensure no lingering graphics devices before opening PDF
  while (dev.cur() > 1) {
    dev.off()
    cat("Closed lingering graphics device\n")
  }
  
  pdf(pdf_path, width = 8.5, height = 11, family = "Helvetica", onefile = TRUE)
  
  first_page <- TRUE
  
  for (compound_id in names(compound_plots)) {
    compound <- compound_plots[[compound_id]]
    plots <- compound$plots
    
    if (length(plots) == 0) next
    
    cat(sprintf("  Adding %s (%s): %d plots\n", compound$short_name, compound$id, length(plots)))
    
    # Sort plot labels
    plot_labels <- names(plots)
    plot_labels <- plot_labels[order(
      as.numeric(gsub("F([0-9]+)_S([0-9]+)", "\\1", plot_labels)),
      as.numeric(gsub("F([0-9]+)_S([0-9]+)", "\\2", plot_labels))
    )]
    
    # Modify plots for PDF: optionally add plot tags to subtitles
    pdf_plots <- lapply(plot_labels, function(label) {
      plot_info <- plots[[label]]
      p <- plot_info$plot
      
      if (add_plot_tags) {
        sample_id <- plot_info$sample_id
        standard_file <- plot_info$standard_file
        plot_tag <- plot_info$plot_tag
        rt_range <- plot_info$rt_range
        
        # Add plot tag to subtitle
        new_subtitle <- sprintf("Sample: %s  |  Standard: %s  |  RT = %.3f min  |  %s",
                               sample_id, standard_file, mean(rt_range), plot_tag)
        
        p <- p +
          ggplot2::labs(subtitle = new_subtitle) +
          ggplot2::theme(
            plot.subtitle = ggtext::element_markdown(hjust = 0.5, face = "italic", size = 6,
                                                    color = "black", lineheight = 1.2, margin = ggplot2::margin(0, 0, 3, 0))
          )
      }
      
      p
    })
    
    # Create pages: 6 plots per page (3 rows x 2 columns)
    plots_per_page <- 6
    n_pages <- ceiling(length(pdf_plots) / plots_per_page)
    
    for (page_num in 1:n_pages) {
      start_idx <- (page_num - 1) * plots_per_page + 1
      end_idx <- min(page_num * plots_per_page, length(pdf_plots))
      page_plots <- pdf_plots[start_idx:end_idx]
      
      # Create title
      if (n_pages == 1) {
        title_text <- sprintf("%s (%s)", compound$short_name, compound$id)
      } else {
        title_text <- sprintf("%s (%s) - page %d/%d", compound$short_name, compound$id, page_num, n_pages)
      }
      
      title_grob <- grid::textGrob(
        title_text,
        gp = grid::gpar(fontsize = 14, fontface = "bold"),
        x = 0.5, y = 0.95, just = "top"
      )
      
      # Create grid: 2 columns x 3 rows
      grid_plot <- gridExtra::arrangeGrob(
        grobs = page_plots,
        ncol = 2,
        top = title_grob
      )
      
      tryCatch({
        if (!first_page) {
          grid::grid.newpage()
        }
        first_page <- FALSE
        grid::grid.draw(grid_plot)
      }, error = function(e) {
        cat(sprintf("  Error drawing compound %s page %d: %s\n", compound$id, page_num, e$message))
      })
    }
  }
  
  # Ensure PDF device is properly closed
  pdf_closed <- tryCatch({
    dev.off()
    cat(sprintf("\n✓ PDF saved to: %s\n", pdf_path))
    TRUE
  }, error = function(e) {
    cat(sprintf("Error closing PDF: %s\n", e$message))
    # Force close any remaining devices
    while (dev.cur() > 1) {
      try(dev.off(), silent = TRUE)
    }
    FALSE
  })
  
  if (!pdf_closed) {
    warning("PDF may not have been created successfully")
    return(invisible(NULL))
  }
  
  # Copy PDF to OneDrive if configured
  if (exists("config", envir = .GlobalEnv)) {
    config <- get("config", envir = .GlobalEnv)
    if (!is.null(config$paths$validation_plot_directory_onedrive)) {
      onedrive_pdf_path <- file.path(config$paths$validation_plot_directory_onedrive, pdf_name)
      cat(sprintf("\n📤 Copying PDF to OneDrive...\n"))
      
      copy_success <- tryCatch({
        file.copy(pdf_path, onedrive_pdf_path, overwrite = TRUE)
        TRUE
      }, error = function(e) {
        FALSE
      })
      
      if (copy_success && file.exists(onedrive_pdf_path)) {
        cat(sprintf("✓ PDF backed up to OneDrive: %s\n", onedrive_pdf_path))
      } else {
        warning(sprintf("Failed to copy PDF to OneDrive: %s", onedrive_pdf_path))
      }
    }
  }
  
  invisible(pdf_path)
}
