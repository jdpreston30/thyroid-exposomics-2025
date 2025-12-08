# Required libraries: mzR, ggplot2, dplyr, gridExtra, ggtext
library(mzR)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(ggtext)

#' Batch Validation Chromatogram PDF Generator
#'
#' Iterates through a validation list tibble and creates mirrored chromatogram plots
#' for each sample/standard combination, compiling them into a single PDF with 
#' S1/S2 grid layout and labeled plot identifiers.
#'
#' @param validation_list Tibble with columns: order, id, short_name, monoisotopic, compound_rt_range, 
#'   mz0-mz3, standards (comma-separated), file1-file6, asterisk, f1_rt-f6_rt, s1_rt-s2_rt
#' @param study Character string: "tumor" or "cadaver" to determine directory (default: "tumor")
#' @param iterate_through Integer: how many of the top sample files to process (default: 5, dynamic based on sample_top_files list)
#' @param output_dir Required character string: output directory for final PDF
#' @param pdf_name Character string: name of output PDF file (default: "validation_spectra.pdf")
#' @param ppm_tolerance Numeric mass tolerance in ppm (default: 5)
#' @param rt_lookup Character: "range" uses compound_rt_range (default), "sample" uses file-specific RT ranges
#' @param stick Logical, whether to plot as vertical sticks (default: FALSE)
#' @param max_i Logical, whether to use maximum intensity (default: FALSE)
#' @param plot_width Plot width in inches (default: 3.9)
#' @param plot_height Plot height in inches (default: 3.25)
#'
#' @return Named list of all plots (invisibly)
#'
#' @export
rtx <- function(validation_list,
                 study = "tumor",
                 iterate_through = 5,
                 output_dir,
                 pdf_name = "validation_spectra.pdf",
                 ppm_tolerance = 5,
                 rt_lookup = "range",
                 stick = FALSE,
                 max_i = FALSE,
                 plot_width = 3.9,
                 plot_height = 3.25) {
  
  # Check required parameter
  if (missing(output_dir)) {
    stop("output_dir is required and must be specified")
  }
  
  # Validate iterate_through parameter
  if (iterate_through < 1) {
    stop("iterate_through must be a positive integer")
  }
  
  # Validate rt_lookup parameter
  if (!rt_lookup %in% c("range", "sample")) {
    stop("rt_lookup must be either 'range' or 'sample'")
  }
  
  # Get directory paths from config
  tumor_raw_dir <- config$paths$tumor_raw_dir
  cadaver_raw_dir <- config$paths$cadaver_raw_dir
  
  # Determine directory based on study
  raw_dir <- if (study == "tumor") {
    tumor_raw_dir
  } else if (study == "cadaver") {
    cadaver_raw_dir
  } else {
    stop("study must be 'tumor' or 'cadaver'")
  }
  
  mzml_dir <- file.path(raw_dir, "mzML_validation")
  
  # Helper function to extract chromatogram data
  extract_chromatogram <- function(file_name, target_mzs, rt_range, ppm_tol, use_max) {
    mzml_file <- paste0(file_name, ".mzML")
    mzml_path <- file.path(mzml_dir, mzml_file)
    
    # Check if file exists
    if (!file.exists(mzml_path)) {
      warning(sprintf("mzML file not found: %s - skipping", mzml_path))
      return(NULL)
    }
    
    # Read mzML file
    ms_data <- mzR::openMSfile(mzml_path)
    header_info <- mzR::header(ms_data)
    
    # Function to extract XIC
    extract_xic <- function(ms_data, header_info, target_mz, ppm_tol, rt_range, use_max) {
      ms1_scans <- header_info[header_info$msLevel == 1, ]
      ms1_scans <- ms1_scans[ms1_scans$retentionTime >= rt_range[1] * 60 & 
                             ms1_scans$retentionTime <= rt_range[2] * 60, ]
      
      # Calculate m/z tolerance in Daltons from ppm
      mz_tol_da <- target_mz * ppm_tol / 1e6
      
      xic_data <- lapply(ms1_scans$seqNum, function(scan_num) {
        spectrum <- mzR::peaks(ms_data, scan_num)
        mz_match <- abs(spectrum[, 1] - target_mz) <= mz_tol_da
        
        # Aggregate intensity based on use_max parameter
        if (sum(mz_match) > 0) {
          matched_intensities <- spectrum[mz_match, 2]
          intensity <- if (use_max) max(matched_intensities) else mean(matched_intensities)
        } else {
          intensity <- 0
        }
        
        data.frame(
          scan = scan_num,
          rt = ms1_scans$retentionTime[ms1_scans$seqNum == scan_num] / 60,
          intensity = intensity,
          mz = target_mz
        )
      })
      
      do.call(rbind, xic_data)
    }
    
    # Extract XICs for all target m/z values
    all_xics <- lapply(seq_along(target_mzs), function(i) {
      xic <- extract_xic(ms_data, header_info, target_mzs[i], ppm_tol, rt_range, use_max)
      xic$mz_index <- i - 1
      xic
    })
    
    chromatogram_data <- do.call(rbind, all_xics)
    
    # Close mzML file
    mzR::close(ms_data)
    
    return(chromatogram_data)
  }
  
  # Initialize plot storage list organized by compound
  compound_plots <- list()
  
  # Main iteration loop through validation_list
  cat(sprintf("\nStarting batch validation for %d compounds...\n", nrow(validation_list)))
  
  for (row_idx in 1:nrow(validation_list)) {
    row <- validation_list[row_idx, ]
    
    # Extract row data
    order_num <- row$order
    id_val <- row$id
    short_name <- row$short_name
    
    cat(sprintf("\n[%d/%d] Processing %s (%s)...\n", row_idx, nrow(validation_list), short_name, id_val))
    
    # Initialize list for this compound
    compound_plots[[id_val]] <- list(
      short_name = short_name,
      id = id_val,
      order = order_num,
      plots = list()
    )
    
    # Extract m/z values (mz0, mz1, mz2, mz3)
    target_mzs <- row |>
      select(matches("^mz[0-9]+$")) |>
      unlist() |>
      as.numeric() |>
      na.omit()
    
    if (length(target_mzs) == 0) {
      warning(sprintf("No m/z values found for ID '%s' - skipping", id_val))
      next
    }
    
    # Determine base RT range (used as fallback in "sample" mode or primary in "range" mode)
    rt_range_str <- row$compound_rt_range
    
    if (is.na(rt_range_str)) {
      warning(sprintf("No RT range found for ID '%s' - skipping", id_val))
      next
    }
    
    # Parse RT range string to numeric vector and expand by 0.2 on both ends
    base_rt_range <- eval(parse(text = rt_range_str))
    base_rt_range <- c(base_rt_range[1] - 0.2, base_rt_range[2] + 0.2)
    
    # Parse standards (comma-separated)
    standards <- strsplit(row$standards, ", ")[[1]]
    
    # Get sample files from individual columns (file1-file6)
    all_samples <- c(row$file1, row$file2, row$file3, row$file4, row$file5, row$file6)
    # Remove NAs and limit to iterate_through
    all_samples <- all_samples[!is.na(all_samples)]
    samples_to_process <- all_samples[1:min(iterate_through, length(all_samples))]
    
    # Iterate through sample files
    for (sample_idx in seq_along(samples_to_process)) {
      sample_file <- samples_to_process[sample_idx]
      
      cat(sprintf("  Sample %d/%d: %s\n", sample_idx, length(samples_to_process), sample_file))
      
      # Determine RT range for this sample file
      sample_rt_range <- base_rt_range  # Default to base range
      if (rt_lookup == "sample") {
        # Look up file-specific RT range
        rt_col_name <- paste0("f", sample_idx, "_rt")
        if (rt_col_name %in% names(row)) {
          rt_range_str <- row[[rt_col_name]]
          if (!is.na(rt_range_str)) {
            # Parse the stored RT range (already has buffer from get_rt_range)
            sample_rt_range <- eval(parse(text = rt_range_str))
            cat(sprintf("    Using file-specific RT range: [%.2f, %.2f]\n", sample_rt_range[1], sample_rt_range[2]))
          } else {
            warning(sprintf("No RT range for sample %d in ID '%s' - using base range", sample_idx, id_val))
          }
        } else {
          warning(sprintf("Column '%s' not found for ID '%s' - using base range", rt_col_name, id_val))
        }
      } else {
        cat(sprintf("    Using base RT range: [%.2f, %.2f]\n", sample_rt_range[1], sample_rt_range[2]))
      }
      
      # Extract sample data
      sample_data <- extract_chromatogram(sample_file, target_mzs, sample_rt_range, ppm_tolerance, max_i)
      
      if (is.null(sample_data)) {
        warning(sprintf("Skipping sample %s - file not found", sample_file))
        next
      }
      
      sample_data$type <- "Sample"
      sample_data$plot_intensity <- sample_data$intensity
      
      # Iterate through standards
      for (std_idx in seq_along(standards)) {
        standard_file <- standards[std_idx]
        
        cat(sprintf("    Standard %d/%d: %s\n", std_idx, length(standards), standard_file))
        
        # Standards use the same RT range as the sample for consistent plotting
        # Extract standard data
        standard_data <- extract_chromatogram(standard_file, target_mzs, sample_rt_range, ppm_tolerance, max_i)
        
        if (is.null(standard_data)) {
          warning(sprintf("Skipping standard %s - file not found", standard_file))
          next
        }
        
        standard_data$type <- "Standard"
        standard_data$plot_intensity <- -standard_data$intensity
        
        # Determine y-axis limits (add 5% buffer above max)
        max_sample <- max(abs(sample_data$intensity), na.rm = TRUE)
        max_standard <- max(abs(standard_data$intensity), na.rm = TRUE)
        y_limit <- max(max_sample, max_standard) * 1.05
        
        # Combine data
        combined_data <- bind_rows(sample_data, standard_data)
        
        # Create labels for m/z
        combined_data$mz_label <- sprintf("mz%d: %.4f", 
                                          combined_data$mz_index, 
                                          combined_data$mz)
        
        # Add asterisk to marked mz if it exists in row$asterisk
        if (!is.na(row$asterisk)) {
          marked_mzs <- strsplit(row$asterisk, ", ")[[1]]
          for (marked_mz in marked_mzs) {
            mz_num <- as.numeric(gsub("mz", "", marked_mz))
            combined_data$mz_label <- ifelse(
              combined_data$mz_index == mz_num,
              paste0(combined_data$mz_label, " **\\***"),
              combined_data$mz_label
            )
          }
        }
        
        # Create title and subtitle
        title_text <- sprintf("%s vs. Standard", tools::toTitleCase(study))
        sample_id <- gsub("_[0-9]+$", "", sample_file)
        subtitle_text <- sprintf("Sample: %s  |  Standard: %s\n%s  |  ID: %s  |  %d ppm", 
                                sample_id, standard_file, short_name, id_val, ppm_tolerance)
        
        # Create plot
        if (stick) {
          p <- ggplot(combined_data, aes(x = rt, y = plot_intensity, color = mz_label, group = interaction(mz_label, type))) +
            geom_segment(aes(xend = rt, yend = 0), linewidth = 0.4) +
            geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4)
        } else {
          p <- ggplot(combined_data, aes(x = rt, y = plot_intensity, color = mz_label, group = interaction(mz_label, type))) +
            geom_line(linewidth = 0.4) +
            geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4)
        }
        
        p <- p +
          scale_color_brewer(palette = "Set1") +
          scale_x_continuous(limits = sample_rt_range, expand = expansion(mult = c(0.05, 0.05), add = 0)) +
          scale_y_continuous(
            expand = c(0, 0),
            limits = c(-y_limit, y_limit),
            labels = function(x) abs(x),
            n.breaks = 8
          ) +
          labs(
            title = title_text,
            subtitle = subtitle_text,
            x = "Retention Time (minutes)",
            y = sprintf("\u2190 Standard  |  %s \u2192", tools::toTitleCase(study)),
            color = NULL
          ) +
          coord_cartesian(clip = "off") +
          theme_classic(base_size = 12) +
          theme(
            plot.margin = margin(30, 20, 20, 20),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            legend.position = "top",
            legend.justification = "left",
            legend.direction = "horizontal",
            legend.text = ggtext::element_markdown(size = 5),
            legend.title = element_blank(),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.key = element_rect(fill = "transparent", color = NA),
            legend.key.size = unit(0.25, "cm"),
            legend.key.width = unit(0.25, "cm"),
            legend.spacing.x = unit(0.05, "cm"),
            legend.box.margin = margin(0, 0, 2, 0),
            legend.margin = margin(0, 0, 0, 0),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 10, margin = margin(0, 0, 3, 0)),
            plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 6, color = "black", lineheight = 1.2, margin = margin(0, 0, 3, 0)),
            axis.text.x = element_text(face = "bold", color = "black", size = 8),
            axis.text.y = element_text(face = "bold", color = "black", size = 8),
            axis.title.x = element_text(face = "bold", color = "black", size = 10),
            axis.title.y = element_text(face = "bold", color = "black", size = 10, margin = margin(r = 10)),
            axis.ticks.length = unit(0.15, "cm"),
            axis.line = element_line(color = "black", linewidth = 0.4),
            axis.ticks = element_line(color = "black", linewidth = 0.4)
          ) +
          guides(color = guide_legend(override.aes = list(linewidth = 0.4)))
        
        # Store plot in compound structure with label (for PDF layout)
        plot_label <- sprintf("F%d_S%d", sample_idx, std_idx)
        # Add study prefix (T_ for tumor, C_ for cadaver)
        study_prefix <- if (study == "cadaver") "C_" else "T_"
        plot_tag <- sprintf("%sF%d_S%d_%s", study_prefix, sample_idx, std_idx, id_val)
        
        # Store both the plot and metadata for PDF modification
        compound_plots[[id_val]]$plots[[plot_label]] <- list(
          plot = p,
          sample_id = sample_id,
          standard_file = standard_file,
          plot_tag = plot_tag
        )
        
        cat(sprintf("      Created plot: %s\n", plot_label))
      }
    }
  }
  
  # Compile all compounds into single PDF
  cat(sprintf("\nCompiling %d compounds into PDF...\n", length(compound_plots)))
  
  # Create output path
  pdf_path <- file.path(output_dir, pdf_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Open PDF device with Helvetica font family (closest to Arial for PDF)
  pdf(pdf_path, width = 8.5, height = 11, family = "Helvetica", onefile = TRUE)
  
  # Track if first page to avoid initial blank
  first_page <- TRUE
  
  # Process each compound
  for (compound_id in names(compound_plots)) {
    compound <- compound_plots[[compound_id]]
    plots <- compound$plots
    
    if (length(plots) == 0) next
    
    cat(sprintf("  Adding %s (%s): %d plots\n", compound$short_name, compound$id, length(plots)))
    
    # Sort plot labels by file then standard (F1_S1, F1_S2, F2_S1, F2_S2, ...)
    plot_labels <- names(plots)
    plot_labels <- plot_labels[order(
      as.numeric(gsub("F([0-9]+)_S([0-9]+)", "\\1", plot_labels)),
      as.numeric(gsub("F([0-9]+)_S([0-9]+)", "\\2", plot_labels))
    )]
    
    # Determine grid layout: 2 columns (S1, S2) x N rows (files)
    n_files <- length(unique(gsub("F([0-9]+)_S([0-9]+)", "\\1", plot_labels)))
    
    # Modify plots for PDF: remove title, remove second subtitle line, add red tag to first subtitle
    pdf_plots <- lapply(plot_labels, function(label) {
      plot_info <- plots[[label]]
      p <- plot_info$plot
      sample_id <- plot_info$sample_id
      standard_file <- plot_info$standard_file
      plot_tag <- plot_info$plot_tag
      
      # Create new subtitle with red tag at the end
      new_subtitle <- sprintf("Sample: %s  |  Standard: %s  |  <span style='color:red;'>%s</span>",
                              sample_id, standard_file, plot_tag)
      
      # Modify the plot: remove title and update subtitle
      p_modified <- p + 
        labs(
          title = NULL,
          subtitle = new_subtitle
        ) +
        theme(
          plot.subtitle = ggtext::element_markdown(hjust = 0.5, size = 6, color = "black", lineheight = 1.2, margin = margin(0, 0, 10, 0))
        )
      
      return(ggplotGrob(p_modified))
    })
    
    # Split plots into pages (max 3 rows = 6 plots per page, since 2 columns)
    plots_per_page <- 6  # 3 rows × 2 columns
    n_pages <- ceiling(length(pdf_plots) / plots_per_page)
    
    for (page_num in 1:n_pages) {
      # Get plots for this page
      start_idx <- (page_num - 1) * plots_per_page + 1
      end_idx <- min(page_num * plots_per_page, length(pdf_plots))
      page_plots <- pdf_plots[start_idx:end_idx]
      
      # Create title grob for compound header (add "continued" if not first page)
      title_text <- if (page_num == 1) {
        sprintf("%s (%s)", compound$short_name, compound$id)
      } else {
        sprintf("%s (%s) - continued", compound$short_name, compound$id)
      }
      
      title_grob <- grid::textGrob(
        title_text,
        gp = grid::gpar(fontsize = 14, fontface = "bold"),
        x = 0.5, y = 0.95, just = "top"
      )
      
      # Create grid layout: 2 columns x up to 3 rows
      grid_plot <- gridExtra::arrangeGrob(
        grobs = page_plots,
        ncol = 2,
        top = title_grob
      )
      
      # Draw to PDF (only call grid.newpage() after first page)
      tryCatch({
        if (!first_page) {
          grid::grid.newpage()
        }
        grid::grid.draw(grid_plot)
        first_page <- FALSE
      }, error = function(e) {
        cat(sprintf("  Error drawing compound %s page %d: %s\n", compound$id, page_num, e$message))
      })
    }
  }
  
  # Close PDF
  tryCatch({
    dev.off()
  }, error = function(e) {
    cat(sprintf("Error closing PDF: %s\n", e$message))
  })
  
  cat(sprintf("\nPDF saved: %s\n", pdf_path))
  cat(sprintf("Total plots generated: %d\n", sum(sapply(compound_plots, function(x) length(x$plots)))))
  
  return(invisible(compound_plots))
}

