# Required libraries: mzR, ggplot2, dplyr, gridExtra, ggtext
library(mzR)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(ggtext)

#' Combined Batch Validation: Chromatogram + Spectrum PDF Generator
#'
#' Creates side-by-side plots of RT chromatograms (left) and m/z spectra (right)
#' for each sample/standard combination, compiling them into a single PDF.
#'
#' @param validation_list Tibble with columns: order, id, short_name, monoisotopic, compound_rt_range, 
#'   mz0-mz3, standards (comma-separated), file1-file6, asterisk, f1_rt-f6_rt, s1_rt-s2_rt
#' @param study Character string: "tumor" or "cadaver" to determine directory (default: "tumor")
#' @param iterate_through Integer: how many of the top sample files to process (default: 5)
#' @param output_dir Required character string: output directory for final PDF
#' @param pdf_name Character string: name of output PDF file (default: "validation_combined.pdf")
#' @param ppm_tolerance Numeric mass tolerance in ppm for both RTX and MZX (default: 5)
#' @param rt_tolerance Numeric: RT window in minutes for MZX spectrum extraction (default: 0.1667 = 10 seconds)
#' @param rt_lookup Character: "range" uses compound_rt_range (default), "sample" uses file-specific RT ranges
#' @param display_mode Character: "max" uses maximum intensity in RT window, "avg" uses average (default: "max")
#' @param stick Logical, whether to plot RTX as vertical sticks (default: FALSE)
#' @param max_i Logical, whether to use maximum intensity for RTX (default: FALSE)
#' @param plot_width Plot width in inches (default: 3.9)
#' @param plot_height Plot height in inches (default: 3.25)
#'
#' @return Named list of all plots (invisibly)
#'
#' @export
rtx_mzx <- function(validation_list,
                    study = "tumor",
                    iterate_through = 5,
                    output_dir,
                    pdf_name = "validation_combined.pdf",
                    ppm_tolerance = 5,
                    rt_tolerance = 10/60,
                    rt_lookup = "range",
                    display_mode = "max",
                    stick = FALSE,
                    max_i = FALSE,
                    plot_width = 3.9,
                    plot_height = 3.25) {
  
  # Check required parameter
  if (missing(output_dir)) {
    stop("output_dir is required and must be specified")
  }
  
  # Validate parameters
  if (iterate_through < 1) {
    stop("iterate_through must be a positive integer")
  }
  
  if (!rt_lookup %in% c("range", "sample")) {
    stop("rt_lookup must be either 'range' or 'sample'")
  }
  
  if (!display_mode %in% c("max", "avg")) {
    stop("display_mode must be either 'max' or 'avg'")
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
  
  # Helper function to extract chromatogram data (for RTX)
  extract_chromatogram <- function(file_name, target_mzs, rt_range, ppm_tol, use_max) {
    mzml_file <- paste0(file_name, ".mzML")
    mzml_path <- file.path(mzml_dir, mzml_file)
    
    if (!file.exists(mzml_path)) {
      warning(sprintf("mzML file not found: %s - skipping", mzml_path))
      return(NULL)
    }
    
    ms_data <- mzR::openMSfile(mzml_path)
    header_info <- mzR::header(ms_data)
    
    extract_xic <- function(ms_data, header_info, target_mz, ppm_tol, rt_range, use_max) {
      ms1_scans <- header_info[header_info$msLevel == 1, ]
      ms1_scans <- ms1_scans[ms1_scans$retentionTime >= rt_range[1] * 60 & 
                             ms1_scans$retentionTime <= rt_range[2] * 60, ]
      
      mz_tol_da <- target_mz * ppm_tol / 1e6
      
      xic_data <- lapply(ms1_scans$seqNum, function(scan_num) {
        spectrum <- mzR::peaks(ms_data, scan_num)
        mz_match <- abs(spectrum[, 1] - target_mz) <= mz_tol_da
        
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
    
    all_xics <- lapply(seq_along(target_mzs), function(i) {
      xic <- extract_xic(ms_data, header_info, target_mzs[i], ppm_tol, rt_range, use_max)
      xic$mz_index <- i - 1
      xic
    })
    
    chromatogram_data <- do.call(rbind, all_xics)
    mzR::close(ms_data)
    
    return(chromatogram_data)
  }
  
  # Helper function to extract spectrum at RT (for MZX)
  extract_spectrum_at_rt <- function(file_name, target_rt, rt_tol, target_mzs, ppm_filter, display_mode) {
    mzml_file <- paste0(file_name, ".mzML")
    mzml_path <- file.path(mzml_dir, mzml_file)
    
    if (!file.exists(mzml_path)) {
      warning(sprintf("mzML file not found: %s - skipping", mzml_path))
      return(NULL)
    }
    
    ms_data <- mzR::openMSfile(mzml_path)
    header_info <- mzR::header(ms_data)
    
    ms1_scans <- header_info[header_info$msLevel == 1, ]
    rt_window <- c((target_rt - rt_tol) * 60, (target_rt + rt_tol) * 60)
    scans_in_window <- ms1_scans[ms1_scans$retentionTime >= rt_window[1] & 
                                   ms1_scans$retentionTime <= rt_window[2], ]
    
    if (nrow(scans_in_window) == 0) {
      warning(sprintf("No scans found in RT window [%.3f, %.3f] for %s", 
                      target_rt - rt_tol, target_rt + rt_tol, file_name))
      mzR::close(ms_data)
      return(NULL)
    }
    
    all_spectra <- lapply(scans_in_window$seqNum, function(scan_num) {
      mzR::peaks(ms_data, scan_num)
    })
    
    combined_mz <- unlist(lapply(all_spectra, function(s) s[, 1]))
    combined_intensity <- unlist(lapply(all_spectra, function(s) s[, 2]))
    
    if (display_mode == "max") {
      spectrum_df <- data.frame(mz = combined_mz, intensity = combined_intensity) |>
        group_by(mz = round(mz, 4)) |>
        summarize(intensity = max(intensity), .groups = "drop")
    } else {
      spectrum_df <- data.frame(mz = combined_mz, intensity = combined_intensity) |>
        group_by(mz = round(mz, 4)) |>
        summarize(intensity = mean(intensity), .groups = "drop")
    }
    
    filtered_spectrum <- do.call(rbind, lapply(seq_along(target_mzs), function(i) {
      target_mz <- target_mzs[i]
      mz_window <- target_mz * ppm_filter / 1e6
      mz_range <- c(target_mz - mz_window, target_mz + mz_window)
      
      in_range <- spectrum_df$mz >= mz_range[1] & spectrum_df$mz <= mz_range[2]
      if (sum(in_range) == 0) return(NULL)
      
      df <- spectrum_df[in_range, ]
      df$mz_index <- i - 1
      df$target_mz <- target_mz
      df
    }))
    
    mzR::close(ms_data)
    return(filtered_spectrum)
  }
  
  # Initialize plot storage
  compound_plots <- list()
  
  # Main iteration loop
  cat(sprintf("\nStarting combined RTX+MZX validation for %d compounds...\n", nrow(validation_list)))
  
  for (row_idx in 1:nrow(validation_list)) {
    row <- validation_list[row_idx, ]
    
    order_num <- row$order
    id_val <- row$id
    short_name <- row$short_name
    
    cat(sprintf("\n[%d/%d] Processing %s (%s)...\n", row_idx, nrow(validation_list), short_name, id_val))
    
    compound_plots[[id_val]] <- list(
      short_name = short_name,
      id = id_val,
      order = order_num,
      plots = list()
    )
    
    # Extract m/z values
    target_mzs <- row |>
      select(matches("^mz[0-9]+$")) |>
      unlist() |>
      as.numeric() |>
      na.omit()
    
    if (length(target_mzs) == 0) {
      warning(sprintf("No m/z values found for ID '%s' - skipping", id_val))
      next
    }
    
    # Get base RT range and RT
    rt_range_str <- row$compound_rt_range
    if (is.na(rt_range_str)) {
      warning(sprintf("No RT range found for ID '%s' - skipping", id_val))
      next
    }
    
    base_rt_range <- eval(parse(text = rt_range_str))
    base_rt_range_expanded <- c(base_rt_range[1] - 0.2, base_rt_range[2] + 0.2)
    base_rt <- mean(base_rt_range)
    
    # Parse standards
    standards <- strsplit(row$standards, ", ")[[1]]
    
    # Get sample files
    all_samples <- c(row$file1, row$file2, row$file3, row$file4, row$file5, row$file6)
    all_samples <- all_samples[!is.na(all_samples)]
    samples_to_process <- all_samples[1:min(iterate_through, length(all_samples))]
    
    # Iterate through sample files
    for (sample_idx in seq_along(samples_to_process)) {
      sample_file <- samples_to_process[sample_idx]
      
      cat(sprintf("  Sample %d/%d: %s\n", sample_idx, length(samples_to_process), sample_file))
      
      # Determine RT range and RT point for this sample
      sample_rt_range <- base_rt_range_expanded
      sample_rt <- base_rt
      
      if (rt_lookup == "sample") {
        rt_col_name <- paste0("f", sample_idx, "_rt")
        rt_range_col_name <- paste0("f", sample_idx, "_rt_range")
        
        if (rt_col_name %in% names(row) && !is.na(row[[rt_col_name]])) {
          sample_rt <- row[[rt_col_name]]
        }
        
        if (rt_range_col_name %in% names(row) && !is.na(row[[rt_range_col_name]])) {
          rt_range_str <- row[[rt_range_col_name]]
          sample_rt_range <- eval(parse(text = rt_range_str))
          cat(sprintf("    Using file-specific RT: %.3f min, range: [%.2f, %.2f]\n", 
                      sample_rt, sample_rt_range[1], sample_rt_range[2]))
        }
      }
      
      # Extract BOTH chromatogram and spectrum for sample
      sample_chrom <- extract_chromatogram(sample_file, target_mzs, sample_rt_range, ppm_tolerance, max_i)
      sample_spec <- extract_spectrum_at_rt(sample_file, sample_rt, rt_tolerance, target_mzs, ppm_tolerance, display_mode)
      
      if (is.null(sample_chrom) || is.null(sample_spec)) {
        warning(sprintf("Skipping sample %s - data extraction failed", sample_file))
        next
      }
      
      sample_chrom$type <- "Sample"
      sample_chrom$plot_intensity <- sample_chrom$intensity
      
      sample_spec$type <- "Sample"
      sample_spec$plot_intensity <- sample_spec$intensity
      
      # Iterate through standards
      for (std_idx in seq_along(standards)) {
        standard_file <- standards[std_idx]
        
        cat(sprintf("    Standard %d/%d: %s\n", std_idx, length(standards), standard_file))
        
        # Extract BOTH chromatogram and spectrum for standard
        standard_chrom <- extract_chromatogram(standard_file, target_mzs, sample_rt_range, ppm_tolerance, max_i)
        standard_spec <- extract_spectrum_at_rt(standard_file, sample_rt, rt_tolerance, target_mzs, ppm_tolerance, display_mode)
        
        if (is.null(standard_chrom) || is.null(standard_spec)) {
          warning(sprintf("Skipping standard %s - data extraction failed", standard_file))
          next
        }
        
        standard_chrom$type <- "Standard"
        standard_chrom$plot_intensity <- -standard_chrom$intensity
        
        standard_spec$type <- "Standard"
        standard_spec$plot_intensity <- -standard_spec$intensity
        
        # === CREATE RTX PLOT (Chromatogram) ===
        max_sample_chrom <- max(abs(sample_chrom$intensity), na.rm = TRUE)
        max_standard_chrom <- max(abs(standard_chrom$intensity), na.rm = TRUE)
        y_limit_chrom <- max(max_sample_chrom, max_standard_chrom) * 1.05
        
        combined_chrom <- bind_rows(sample_chrom, standard_chrom)
        
        # Filter to only rows with non-zero intensity (detected fragments)
        combined_chrom <- combined_chrom |>
          filter(intensity > 0)
        
        if (nrow(combined_chrom) == 0) {
          warning(sprintf("No fragments detected for %s in %s vs %s", id_val, sample_file, standard_file))
          next
        }
        
        combined_chrom$mz_label <- sprintf("mz%d: %.4f", combined_chrom$mz_index, combined_chrom$mz)
        
        # Add asterisks
        if (!is.na(row$asterisk)) {
          marked_mzs <- strsplit(row$asterisk, ", ")[[1]]
          for (marked_mz in marked_mzs) {
            mz_num <- as.numeric(gsub("mz", "", marked_mz))
            combined_chrom$mz_label <- ifelse(
              combined_chrom$mz_index == mz_num,
              paste0(combined_chrom$mz_label, " **\\***"),
              combined_chrom$mz_label
            )
          }
        }
        
        sample_id <- gsub("_[0-9]+$", "", sample_file)
        
        if (stick) {
          p_rtx <- ggplot(combined_chrom, aes(x = rt, y = plot_intensity, color = mz_label, group = interaction(mz_label, type))) +
            geom_segment(aes(xend = rt, yend = 0), linewidth = 0.4) +
            geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4)
        } else {
          p_rtx <- ggplot(combined_chrom, aes(x = rt, y = plot_intensity, color = mz_label, group = interaction(mz_label, type))) +
            geom_line(linewidth = 0.4) +
            geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4)
        }
        
        p_rtx <- p_rtx +
          scale_color_viridis_d(option = "turbo", end = 0.9) +
          scale_x_continuous(limits = sample_rt_range, expand = expansion(mult = c(0.05, 0.05), add = 0)) +
          scale_y_continuous(
            expand = c(0, 0),
            limits = c(-y_limit_chrom, y_limit_chrom),
            labels = function(x) abs(x),
            n.breaks = 8
          ) +
          labs(
            title = "RTX (Chromatogram)",
            subtitle = sprintf("Sample: %s  |  Standard: %s", sample_id, standard_file),
            x = "Retention Time (min)",
            y = sprintf("\u2190 Std  |  %s \u2192", tools::toTitleCase(study)),
            color = NULL
          ) +
          coord_cartesian(clip = "off") +
          theme_classic(base_size = 12) +
          theme(
            plot.margin = margin(20, 10, 20, 20),
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
            plot.title = element_text(hjust = 0.5, face = "bold", size = 9, margin = margin(0, 0, 2, 0)),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5, face = "italic", size = 6, 
                                                     color = "black", lineheight = 1.2, margin = margin(0, 0, 3, 0)),
            axis.text.x = element_text(face = "bold", color = "black", size = 7),
            axis.text.y = element_text(face = "bold", color = "black", size = 7),
            axis.title.x = element_text(face = "bold", color = "black", size = 9),
            axis.title.y = element_text(face = "bold", color = "black", size = 9, margin = margin(r = 8)),
            axis.ticks.length = unit(0.15, "cm"),
            axis.line = element_line(color = "black", linewidth = 0.4),
            axis.ticks = element_line(color = "black", linewidth = 0.4)
          ) +
          guides(color = guide_legend(override.aes = list(linewidth = 0.4)))
        
        # === CREATE MZX PLOT (Spectrum) ===
        max_sample_spec <- max(abs(sample_spec$intensity), na.rm = TRUE)
        max_standard_spec <- max(abs(standard_spec$intensity), na.rm = TRUE)
        y_limit_spec <- max(max_sample_spec, max_standard_spec) * 1.05
        
        combined_spec <- bind_rows(sample_spec, standard_spec)
        
        # Filter to only rows with non-zero intensity (detected fragments)
        combined_spec <- combined_spec |>
          filter(intensity > 0)
        
        if (nrow(combined_spec) == 0) {
          warning(sprintf("No fragments detected in spectrum for %s in %s vs %s", id_val, sample_file, standard_file))
          next
        }
        
        combined_spec$mz_label <- sprintf("mz%d: %.4f", combined_spec$mz_index, combined_spec$target_mz)
        
        # Add asterisks
        if (!is.na(row$asterisk)) {
          marked_mzs <- strsplit(row$asterisk, ", ")[[1]]
          for (marked_mz in marked_mzs) {
            mz_num <- as.numeric(gsub("mz", "", marked_mz))
            combined_spec$mz_label <- ifelse(
              combined_spec$mz_index == mz_num,
              paste0(combined_spec$mz_label, " **\\***"),
              combined_spec$mz_label
            )
          }
        }
        
        p_mzx <- ggplot(combined_spec, aes(x = mz, y = plot_intensity, color = mz_label, group = interaction(mz_label, type))) +
          geom_segment(aes(xend = mz, yend = 0), linewidth = 0.4) +
          geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
          scale_color_viridis_d(option = "turbo", end = 0.9) +
          scale_x_continuous(expand = expansion(mult = c(0.05, 0.05), add = 0)) +
          scale_y_continuous(
            expand = c(0, 0),
            limits = c(-y_limit_spec, y_limit_spec),
            labels = function(x) abs(x),
            n.breaks = 8
          ) +
          labs(
            title = "MZX (Spectrum)",
            subtitle = sprintf("RT: %.3f min  |  %d ppm", sample_rt, ppm_tolerance),
            x = "m/z",
            y = sprintf("\u2190 Std  |  %s \u2192", tools::toTitleCase(study)),
            color = NULL
          ) +
          coord_cartesian(clip = "off") +
          theme_classic(base_size = 12) +
          theme(
            plot.margin = margin(20, 20, 20, 10),
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
            plot.title = element_text(hjust = 0.5, face = "bold", size = 9, margin = margin(0, 0, 2, 0)),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5, face = "italic", size = 6, 
                                                     color = "black", lineheight = 1.2, margin = margin(0, 0, 3, 0)),
            axis.text.x = element_text(face = "bold", color = "black", size = 7),
            axis.text.y = element_text(face = "bold", color = "black", size = 7),
            axis.title.x = element_text(face = "bold", color = "black", size = 9),
            axis.title.y = element_text(face = "bold", color = "black", size = 9, margin = margin(r = 8)),
            axis.ticks.length = unit(0.15, "cm"),
            axis.line = element_line(color = "black", linewidth = 0.4),
            axis.ticks = element_line(color = "black", linewidth = 0.4)
          ) +
          guides(color = guide_legend(override.aes = list(linewidth = 0.4)))
        
        # Store both plots
        plot_label <- sprintf("F%d_S%d", sample_idx, std_idx)
        study_prefix <- if (study == "cadaver") "C_" else "T_"
        plot_tag_rtx <- sprintf("%sF%d_S%d_%s", study_prefix, sample_idx, std_idx, id_val)
        plot_tag_mzx <- sprintf("%sF%d_S%d_%s_spectra", study_prefix, sample_idx, std_idx, id_val)
        
        compound_plots[[id_val]]$plots[[plot_label]] <- list(
          plot_rtx = p_rtx,
          plot_mzx = p_mzx,
          sample_id = sample_id,
          standard_file = standard_file,
          plot_tag_rtx = plot_tag_rtx,
          plot_tag_mzx = plot_tag_mzx
        )
        
        cat(sprintf("      Created plots: %s\n", plot_label))
      }
    }
  }
  
  # Compile into PDF
  cat(sprintf("\nCompiling %d compounds into PDF...\n", length(compound_plots)))
  
  pdf_path <- file.path(output_dir, pdf_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  pdf(pdf_path, width = 8.5, height = 11, family = "Helvetica", onefile = TRUE)
  
  first_page <- TRUE
  
  for (compound_id in names(compound_plots)) {
    compound <- compound_plots[[compound_id]]
    plots <- compound$plots
    
    if (length(plots) == 0) next
    
    cat(sprintf("  Adding %s (%s): %d plot pairs\n", compound$short_name, compound$id, length(plots)))
    
    # Sort plot labels
    plot_labels <- names(plots)
    plot_labels <- plot_labels[order(
      as.numeric(gsub("F([0-9]+)_S([0-9]+)", "\\1", plot_labels)),
      as.numeric(gsub("F([0-9]+)_S([0-9]+)", "\\2", plot_labels))
    )]
    
    # Modify plots for PDF: add red tags to both RTX and MZX subtitles
    pdf_plot_pairs <- lapply(plot_labels, function(label) {
      plot_info <- plots[[label]]
      p_rtx <- plot_info$plot_rtx
      p_mzx <- plot_info$plot_mzx
      sample_id <- plot_info$sample_id
      standard_file <- plot_info$standard_file
      plot_tag_rtx <- plot_info$plot_tag_rtx
      plot_tag_mzx <- plot_info$plot_tag_mzx
      
      # Modify RTX plot: add red tag (chromatogram tag)
      new_subtitle_rtx <- sprintf("Sample: %s  |  Standard: %s  |  <span style='color:red;'>%s</span>",
                                  sample_id, standard_file, plot_tag_rtx)
      
      p_rtx <- p_rtx +
        labs(subtitle = new_subtitle_rtx) +
        theme(
          plot.subtitle = ggtext::element_markdown(hjust = 0.5, face = "italic", size = 6,
                                                   color = "black", lineheight = 1.2, margin = margin(0, 0, 3, 0))
        )
      
      # Modify MZX plot: add red tag (spectrum tag with _spectra suffix)
      # Extract RT from existing subtitle
      existing_subtitle <- p_mzx$labels$subtitle
      new_subtitle_mzx <- sprintf("%s  |  <span style='color:red;'>%s</span>",
                                  existing_subtitle, plot_tag_mzx)
      
      p_mzx <- p_mzx +
        labs(subtitle = new_subtitle_mzx) +
        theme(
          plot.subtitle = ggtext::element_markdown(hjust = 0.5, face = "italic", size = 6,
                                                   color = "black", lineheight = 1.2, margin = margin(0, 0, 3, 0))
        )
      
      list(rtx = p_rtx, mzx = p_mzx)
    })
    
    # Create pages: 3 rows per page, each row has RTX (left) and MZX (right)
    pairs_per_page <- 3
    n_pages <- ceiling(length(pdf_plot_pairs) / pairs_per_page)
    
    for (page_num in 1:n_pages) {
      start_idx <- (page_num - 1) * pairs_per_page + 1
      end_idx <- min(page_num * pairs_per_page, length(pdf_plot_pairs))
      page_pairs <- pdf_plot_pairs[start_idx:end_idx]
      
      # Flatten into alternating RTX, MZX for grid (row1: RTX, MZX; row2: RTX, MZX; row3: RTX, MZX)
      page_plots <- unlist(lapply(page_pairs, function(pair) list(pair$rtx, pair$mzx)), recursive = FALSE)
      
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
      
      # Create grid: 2 columns (RTX, MZX) x up to 3 rows
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
  
  tryCatch({
    dev.off()
    cat(sprintf("\nPDF saved to: %s\n", pdf_path))
  }, error = function(e) {
    cat(sprintf("Error closing PDF: %s\n", e$message))
  })
  
  invisible(compound_plots)
}
