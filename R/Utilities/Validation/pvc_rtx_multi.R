#' Plot Validation Chromatogram Multi-Panel (Stacked Individual Plots)
#'
#' Creates separate stacked chromatogram plots for each m/z trace with independent y-axes.
#' All plots are combined into a single PNG output using cowplot.
#'
#' @param id Character string matching the 'id' column in reference_table
#' @param file_name_sample Character string of the sample raw file name (without .mzML extension)
#' @param file_name_standard Character string of the standard raw file name (without .mzML extension)
#' @param reference_table Tibble with columns: id, cas, mz0-mz3, tumor_rt_range, cadaver_rt_range (default: mz_reference_table)
#' @param study Character string: "tumor" or "cadaver" to determine directory and RT range (default: "tumor")
#' @param file_list Tibble with columns: ID, type, files (for sample metadata) (default: file_list)
#' @param selected_features Tibble with columns: source, id, subid, name, tmz, trt, cas (default: selected_gc2_features)
#' @param source Character string: "IARC" or "quant" to resolve which tmz/trt to use if multiple exist (default: NULL uses first)
#' @param ppm_tolerance Numeric mass tolerance in ppm (default: 5)
#' @param rtr Optional numeric vector c(min, max) to override dynamic RT range
#' @param standard Logical, whether to include standard plot (default: TRUE). If FALSE, only plots sample chromatogram
#' @param show_lib_rt Logical, whether to show library RT in subtitle (default: TRUE). If FALSE, omits "Lib RT = X.XX" from subtitle
#' @param stick Logical, whether to plot as vertical sticks (default: FALSE). If TRUE, uses geom_segment instead of geom_line
#' @param max_i Logical, whether to use maximum intensity (default: FALSE). If TRUE, uses max instead of mean when aggregating peaks within ppm tolerance
#' @param mz_manual Numeric vector of m/z values to plot (default: NULL). If provided, overrides reference table m/z lookup. Use c(132.0453, 144.0433) to plot specific m/z values
#' @param mz_display Character "all" or numeric vector specifying which m/z indices to display (default: "all"). Use c(0,1,2) to show only mz0, mz1, mz2
#' @param spectrum Logical, whether to include full mass spectrum plot below chromatograms (default: FALSE). Shows all peaks in RT range as sticks from m/z 85-850
#' @param rt Numeric retention time to use for spectrum extraction when spectrum = TRUE (default: NULL). If NULL, uses center of rtr range. Units specified by rt_unit
#' @param rt_unit Character string specifying units for rt parameter: "min" for minutes (default) or "s" for seconds
#' @param spectra_labels Logical, whether to label peaks in spectrum plot (default: FALSE). If TRUE, labels top peaks based on label_tol
#' @param label_tol Numeric value between 0 and 1 specifying intensity percentile threshold for labeling (default: 0.99). Value of 0.99 labels top 1% of peaks
#' @param png_name Optional custom filename for PNG output (without .png extension). If NULL (default), no file is saved
#' @param output_dir Directory to save PNG file (default: "Outputs/Spectra")
#' @param width Plot width in inches (default: 6)
#' @param height Plot height per panel in inches (default: 2)
#'
#' @return Combined plot object (saves to PNG only if png_name is provided)
#'
#' @export
pvc_rtx_multi <- function(id,
                          file_name_sample,
                          file_name_standard = NULL,
                          reference_table = NULL,
                          study = "tumor",
                          file_list = NULL,
                          selected_features = NULL,
                          source = NULL,
                          ppm_tolerance = 5,
                          rtr = NULL,
                          standard = TRUE,
                          show_lib_rt = TRUE,
                          stick = FALSE,
                          max_i = FALSE,
                          mz_manual = NULL,
                          mz_display = "all",
                          spectrum = FALSE,
                          rt = NULL,
                          rt_unit = "min",
                          spectra_labels = FALSE,
                          label_tol = 0.99,
                          png_name = NULL,
                          output_dir = "Outputs/Spectra/rtx",
                          width = 6,
                          height = 2) {
  
  # Use defaults from environment if not provided
  if (is.null(reference_table)) {
    reference_table <- mz_reference_table
  }
  if (is.null(file_list)) {
    file_list <- get("file_list", envir = parent.frame())
  }
  if (is.null(selected_features)) {
    selected_features <- get("selected_gc2_features", envir = parent.frame())
  }
  
  # Get directory paths from config
  tumor_raw_dir <- config$paths$tumor_raw_dir
  cadaver_raw_dir <- config$paths$cadaver_raw_dir
  
  # Look up target m/z and retention time from selected_features
  feature_match <- selected_features |>
    filter(id == !!id)
  
  # If source specified, filter by it (handles duplicates like o-Toluidine)
  if (!is.null(source)) {
    feature_match <- feature_match |>
      filter(source == !!source)
  }
  
  if (nrow(feature_match) == 0) {
    warning(sprintf("ID '%s' not found in selected_features. No target m/z or RT marker will be added.", id))
    target_tmz <- NULL
    library_trt <- NULL
  } else {
    # Use first match if multiple remain
    target_tmz <- feature_match$tmz[1]
    library_trt <- feature_match$trt[1]
  }
  
  # Determine target m/z values - use manual if provided, otherwise look up from reference table
  if (!is.null(mz_manual)) {
    # Use manually provided m/z values
    if (!is.numeric(mz_manual)) {
      stop("mz_manual must be a numeric vector")
    }
    target_mzs <- mz_manual
    ref_row <- NULL  # Set to NULL to indicate manual mode
  } else {
    # Get reference data for this feature
    ref_row <- reference_table |>
      filter(id == !!id)
    
    if (nrow(ref_row) == 0) {
      stop(sprintf("ID '%s' not found in reference_table", id))
    }
    
    # Extract target m/z values (mz0, mz1, mz2, mz3)
    target_mzs <- ref_row |>
      select(matches("^mz[0-9]+$")) |>
      unlist() |>
      as.numeric() |>
      na.omit()
    
    if (length(target_mzs) == 0) {
      stop(sprintf("No m/z values found for ID '%s'", id))
    }
  }
  
  # Determine RT range - use manual if provided, otherwise dynamic
  manual_rt_override <- !is.null(rtr)
  
  if (!is.null(rtr)) {
    rt_range <- rtr
  } else {
    # If using manual m/z values, RT range must be provided
    if (!is.null(mz_manual)) {
      stop("rtr (RT range) must be provided when using mz_manual")
    }
    
    # Determine RT range based on study
    rt_range_str <- if (study == "tumor") {
      ref_row$tumor_rt_range
    } else if (study == "cadaver") {
      ref_row$cadaver_rt_range
    } else {
      stop("study must be 'tumor' or 'cadaver'")
    }
    
    if (is.na(rt_range_str)) {
      stop(sprintf("No %s RT range found for ID '%s'", study, id))
    }
    
    # Parse RT range string to numeric vector and expand by 0.2 on both ends
    rt_range <- eval(parse(text = rt_range_str))
    rt_range <- c(rt_range[1] - 0.2, rt_range[2] + 0.2)
  }
  
  # Determine directory based on study
  raw_dir <- if (study == "tumor") {
    tumor_raw_dir
  } else {
    cadaver_raw_dir
  }
  
  mzml_dir <- file.path(raw_dir, "mzML_validation")
  
  # Helper function to extract chromatogram data
  extract_chromatogram <- function(file_name) {
    mzml_file <- paste0(file_name, ".mzML")
    mzml_path <- file.path(mzml_dir, mzml_file)
    
    # Check if file exists
    if (!file.exists(mzml_path)) {
      stop(sprintf("mzML file not found: %s", mzml_path))
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
      xic <- extract_xic(ms_data, header_info, target_mzs[i], ppm_tolerance, rt_range, max_i)
      xic$mz_index <- i - 1
      xic
    })
    
    chromatogram_data <- do.call(rbind, all_xics)
    
    # Close mzML file
    mzR::close(ms_data)
    
    return(chromatogram_data)
  }
  
  # Helper function to extract full spectrum data (if spectrum = TRUE)
  extract_full_spectrum <- function(file_name, target_rt, rt_unit) {
    mzml_file <- paste0(file_name, ".mzML")
    mzml_path <- file.path(mzml_dir, mzml_file)
    
    if (!file.exists(mzml_path)) {
      stop(sprintf("mzML file not found: %s", mzml_path))
    }
    
    ms_data <- mzR::openMSfile(mzml_path)
    header_info <- mzR::header(ms_data)
    
    # Convert target RT to seconds (internal storage is always seconds)
    target_rt_seconds <- if (rt_unit == "s") target_rt else target_rt * 60
    
    # Get MS1 scans and find the one closest to the target RT
    ms1_scans <- header_info[header_info$msLevel == 1, ]
    ms1_scans <- ms1_scans[ms1_scans$retentionTime >= rt_range[1] * 60 & 
                           ms1_scans$retentionTime <= rt_range[2] * 60, ]
    
    # Find scan closest to target RT
    closest_scan_idx <- which.min(abs(ms1_scans$retentionTime - target_rt_seconds))
    scan_num <- ms1_scans$seqNum[closest_scan_idx]
    actual_rt <- ms1_scans$retentionTime[closest_scan_idx] / 60
    
    # Extract spectrum from the single scan at midpoint
    spectrum <- mzR::peaks(ms_data, scan_num)
    spectrum_data <- data.frame(mz = spectrum[, 1], intensity = spectrum[, 2]) |>
      filter(mz >= 85, mz <= 850)
    
    mzR::close(ms_data)
    return(list(data = spectrum_data, actual_rt = actual_rt))
  }
  
  # Extract data for sample
  sample_data <- extract_chromatogram(file_name_sample)
  sample_data$type <- "Sample"
  sample_data$plot_intensity <- sample_data$intensity  # Keep positive
  
  # Handle standard data based on standard parameter
  if (standard) {
    if (is.null(file_name_standard)) {
      stop("file_name_standard must be provided when standard = TRUE")
    }
    standard_data <- extract_chromatogram(file_name_standard)
    standard_data$type <- "Standard"
    standard_data$plot_intensity <- -standard_data$intensity  # Flip to negative for plotting
    
    # Combine data
    combined_data <- bind_rows(sample_data, standard_data)
  } else {
    # Only sample data
    combined_data <- sample_data
  }
  
  # Filter m/z indices based on mz_display parameter
  if (!identical(mz_display, "all")) {
    if (!is.numeric(mz_display)) {
      stop("mz_display must be either 'all' or a numeric vector")
    }
    combined_data <- combined_data |>
      filter(mz_index %in% mz_display)
  }
  
  # Check if any data remains after filtering
  if (nrow(combined_data) == 0) {
    stop("No data to plot after applying mz_display filter")
  }
  
  # Create labels with asterisk for target m/z
  combined_data$mz_label <- sprintf("mz%d: %.4f", 
                                    combined_data$mz_index, 
                                    combined_data$mz)
  
  # Add asterisk to target m/z if it exists
  if (!is.null(target_tmz)) {
    combined_data$mz_label <- ifelse(
      abs(combined_data$mz - target_tmz) < 0.01,  # Match with 0.01 tolerance
      paste0(combined_data$mz_label, "*"),
      combined_data$mz_label
    )
  }
  
  # Get sample metadata
  sample_metadata <- file_list |>
    filter(grepl(file_name_sample, files, fixed = TRUE))
  
  if (nrow(sample_metadata) == 0) {
    sample_id <- gsub("_[0-9]+$", "", file_name_sample)
    sample_type <- "Unknown"
  } else {
    sample_id <- sample_metadata$ID[1]
    sample_type <- sample_metadata$type[1]
    sample_type <- tools::toTitleCase(sample_type)
  }
  
  # Create dynamic title and subtitle based on whether standard is included
  study_label <- tools::toTitleCase(study)
  
  # Get display name - use from ref_row if available, otherwise use "Custom m/z"
  display_name <- if (!is.null(ref_row)) ref_row$short_display_name else "Custom m/z"
  
  if (standard) {
    # Get standard ID
    standard_id <- gsub("_[0-9]+$", "", file_name_standard)
    title_text <- sprintf("%s vs. Standard", study_label)
    
    # Create subtitle with metadata (add library RT if available and show_lib_rt is TRUE)
    if (!is.null(library_trt) && show_lib_rt) {
      subtitle_text <- sprintf("Sample: %s (%s)  |  Standard: %s\n%s  |  ID: %s  |  %d ppm  |  Lib RT = %.2f", 
                              sample_id, file_name_sample, file_name_standard, display_name, id, ppm_tolerance, library_trt)
    } else {
      subtitle_text <- sprintf("Sample: %s (%s)  |  Standard: %s\n%s  |  ID: %s  |  %d ppm", 
                              sample_id, file_name_sample, file_name_standard, display_name, id, ppm_tolerance)
    }
  } else {
    title_text <- sprintf("%s Chromatogram", study_label)
    
    # Create subtitle with metadata (add library RT if available and show_lib_rt is TRUE)
    if (!is.null(library_trt) && show_lib_rt) {
      subtitle_text <- sprintf("Sample: %s (%s)\n%s  |  ID: %s  |  %d ppm  |  Lib RT = %.2f", 
                              sample_id, file_name_sample, display_name, id, ppm_tolerance, library_trt)
    } else {
      subtitle_text <- sprintf("Sample: %s (%s)\n%s  |  ID: %s  |  %d ppm", 
                              sample_id, file_name_sample, display_name, id, ppm_tolerance)
    }
  }
  
  # Get unique m/z indices for creating separate plots
  mz_indices <- unique(combined_data$mz_index)
  
  # Define color palette for different m/z traces (use RColorBrewer Set1 colors)
  color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
  
  # Create individual plots for each m/z
  plot_list <- lapply(seq_along(mz_indices), function(i) {
    mz_idx <- mz_indices[i]
    
    # Filter data for this m/z
    mz_data <- combined_data |>
      filter(mz_index == mz_idx)
    
    # Get the label for this m/z
    mz_label_text <- unique(mz_data$mz_label)[1]
    
    # Assign color based on position in mz_indices (cycles through palette)
    plot_color <- color_palette[((i - 1) %% length(color_palette)) + 1]
    
    # Determine y-axis limits for this specific m/z
    if (standard) {
      max_sample <- max(abs(mz_data$intensity[mz_data$type == "Sample"]), na.rm = TRUE)
      max_standard <- max(abs(mz_data$intensity[mz_data$type == "Standard"]), na.rm = TRUE)
      y_max_raw <- max(max_sample, max_standard)
    } else {
      y_max_raw <- max(abs(mz_data$intensity), na.rm = TRUE)
    }
    
    # Round y-axis limit to nice round number
    # Get the order of magnitude
    magnitude <- 10^floor(log10(y_max_raw))
    # Round up to the nearest "nice" value (1, 2, 4, 5, or 10 times magnitude)
    nice_values <- c(1, 2, 4, 5, 10) * magnitude
    y_limit <- min(nice_values[nice_values >= y_max_raw])
    
    # Create plot with conditional formatting based on standard and stick parameters
    if (standard) {
      # Mirrored plot
      if (stick) {
        p <- ggplot(mz_data, aes(x = rt, y = plot_intensity, group = type)) +
          geom_segment(aes(xend = rt, yend = 0), linewidth = 0.8, color = plot_color) +
          geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8)
      } else {
        p <- ggplot(mz_data, aes(x = rt, y = plot_intensity, group = type)) +
          geom_line(linewidth = 0.8, color = plot_color) +
          geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8)
      }
      p <- p +
        {if (manual_rt_override) {
          scale_x_continuous(limits = rt_range, breaks = seq(floor(rt_range[1]*200)/200, ceiling(rt_range[2]*200)/200, by = 0.005), expand = expansion(mult = c(0.05, 0.05), add = 0))
        } else {
          scale_x_continuous(breaks = seq(floor(min(mz_data$rt)*200)/200, ceiling(max(mz_data$rt)*200)/200, by = 0.005), expand = expansion(mult = c(0.05, 0.05), add = 0))
        }} +
        scale_y_continuous(
          expand = expansion(mult = c(0, 0.02), add = 0),
          limits = c(-y_limit, y_limit),
          labels = function(x) abs(x),
          n.breaks = 6
        ) +
        labs(
          title = mz_label_text,
          x = if (mz_idx == max(mz_indices)) "Retention Time (minutes)" else NULL,
          y = sprintf("\u2190 Std  |  %s \u2192", study_label)
        )
    } else {
      # Single chromatogram plot
      if (stick) {
        p <- ggplot(mz_data, aes(x = rt, y = plot_intensity)) +
          geom_segment(aes(xend = rt, yend = 0), linewidth = 0.8, color = plot_color)
      } else {
        p <- ggplot(mz_data, aes(x = rt, y = plot_intensity)) +
          geom_line(linewidth = 0.8, color = plot_color)
      }
      p <- p +
        {if (manual_rt_override) {
          scale_x_continuous(limits = rt_range, breaks = seq(floor(rt_range[1]*200)/200, ceiling(rt_range[2]*200)/200, by = 0.005), expand = expansion(mult = c(0.05, 0.05), add = 0))
        } else {
          scale_x_continuous(breaks = seq(floor(min(mz_data$rt)*200)/200, ceiling(max(mz_data$rt)*200)/200, by = 0.005), expand = expansion(mult = c(0.05, 0.05), add = 0))
        }} +
        scale_y_continuous(
          expand = expansion(mult = c(0, 0.02), add = 0),
          limits = c(0, y_limit),
          n.breaks = 6
        ) +
        labs(
          title = mz_label_text,
          x = if (mz_idx == max(mz_indices)) "Retention Time (minutes)" else NULL,
          y = "Intensity"
        )
    }
    
    p <- p +
      coord_cartesian(clip = "off") +
      theme_classic(base_size = 10, base_family = "Arial") +
      theme(
        plot.margin = margin(10, 15, 10, 15),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 9, margin = margin(0, 0, 3, 0)),
        axis.text.x = element_text(face = "bold", color = "black", size = 7),
        axis.text.y = element_text(face = "bold", color = "black", size = 7),
        axis.title.x = element_text(face = "bold", color = "black", size = 8),
        axis.title.y = element_text(face = "bold", color = "black", size = 8, margin = margin(r = 8)),
        axis.ticks.length = unit(0.1, "cm"),
        axis.line = element_line(color = "black", linewidth = 0.6),
        axis.ticks = element_line(color = "black", linewidth = 0.6)
      )
    
    # Add library RT marker as red tick on x-axis if available
    if (!is.null(library_trt)) {
      if (standard) {
        p <- p + 
          annotate("segment", x = library_trt, xend = library_trt, 
                   y = -y_limit * 1.05, yend = -y_limit * 0.95, 
                   color = "red", linewidth = 1.2)
      } else {
        p <- p + 
          geom_vline(xintercept = library_trt, linetype = "dashed", color = "red", linewidth = 0.8)
      }
    }
    
    return(p)
  })
  
  # Create title plot (just text, no axes)
  title_plot <- ggplot() +
    theme_void() +
    labs(title = title_text, subtitle = subtitle_text) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11, margin = margin(0, 0, 3, 0)),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 7, color = "black", lineheight = 1.2, margin = margin(0, 0, 5, 0)),
      plot.margin = margin(10, 10, 5, 10)
    )
  
  # Create spectrum plot if requested
  if (spectrum) {
    # Validate rt_unit
    if (!rt_unit %in% c("min", "s")) {
      stop("rt_unit must be either 'min' or 's'")
    }
    
    # Determine RT to use: user-specified or center of rtr
    spectrum_rt <- if (!is.null(rt)) rt else mean(rt_range)
    spectrum_result <- extract_full_spectrum(file_name_sample, spectrum_rt, rt_unit)
    spectrum_data <- spectrum_result$data
    actual_rt_min <- spectrum_result$actual_rt
    
    # Create title with RT info
    if (rt_unit == "s") {
      spectrum_title <- sprintf("Spectrum (%.1f s)", actual_rt_min * 60)
    } else {
      spectrum_title <- sprintf("Spectrum (%.2f min)", actual_rt_min)
    }
    
    # Build spectrum plot
    spectrum_plot <- ggplot(spectrum_data, aes(x = mz, y = intensity)) +
      geom_segment(aes(xend = mz, yend = 0), linewidth = 0.3, color = "black")
    
    # Add peak labels if requested
    if (spectra_labels) {
      # Validate label_tol
      if (label_tol < 0 || label_tol > 1) {
        stop("label_tol must be between 0 and 1")
      }
      
      # Identify top peaks for labeling based on label_tol
      intensity_threshold <- quantile(spectrum_data$intensity, label_tol)
      top_peaks <- spectrum_data |>
        filter(intensity >= intensity_threshold) |>
        arrange(desc(intensity))
      
      spectrum_plot <- spectrum_plot +
        geom_text(data = top_peaks, 
                  aes(label = sprintf("%.1f", mz)), 
                  vjust = -0.3, 
                  size = 1.5, 
                  color = "black",
                  fontface = "bold")
    }
    
    spectrum_plot <- spectrum_plot +
      scale_x_continuous(
        limits = c(85, 850),
        breaks = seq(100, 850, by = 50),
        minor_breaks = seq(90, 850, by = 10),
        expand = expansion(mult = c(0.01, 0.01))
      ) +
      scale_y_continuous(
        expand = expansion(mult = c(0, if (spectra_labels) 0.08 else 0.02)),
        n.breaks = 6
      ) +
      labs(
        title = spectrum_title,
        x = "m/z",
        y = "Intensity"
      ) +
      coord_cartesian(clip = "off") +
      theme_classic(base_size = 10, base_family = "Arial") +
      theme(
        plot.margin = margin(10, 15, 10, 15),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 9, margin = margin(0, 0, 3, 0)),
        axis.text.x = element_text(face = "bold", color = "black", size = 7),
        axis.text.y = element_text(face = "bold", color = "black", size = 7),
        axis.title.x = element_text(face = "bold", color = "black", size = 8),
        axis.title.y = element_text(face = "bold", color = "black", size = 8, margin = margin(r = 8)),
        axis.ticks.length = unit(0.1, "cm"),
        axis.line = element_line(color = "black", linewidth = 0.6),
        axis.ticks = element_line(color = "black", linewidth = 0.6),
        gaxis.ticks.length.x = unit(0.15, "cm")
      )
    
    # Add to plot list
    plot_list <- c(plot_list, list(spectrum_plot))
  }
  
  # Combine plots using cowplot
  combined_plot <- cowplot::plot_grid(
    title_plot,
    cowplot::plot_grid(plotlist = plot_list, ncol = 1, align = "v", axis = "lr"),
    ncol = 1,
    rel_heights = c(0.12, 0.88)
  )
  
  # Save plot to PNG only if png_name is provided
  if (!is.null(png_name)) {
    # Ensure filename has .png extension
    if (!grepl("\\.png$", png_name, ignore.case = TRUE)) {
      output_filename <- paste0(png_name, ".png")
    } else {
      output_filename <- png_name
    }
    
    # Calculate total height based on number of panels
    # Account for spectrum plot being added after mz_indices count
    n_chromatogram_panels <- length(unique(combined_data$mz_index))
    n_spectrum_panels <- if (spectrum) 1 else 0
    total_height <- height * (n_chromatogram_panels + n_spectrum_panels) + 1  # Add 1 inch for title
    
    # Save plot to PNG
    print_to_png(combined_plot, output_filename, width = width, height = total_height, dpi = 300, output_dir = output_dir)
  }
  
  return(invisible(combined_plot))
}
