#' Plot Validation Chromatogram Double (Sample + Standard Mirror)
#'
#' Creates a mirrored chromatogram plot with sample on top (positive y-axis)
#' and standard below (negative y-axis, flipped) for direct RT comparison.
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
#' @param png_name Optional custom filename for PNG output (without .png extension). If NULL, uses default naming: id_sample_replicate_std_replicate.png
#' @param output_dir Directory to save PNG file (default: "Outputs/Spectra")
#' @param width Plot width in inches (default: 6)
#' @param height Plot height in inches (default: 6)
#'
#' @return ggplot object (also saves to PNG)
#'
#' @export
pvcd <- function(id,
                 file_name_sample,
                 file_name_standard,
                 reference_table = NULL,
                 study = "tumor",
                 file_list = NULL,
                 selected_features = NULL,
                 source = NULL,
                 ppm_tolerance = 5,
                 rtr = NULL,
                 png_name = NULL,
                 output_dir = "Outputs/Spectra",
                 width = 7,
                 height = 6) {
  
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
  
  # Determine RT range - use manual if provided, otherwise dynamic
  # Track whether manual override was used for later axis control
  manual_rt_override <- !is.null(rtr)
  
  if (!is.null(rtr)) {
    rt_range <- rtr
  } else {
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
    extract_xic <- function(ms_data, header_info, target_mz, ppm_tol, rt_range) {
      ms1_scans <- header_info[header_info$msLevel == 1, ]
      ms1_scans <- ms1_scans[ms1_scans$retentionTime >= rt_range[1] * 60 & 
                             ms1_scans$retentionTime <= rt_range[2] * 60, ]
      
      # Calculate m/z tolerance in Daltons from ppm
      # ppm = (delta_mz / target_mz) * 1e6
      # So: delta_mz = target_mz * ppm / 1e6
      mz_tol_da <- target_mz * ppm_tol / 1e6
      
      xic_data <- lapply(ms1_scans$seqNum, function(scan_num) {
        spectrum <- mzR::peaks(ms_data, scan_num)
        mz_match <- abs(spectrum[, 1] - target_mz) <= mz_tol_da
        intensity <- if (sum(mz_match) > 0) sum(spectrum[mz_match, 2]) else 0
        
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
      xic <- extract_xic(ms_data, header_info, target_mzs[i], ppm_tolerance, rt_range)
      xic$mz_index <- i - 1
      xic
    })
    
    chromatogram_data <- do.call(rbind, all_xics)
    
    # Close mzML file
    mzR::close(ms_data)
    
    return(chromatogram_data)
  }
  
  # Extract data for both files
  sample_data <- extract_chromatogram(file_name_sample)
  sample_data$type <- "Sample"
  sample_data$plot_intensity <- sample_data$intensity  # Keep positive
  
  standard_data <- extract_chromatogram(file_name_standard)
  standard_data$type <- "Standard"
  standard_data$plot_intensity <- -standard_data$intensity  # Flip to negative for plotting
  
  # Determine y-axis limits based on max intensity from either side
  max_sample <- max(abs(sample_data$intensity), na.rm = TRUE)
  max_standard <- max(abs(standard_data$intensity), na.rm = TRUE)
  y_limit <- max(max_sample, max_standard)
  
  # Combine data
  combined_data <- bind_rows(sample_data, standard_data)
  
  # Create labels with asterisk for target m/z
  combined_data$mz_label <- sprintf("mz%d: %.2f", 
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
  
  # Get standard ID
  standard_id <- gsub("_[0-9]+$", "", file_name_standard)
  
  # Create dynamic title based on study type
  study_label <- tools::toTitleCase(study)
  title_text <- sprintf("%s vs Standard", study_label)
  
  # Create subtitle with metadata (add library RT if available)
  if (!is.null(library_trt)) {
    subtitle_text <- sprintf("Sample: %s (%s)  |  Standard: %s  |  ID: %s  |  %d ppm  |  Lib RT = %.2f", 
                            sample_id, file_name_sample, file_name_standard, id, ppm_tolerance, library_trt)
  } else {
    subtitle_text <- sprintf("Sample: %s (%s)  |  Standard: %s  |  ID: %s  |  %d ppm", 
                            sample_id, file_name_sample, file_name_standard, id, ppm_tolerance)
  }
  
  # Create plot
  p <- ggplot(combined_data, aes(x = rt, y = plot_intensity, color = mz_label, group = interaction(mz_label, type))) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    scale_color_brewer(palette = "Set1") +
    {if (manual_rt_override) {
      scale_x_continuous(limits = rt_range, expand = expansion(mult = c(0.01, 0.01), add = 0))
    } else {
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.05), add = 0))
    }} +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.05), add = 0),
      limits = c(-y_limit, y_limit),
      labels = function(x) abs(x),
      n.breaks = 8
    ) +
    labs(
      title = sprintf("%s (%s)", ref_row$short_display_name, title_text),
      subtitle = subtitle_text,
      x = "Retention Time (minutes)",
      y = sprintf("\u2190 Standard Intensity              %s Intensity \u2192", study_label),
      color = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      plot.margin = margin(20, 20, 20, 20),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.position = "top",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.text = element_text(size = 7.5, face = "plain", family = "Arial"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.key.size = unit(0.6, "cm"),
      legend.spacing.x = unit(0.3, "cm"),
      legend.box.margin = margin(-5, 0, -10, 0),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 8, color = "black", lineheight = 1.2, margin = margin(0, 0, 0, 0)),
      axis.text.x = element_text(face = "bold", color = "black", size = 10),
      axis.text.y = element_text(face = "bold", color = "black", size = 10),
      axis.title.x = element_text(face = "bold", color = "black", size = 12),
      axis.title.y = element_text(face = "bold", color = "black", size = 12, margin = margin(r = 10)),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8)
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 0.8)))
  
  # Add library RT marker as red tick on x-axis if available
  if (!is.null(library_trt)) {
    p <- p + 
      annotate("segment", x = library_trt, xend = library_trt, 
               y = -y_limit * 1.05, yend = -y_limit * 0.95, 
               color = "red", linewidth = 1.2)
  }
  
  # Determine output filename
  if (!is.null(png_name)) {
    # Use custom name provided by user
    output_filename <- paste0(png_name, ".png")
  } else {
    # Use default naming convention
    sample_replicate <- sub(".*_([0-9]+)$", "\\1", file_name_sample)
    standard_replicate <- sub(".*_([0-9]+)$", "\\1", file_name_standard)
    output_filename <- sprintf("%s_%s_%s_std_%s.png", id, sample_id, sample_replicate, standard_replicate)
  }
  
  # Save plot to PNG
  print_to_png(p, output_filename, width = width, height = height, output_dir = output_dir)
  
  return(invisible(p))
}
