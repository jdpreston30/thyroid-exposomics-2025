#' Plot Validation Chromatogram for Manual Spectral Review
#'
#' Creates an extracted ion chromatogram (XIC) plot for a specific sample and feature
#' using target m/z values from the reference table. Adds sample metadata as text annotation.
#'
#' @param id Character string matching the 'id' column in reference_table
#' @param file_name Character string of the raw file name (without .mzML extension)
#' @param reference_table Tibble with columns: id, cas, mz0-mz3, tumor_rt_range, cadaver_rt_range (default: mz_reference_table)
#' @param study Character string: "tumor" or "cadaver" to determine directory and RT range (default: "tumor")
#' @param file_list Tibble with columns: ID, type, files (for sample metadata) (default: file_list)
#' @param mz_tolerance Numeric m/z tolerance in Daltons (default: 0.5)
#' @param rtr Optional numeric vector c(min, max) to override dynamic RT range
#' @param output_dir Directory to save PNG file (default: "Outputs/Spectra")
#' @param width Plot width in inches (default: 6)
#' @param height Plot height in inches (default: 6)
#'
#' @return ggplot object (also saves to PNG)
#'
#' @export
pvc <- function(id,
                file_name,
                reference_table = NULL,
                study = "tumor",
                file_list = NULL,
                mz_tolerance = 0.5,
                rtr = NULL,
                output_dir = "Outputs/Spectra",
                width = 6,
                height = 6) {
  
  # Use defaults from environment if not provided
  if (is.null(reference_table)) {
    reference_table <- mz_reference_table
  }
  if (is.null(file_list)) {
    file_list <- get("file_list", envir = parent.frame())
  }
  
  # Get directory paths from config
  tumor_raw_dir <- config$paths$tumor_raw_dir
  cadaver_raw_dir <- config$paths$cadaver_raw_dir
  
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
  mzml_file <- paste0(file_name, ".mzML")
  mzml_path <- file.path(mzml_dir, mzml_file)
  
  # Check if file exists
  if (!file.exists(mzml_path)) {
    stop(sprintf("mzML file not found: %s", mzml_path))
  }
  
  # Get sample metadata from file_list by matching file_name in files column
  sample_metadata <- file_list |>
    filter(grepl(file_name, files, fixed = TRUE))
  
  if (nrow(sample_metadata) == 0) {
    sample_id <- gsub("_[0-9]+$", "", file_name)
    sample_type <- "Unknown"
  } else {
    sample_id <- sample_metadata$ID[1]
    sample_type <- sample_metadata$type[1]
    
    # Capitalize type for display
    sample_type <- tools::toTitleCase(sample_type)
  }
  
  # Read mzML file
  ms_data <- mzR::openMSfile(mzml_path)
  header_info <- mzR::header(ms_data)
  
  # Function to extract XIC (extracted ion chromatogram)
  extract_xic <- function(ms_data, header_info, target_mz, mz_tol, rt_range) {
    
    # Filter MS1 scans only
    ms1_scans <- header_info[header_info$msLevel == 1, ]
    
    # Apply RT filter
    ms1_scans <- ms1_scans[ms1_scans$retentionTime >= rt_range[1] * 60 & 
                           ms1_scans$retentionTime <= rt_range[2] * 60, ]
    
    # Extract intensity for target m/z across all scans
    xic_data <- lapply(ms1_scans$seqNum, function(scan_num) {
      spectrum <- mzR::peaks(ms_data, scan_num)
      mz_match <- abs(spectrum[, 1] - target_mz) <= mz_tol
      
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
    xic <- extract_xic(ms_data, header_info, target_mzs[i], mz_tolerance, rt_range)
    xic$mz_index <- i - 1  # 0-indexed to match mz0, mz1, etc.
    xic
  })
  
  chromatogram_data <- do.call(rbind, all_xics)
  chromatogram_data$mz_label <- sprintf("mz%d: %.2f", 
                                        chromatogram_data$mz_index, 
                                        chromatogram_data$mz)
  
  # Close mzML file
  mzR::close(ms_data)
  
  # Create subtitle with metadata
  subtitle_text <- sprintf("ID: %s  |  CAS: %s  |  Sample: %s  |  File: %s", 
                          id, toupper(ref_row$cas), sample_id, file_name)
  
  # Create plot
  p <- ggplot(chromatogram_data, aes(x = rt, y = intensity, color = mz_label)) +
    geom_line(linewidth = 0.8) +
    scale_color_brewer(palette = "Set1") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05), add = 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05), add = 0)) +
    labs(
      title = sprintf("%s (%s)", ref_row$short_display_name, sample_type),
      subtitle = subtitle_text,
      x = "Retention Time (minutes)",
      y = "Intensity",
      color = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      plot.margin = margin(20, 20, 20, 20),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.text = element_text(size = 7.5, face = "plain", family = "Arial"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.key.size = unit(0.6, "cm"),
      legend.spacing.y = unit(0.2, "cm"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 8, color = "black", lineheight = 1.2),
      axis.text.x = element_text(face = "bold", color = "black", size = 10),
      axis.text.y = element_text(face = "bold", color = "black", size = 10),
      axis.title.x = element_text(face = "bold", color = "black", size = 12),
      axis.title.y = element_text(face = "bold", color = "black", size = 12, margin = margin(r = 10)),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8)
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 0.8)))
  
  # Extract replicate number from file_name (last digit before extension)
  replicate <- sub(".*_([0-9]+)$", "\\1", file_name)
  
  # Create output filename: id_sample_replicate.png
  output_filename <- sprintf("%s_%s_%s.png", id, sample_id, replicate)
  
  # Save plot to PNG
  print_to_png(p, output_filename, width = width, height = height, output_dir = output_dir)
  
  return(invisible(p))
}
