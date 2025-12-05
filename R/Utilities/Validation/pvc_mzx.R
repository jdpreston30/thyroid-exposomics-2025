#' Plot Validation Mass Spectrum Double (Sample + Standard Mirror)
#'
#' Creates a mirrored mass spectrum plot with sample on top (positive y-axis)
#' and standard below (negative y-axis, flipped) for direct m/z comparison.
#' Extracts spectra within a retention time window and displays as vertical bars.
#'
#' @param id Character string matching the 'id' column in reference_table
#' @param file_name_sample Character string of the sample raw file name (without .mzML extension)
#' @param file_name_standard Character string of the standard raw file name (without .mzML extension)
#' @param mzr Numeric vector c(min, max) for m/z range on x-axis. If NULL (default), calculated from reference table mz0-mz3 values
#' @param rt Numeric value for target retention time in minutes (REQUIRED)
#' @param rt_window Numeric value for +/- window around rt in minutes (default: 0.1)
#' @param ppm_filter Numeric ppm tolerance to filter peaks matching mz0-mz3 from reference table. If NULL (default), shows all peaks
#' @param block_label Character string: "bottom" (default) places clustered labels below standard spectrum, "top" places them above sample spectrum
#' @param source Character string: "IARC" or "quant" to resolve which feature to mark with asterisk if multiple exist (default: NULL uses first)
#' @param reference_table Tibble with columns: id, cas, short_display_name, mz0-mz3 (default: mz_reference_table)
#' @param study Character string: "tumor" or "cadaver" to determine directory and sample color (default: "tumor")
#' @param file_list Tibble with columns: ID, type, files (for sample metadata) (default: file_list)
#' @param selected_features Tibble with columns: source, id, subid, name, tmz, trt, cas (default: selected_gc2_features)
#' @param png_name Optional custom filename for PNG output (without .png extension). If NULL (default), no file is saved
#' @param output_dir Directory to save PNG file (default: "Outputs/Spectra")
#' @param width Plot width in inches (default: 7)
#' @param height Plot height in inches (default: 6)
#'
#' @return ggplot object (saves to PNG only if png_name is provided)
#'
#' @export
pvc_mzx <- function(id,
                    file_name_sample,
                    file_name_standard,
                    mzr = NULL,
                    rt,
                    rt_window = 0.1,
                    ppm_filter = NULL,
                    block_label = "bottom",
                    source = NULL,
                    reference_table = NULL,
                    study = "tumor",
                    file_list = NULL,
                    selected_features = NULL,
                    png_name = NULL,
                    output_dir = "Outputs/Spectra/mzx",
                    width = 3.9,
                    height = 3.25) {
  
  # Check required arguments
  if (missing(rt)) {
    stop("rt argument is required and must be a numeric value (retention time in minutes)")
  }
  
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
  
  # Look up target m/z and retention time from selected_features for asterisk marker
  feature_match <- selected_features |>
    filter(id == !!id)
  
  # If source specified, filter by it (handles duplicates like o-Toluidine)
  if (!is.null(source)) {
    feature_match <- feature_match |>
      filter(source == !!source)
  }
  
  if (nrow(feature_match) == 0) {
    warning(sprintf("ID '%s' not found in selected_features. No target m/z marker will be added.", id))
    target_tmz <- NULL
  } else {
    # Use first match if multiple remain
    target_tmz <- feature_match$tmz[1]
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
    warning(sprintf("No m/z values found for ID '%s'.", id))
    target_mzs <- NULL
  }
  
  # Determine m/z range if not provided
  if (is.null(mzr)) {
    if (!is.null(target_mzs) && length(target_mzs) > 0) {
      # Find min and max m/z values
      min_mz <- min(target_mzs)
      max_mz <- max(target_mzs)
      
      # Calculate span
      span <- max_mz - min_mz
      
      # Choose rounding increment based on span
      # If span < 10, use 10; if span < 50, use 25; otherwise use 100
      increment <- if (span < 10) {
        10
      } else if (span < 50) {
        25
      } else {
        100
      }
      
      # Round down to nearest increment for min, up to nearest increment for max
      mzr_min <- floor(min_mz / increment) * increment
      mzr_max <- ceiling(max_mz / increment) * increment
      
      mzr <- c(mzr_min, mzr_max)
    } else {
      stop("Cannot determine mzr automatically: no m/z values found in reference table. Please provide mzr manually.")
    }
  }
  
  # Determine directory based on study
  raw_dir <- if (study == "tumor") {
    tumor_raw_dir
  } else if (study == "cadaver") {
    cadaver_raw_dir
  } else {
    stop("study must be 'tumor' or 'cadaver'")
  }
  
  mzml_dir <- file.path(raw_dir, "mzML_validation")
  
  # Determine sample color based on study
  sample_color <- if (study == "tumor") {
    "#23744E"  # FV-PTC green
  } else {
    "#BE4E4D"  # Tumor tissue red
  }
  standard_color <- "#294B88"  # Follicular blue
  
  # Helper function to extract mass spectrum data within RT window
  extract_spectrum <- function(file_name) {
    mzml_file <- paste0(file_name, ".mzML")
    mzml_path <- file.path(mzml_dir, mzml_file)
    
    # Check if file exists
    if (!file.exists(mzml_path)) {
      stop(sprintf("mzML file not found: %s", mzml_path))
    }
    
    # Read mzML file
    ms_data <- mzR::openMSfile(mzml_path)
    header_info <- mzR::header(ms_data)
    
    # Find MS1 scans within RT window
    rt_min <- (rt - rt_window) * 60  # Convert to seconds
    rt_max <- (rt + rt_window) * 60
    
    ms1_scans <- header_info[header_info$msLevel == 1, ]
    scans_in_window <- ms1_scans[ms1_scans$retentionTime >= rt_min & 
                                  ms1_scans$retentionTime <= rt_max, ]
    
    if (nrow(scans_in_window) == 0) {
      warning(sprintf("No MS1 scans found in RT window %.2f +/- %.2f for file %s", 
                     rt, rt_window, file_name))
      mzR::close(ms_data)
      return(data.frame(mz = numeric(0), intensity = numeric(0)))
    }
    
    # Extract and aggregate spectra from all scans in window
    all_spectra <- lapply(scans_in_window$seqNum, function(scan_num) {
      spectrum <- mzR::peaks(ms_data, scan_num)
      # Filter to m/z range
      in_range <- spectrum[, 1] >= mzr[1] & spectrum[, 1] <= mzr[2]
      data.frame(
        mz = spectrum[in_range, 1],
        intensity = spectrum[in_range, 2]
      )
    })
    
    # Combine all spectra
    combined_spectrum <- do.call(rbind, all_spectra)
    
    # Average intensities for overlapping m/z values (bin by 0.01 Da)
    combined_spectrum$mz_bin <- round(combined_spectrum$mz, 2)
    spectrum_data <- combined_spectrum |>
      group_by(mz_bin) |>
      summarise(
        mz = mean(mz),
        intensity = mean(intensity),
        .groups = "drop"
      ) |>
      select(mz, intensity)
    
    # Close mzML file
    mzR::close(ms_data)
    
    return(spectrum_data)
  }
  
  # Extract data for both files
  sample_data <- extract_spectrum(file_name_sample)
  sample_data$type <- "Sample"
  sample_data$plot_intensity <- sample_data$intensity  # Keep positive
  sample_data$color_group <- "Sample"
  
  standard_data <- extract_spectrum(file_name_standard)
  standard_data$type <- "Standard"
  standard_data$plot_intensity <- -standard_data$intensity  # Flip to negative for plotting
  standard_data$color_group <- "Standard"
  
  # Combine data
  combined_data <- bind_rows(sample_data, standard_data)
  
  # Filter to only show peaks matching target m/z values if ppm_filter specified
  # Also assign labels (mz0, mz1, etc.) to matching peaks
  mz_labels <- NULL
  if (!is.null(ppm_filter) && !is.null(target_mzs)) {
    # Function to find which target m/z matches (if any) and return label with value
    get_mz_label <- function(mz, targets, ppm) {
      matches <- sapply(seq_along(targets), function(i) {
        if (abs(mz - targets[i]) / targets[i] * 1e6 <= ppm) {
          return(sprintf("mz%d=%.4f", i - 1, targets[i]))
        }
        return(NA_character_)
      })
      matches <- matches[!is.na(matches)]
      if (length(matches) > 0) return(matches[1])
      return(NA_character_)
    }
    
    # Add labels to combined data
    combined_data$mz_label <- sapply(combined_data$mz, get_mz_label, 
                                      targets = target_mzs, ppm = ppm_filter)
    
    # Filter to only matching peaks
    combined_data <- combined_data |>
      filter(!is.na(mz_label))
    
    if (nrow(combined_data) == 0) {
      warning(sprintf("No peaks found within %d ppm of reference m/z values", ppm_filter))
    }
    
    # Create label data based on block_label parameter
    if (block_label == "top") {
      # Place labels at top of sample peaks (positive side)
      mz_labels <- combined_data |>
        filter(type == "Sample") |>
        group_by(mz_label) |>
        slice_max(plot_intensity, n = 1) |>
        ungroup() |>
        select(mz, plot_intensity, mz_label)
    } else {
      # Default: place labels at bottom of standard peaks (negative side)
      mz_labels <- combined_data |>
        filter(type == "Standard") |>
        group_by(mz_label) |>
        slice_min(plot_intensity, n = 1) |>
        ungroup() |>
        select(mz, plot_intensity, mz_label)
    }
    
    # Add asterisk to target m/z AFTER filtering to one spectrum type
    # Extract reference m/z from label string (e.g., "mz0=107.0730" -> 107.0730)
    if (!is.null(target_tmz)) {
      # Extract numeric value from label string
      label_mz_values <- as.numeric(sub("mz[0-9]+=", "", mz_labels$mz_label))
      
      mz_labels$mz_label <- ifelse(
        abs(label_mz_values - target_tmz) < 0.01,
        paste0(mz_labels$mz_label, "*"),
        mz_labels$mz_label
      )
    }
    
    # Check if labels are crowded (any two labels within 5% of x-axis range OR within 10 m/z units)
    # Also identify crowded vs outlier labels
    if (nrow(mz_labels) > 1) {
      mz_range <- diff(range(mz_labels$mz, na.rm = TRUE))
      mz_sorted <- sort(mz_labels$mz)
      min_distance <- min(diff(mz_sorted))
      # Labels are crowded if minimum distance < 5% of range OR < 10 m/z units (absolute)
      labels_crowded <- (min_distance / mz_range) < 0.05 | min_distance < 10
      
      # If crowded, identify which labels are in the cluster vs outliers
      if (labels_crowded) {
        # Use hierarchical clustering to group close peaks
        # Define closeness as within 10 m/z units
        mz_sorted_df <- data.frame(mz = sort(mz_labels$mz))
        mz_sorted_df$cluster_id <- 1
        
        if (nrow(mz_sorted_df) > 1) {
          cluster_counter <- 1
          for (i in 2:nrow(mz_sorted_df)) {
            # If this peak is within 10 units of previous peak, same cluster
            if (mz_sorted_df$mz[i] - mz_sorted_df$mz[i-1] <= 10) {
              mz_sorted_df$cluster_id[i] <- cluster_counter
            } else {
              # Start new cluster
              cluster_counter <- cluster_counter + 1
              mz_sorted_df$cluster_id[i] <- cluster_counter
            }
          }
        }
        
        # Identify largest cluster as the "crowded" group
        cluster_sizes <- table(mz_sorted_df$cluster_id)
        largest_cluster <- as.numeric(names(cluster_sizes)[which.max(cluster_sizes)])
        
        mz_labels <- mz_labels |>
          left_join(mz_sorted_df, by = "mz") |>
          mutate(is_clustered = (cluster_id == largest_cluster) & (max(cluster_sizes) > 1))
        
        # Add nudge directions for clustered labels (angle them)
        if (sum(mz_labels$is_clustered) > 1) {
          clustered_labels <- mz_labels |> filter(is_clustered) |> arrange(mz)
          n_clustered <- nrow(clustered_labels)
          
          # Assign angles: left 90°, center 0°, right 90° (or distribute evenly)
          if (n_clustered == 2) {
            angles <- c(-90, 90)
          } else if (n_clustered == 3) {
            angles <- c(-90, 0, 90)
          } else {
            # For more than 3, distribute evenly
            angles <- seq(-90, 90, length.out = n_clustered)
          }
          
          # Create a mapping of mz to angle based on sorted order
          clustered_labels$angle <- angles
          
          mz_labels <- mz_labels |>
            left_join(clustered_labels |> select(mz, angle), by = "mz") |>
            mutate(angle = if_else(is.na(angle), 0, angle))
        } else {
          mz_labels$angle <- 0
        }
      } else {
        mz_labels$is_clustered <- FALSE
        mz_labels$angle <- 0
      }
    } else {
      labels_crowded <- FALSE
      if (nrow(mz_labels) > 0) {
        mz_labels$is_clustered <- FALSE
        mz_labels$angle <- 0
      }
    }
  } else {
    labels_crowded <- FALSE
  }
  
  # Determine y-axis limits based on max intensity from either side
  max_sample <- max(abs(combined_data$intensity[combined_data$type == "Sample"]), na.rm = TRUE)
  max_standard <- max(abs(combined_data$intensity[combined_data$type == "Standard"]), na.rm = TRUE)
  y_limit <- max(max_sample, max_standard)
  
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
  title_text <- sprintf("%s vs. Standard", study_label)
  
  # Create subtitle with metadata
  subtitle_text <- sprintf("Sample: %s (%s)  |  Standard: %s\n%s  |  ID: %s  |  RT: %.2f ± %.2f min", 
                          sample_id, file_name_sample, file_name_standard, ref_row$short_display_name, id, rt, rt_window)
  
  # Create color scale
  color_values <- c("Sample" = sample_color, "Standard" = standard_color)
  
  # Create plot with vertical bars (geom_segment for spectral lines)
  p <- ggplot(combined_data, aes(x = mz, y = plot_intensity, color = color_group)) +
    geom_segment(aes(xend = mz, yend = 0), linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    scale_color_manual(values = color_values, guide = "none") +
    scale_x_continuous(limits = mzr, expand = expansion(mult = c(0.05, 0.05), add = 0)) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.05), add = 0),
      limits = c(-y_limit, y_limit),
      labels = function(x) abs(x),
      n.breaks = 8
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "m/z",
      y = sprintf("\u2190 Standard  |  %s \u2192", study_label)
    ) +
    # Add m/z labels above sample peaks if ppm_filter was used
    {if (!is.null(mz_labels) && nrow(mz_labels) > 0 && labels_crowded) {
      # Create combined label for clustered peaks and separate labels for outliers
      clustered_data <- mz_labels |> filter(is_clustered)
      outlier_data <- mz_labels |> filter(!is_clustered)
      
      if (nrow(clustered_data) > 0) {
        # Sort by mz and create multi-line label
        clustered_sorted <- clustered_data |> arrange(mz)
        combined_label <- paste(clustered_sorted$mz_label, collapse = "\n")
        
        # Use middle peak position and appropriate intensity based on block_label
        middle_idx <- ceiling(nrow(clustered_sorted) / 2)
        if (block_label == "top") {
          # Place at top (max intensity)
          label_df <- data.frame(
            mz = clustered_sorted$mz[middle_idx],
            plot_intensity = max(clustered_sorted$plot_intensity),
            label = combined_label
          )
        } else {
          # Place at bottom (min intensity)
          label_df <- data.frame(
            mz = clustered_sorted$mz[middle_idx],
            plot_intensity = min(clustered_sorted$plot_intensity),
            label = combined_label
          )
        }
      } else {
        label_df <- data.frame(mz = numeric(0), plot_intensity = numeric(0), label = character(0))
      }
      
      list(
        # Multi-line label for clustered peaks
        if (nrow(label_df) > 0) {
          geom_text(data = label_df, 
                    aes(x = mz, y = plot_intensity, label = label),
                    color = "black", size = 2, fontface = "italic", 
                    vjust = if (block_label == "top") -0.5 else 1.5, 
                    hjust = 0.5, lineheight = 0.9,
                    inherit.aes = FALSE)
        },
        # Regular labels for outlier peaks
        if (nrow(outlier_data) > 0) {
          geom_text(data = outlier_data, 
                    aes(x = mz, y = plot_intensity, label = mz_label),
                    color = "black", size = 2, fontface = "italic", 
                    vjust = if (block_label == "top") -0.5 else 1.5, 
                    hjust = 0.5, inherit.aes = FALSE)
        }
      )
    } else if (!is.null(mz_labels) && nrow(mz_labels) > 0) {
      geom_text(data = mz_labels, 
                aes(x = mz, y = plot_intensity, label = mz_label),
                color = "black", size = 2, fontface = "italic", 
                vjust = if (block_label == "top") -0.5 else 1.5, 
                hjust = 0.5, inherit.aes = FALSE)
    }} +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      plot.margin = margin(20, 20, 20, 20),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10, margin = margin(0, 0, 3, 0)),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 6, color = "black", lineheight = 1.2, margin = margin(0, 0, 5, 0)),
      axis.text.x = element_text(face = "bold", color = "black", size = 8),
      axis.text.y = element_text(face = "bold", color = "black", size = 8),
      axis.title.x = element_text(face = "bold", color = "black", size = 10),
      axis.title.y = element_text(face = "bold", color = "black", size = 10, margin = margin(r = 10)),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8)
    )
  
  # Save plot to PNG only if png_name is provided
  if (!is.null(png_name)) {
    # Ensure filename has .png extension
    if (!grepl("\\.png$", png_name, ignore.case = TRUE)) {
      output_filename <- paste0(png_name, ".png")
    } else {
      output_filename <- png_name
    }
    
    # Save plot to PNG
    print_to_png(p, output_filename, width = width, height = height, dpi = 300, output_dir = output_dir)
  }
  
  return(invisible(p))
}
