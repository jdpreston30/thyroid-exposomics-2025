# Required libraries: mzR, ggplot2, dplyr, gridExtra, ggtext
library(mzR)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(ggtext)

#' Process Single Compound (Worker Function for Parallel Processing)
#' @keywords internal
process_single_compound <- function(row, row_idx, total_rows, mzml_dir, iterate_through,
                                   rt_lookup, window, ppm_tolerance, stick, max_i,
                                   save_rds, rds_save_folder, overwrite_rds, output_dir, study = "tumor",
                                   run_standard = TRUE, fragment_pare = FALSE, force_plot = FALSE, debug = FALSE) {
  
  # Helper function to open and cache mzML files (per-worker cache)
  get_mzml_data_worker <- function(file_name, cache_env) {
    if (exists(file_name, envir = cache_env)) {
      return(cache_env[[file_name]])
    }
    
    mzml_file <- paste0(file_name, ".mzML")
    mzml_path <- file.path(mzml_dir, mzml_file)
    
    if (!file.exists(mzml_path)) {
      return(NULL)
    }
    
    ms_data <- mzR::openMSfile(mzml_path)
    header_info <- mzR::header(ms_data)
    
    cache_env[[file_name]] <- list(
      ms_data = ms_data,
      header_info = header_info
    )
    
    return(cache_env[[file_name]])
  }
  
  # Helper function to extract chromatogram
  extract_chrom_worker <- function(file_name, target_mzs, rt_range, ppm_tol, use_max, cache_env, debug = FALSE) {
    mzml_data <- get_mzml_data_worker(file_name, cache_env)
    if (is.null(mzml_data)) {
      if (debug) cat(sprintf("  [DEBUG] File %s could not be opened\n", file_name))
      return(NULL)
    }
    
    ms_data <- mzml_data$ms_data
    header_info <- mzml_data$header_info
    
    extract_xic <- function(ms_data, header_info, target_mz, ppm_tol, rt_range, use_max) {
      ms1_scans <- header_info[header_info$msLevel == 1, ]
      total_scans <- nrow(ms1_scans)
      ms1_scans <- ms1_scans[ms1_scans$retentionTime >= rt_range[1] * 60 & 
                             ms1_scans$retentionTime <= rt_range[2] * 60, ]
      scans_in_range <- nrow(ms1_scans)
      
      if (debug && scans_in_range == 0) {
        cat(sprintf("  [DEBUG] No scans found in RT range %.2f-%.2f min (total scans: %d, RT range: %.2f-%.2f min)\n", 
                    rt_range[1], rt_range[2], total_scans, 
                    min(header_info$retentionTime)/60, max(header_info$retentionTime)/60))
      }
      
      mz_tol_da <- target_mz * ppm_tol / 1e6
      
      xic_data <- lapply(ms1_scans$seqNum, function(scan_num) {
        spectrum <- mzR::peaks(ms_data, scan_num)
        mz_match <- abs(spectrum[, 1] - target_mz) <= mz_tol_da
        
        if (sum(mz_match) > 0) {
          intensity <- if (use_max) max(spectrum[mz_match, 2]) else sum(spectrum[mz_match, 2])
          rt <- ms1_scans$retentionTime[ms1_scans$seqNum == scan_num] / 60
          data.frame(rt = rt, intensity = intensity, mz = target_mz)
        } else {
          NULL
        }
      })
      
      xic_result <- do.call(rbind, xic_data)
      
      if (debug && (is.null(xic_result) || nrow(xic_result) == 0)) {
        cat(sprintf("  [DEBUG] No ions detected at m/z %.4f (±%.4f Da) in %d scans\n", 
                    target_mz, mz_tol_da, scans_in_range))
      }
      
      xic_result
    }
    
    all_xics <- lapply(seq_along(target_mzs), function(i) {
      mz_val <- target_mzs[i]
      xic <- extract_xic(ms_data, header_info, mz_val, ppm_tol, rt_range, use_max)
      if (!is.null(xic)) {
        xic$mz_index <- i - 1
        xic
      }
    })
    
    result <- do.call(rbind, all_xics)
    
    if (debug && is.null(result)) {
      cat(sprintf("  [DEBUG] extract_chrom_worker returning NULL for file %s (target m/z: %s)\n", 
                  file_name, paste(round(target_mzs, 4), collapse = ", ")))
    }
    
    result
  }
  
  # Initialize worker cache
  worker_cache <- new.env(hash = TRUE)
  
  # Extract compound info
  order_num <- row$order
  id_val <- row$id
  short_name <- row$short_name
  
  # Sanitize short_name: replace Greek letters with text equivalents
  short_name <- gsub("\u03B1", "alpha", short_name)  # α
  short_name <- gsub("\u03B2", "beta", short_name)   # β
  short_name <- gsub("\u03B3", "gamma", short_name)  # γ
  short_name <- gsub("\u03B4", "delta", short_name)  # δ
  short_name <- gsub("\u03BC", "mu", short_name)     # μ
  
  compound_result <- list(
    short_name = short_name,
    id = id_val,
    order = order_num,
    plots = list()
  )
  
  # Extract m/z values
  original_frag_index <- NULL
  if (fragment_pare && "top_frag" %in% names(row) && !is.na(row$top_frag)) {
    # Fragment paring mode: only extract the specific fragment from top_frag column
    original_frag_index <- as.numeric(row$top_frag)
    top_frag_col <- paste0("mz", row$top_frag)
    if (top_frag_col %in% names(row)) {
      target_mzs <- as.numeric(row[[top_frag_col]])
      if (is.na(target_mzs)) target_mzs <- numeric(0)
    } else {
      target_mzs <- numeric(0)
    }
  } else {
    # Default mode: extract all mz values
    target_mzs <- row |>
      dplyr::select(matches("^mz[0-9]+$")) |>
      unlist() |>
      as.numeric() |>
      na.omit()
  }
  
  if (length(target_mzs) == 0) return(compound_result)
  
  # Get RT range
  rt_range_str <- row$compound_rt_range
  if (is.na(rt_range_str)) return(compound_result)
  
  base_rt_range <- eval(parse(text = rt_range_str))
  base_rt_range_expanded <- c(base_rt_range[1] - 0.2, base_rt_range[2] + 0.2)
  
  # Parse standards
  standards <- strsplit(row$standards, ", ")[[1]]
  
  # Get sample files - dynamically collect all file columns that exist
  file_cols <- grep("^file\\d+$", names(row), value = TRUE)
  all_samples <- as.character(row[file_cols])
  names(all_samples) <- NULL
  
  if (force_plot) {
    # When force_plot=TRUE, iterate through ALL slots up to iterate_through
    samples_to_process <- all_samples[1:min(iterate_through, length(all_samples))]
  } else {
    # Normal behavior: only process non-NA files
    all_samples <- all_samples[!is.na(all_samples)]
    samples_to_process <- all_samples[1:min(iterate_through, length(all_samples))]
  }
  
  # Process each sample
  for (sample_idx in seq_along(samples_to_process)) {
    sample_file <- samples_to_process[sample_idx]
    
    # Skip if file is NA and force_plot is FALSE
    if (is.na(sample_file) && !force_plot) next
    
    # Determine RT range
    sample_rt_range <- base_rt_range_expanded
    use_hard_limits <- FALSE
    rt_is_fallback <- FALSE
    
    if (rt_lookup == "sample") {
      rt_range_col_name <- paste0("f", sample_idx, "_rt_range")
      if (rt_range_col_name %in% names(row) && !is.na(row[[rt_range_col_name]])) {
        rt_range_str <- row[[rt_range_col_name]]
        sample_rt_range <- eval(parse(text = rt_range_str))
        sample_rt_range <- c(sample_rt_range[1] - 0.2, sample_rt_range[2] + 0.2)
      }
    } else if (rt_lookup == "window") {
      rt_col_name <- paste0("f", sample_idx, "_rt")
      if (rt_col_name %in% names(row) && !is.na(row[[rt_col_name]])) {
        rt_value <- row[[rt_col_name]]
        half_window <- window / 2
        sample_rt_range <- c(rt_value - half_window, rt_value + half_window)
        use_hard_limits <- TRUE
      } else {
        rt_is_fallback <- TRUE
      }
    }
    
    # Process based on run_standard flag
    if (!run_standard) {
      #+ Sample-only mode (no standard comparison)
      # Extract sample chromatogram
      if (debug) {
        cat(sprintf("\n[DEBUG] Processing %s (%s), sample %d: %s\n", short_name, id_val, sample_idx, sample_file))
        cat(sprintf("[DEBUG] RT range: %.2f-%.2f, target m/z: %s\n", 
                    sample_rt_range[1], sample_rt_range[2], paste(round(target_mzs, 4), collapse = ", ")))
      }
      sample_chrom <- extract_chrom_worker(sample_file, target_mzs, sample_rt_range, ppm_tolerance, max_i, worker_cache, debug)
      
      if (is.null(sample_chrom) && !force_plot) next
      
      # Track whether we have actual data
      has_data <- !is.null(sample_chrom)
      if (debug && !has_data) {
        cat(sprintf("[DEBUG] No data extracted for %s - creating placeholder plot\n", short_name))
      }
      
      # Only process sample_chrom fields if we have data
      if (has_data) {
        # Adjust mz_index if fragment_pare is TRUE to preserve original fragment number
        if (!is.null(original_frag_index)) {
          sample_chrom$mz_index <- original_frag_index
        }
        
        sample_chrom$type <- "Sample"
        
        # Create mz labels with actual m/z values (before stick transformation)
        sample_chrom$mz_label <- sprintf("mz%d: %.4f", sample_chrom$mz_index, sample_chrom$mz)
      }
      
      # Apply stick transformation if needed (only if we have data)
      if (has_data && stick) {
        sample_chrom <- sample_chrom |>
          group_by(rt, mz_label) |>
          summarize(intensity = max(intensity), .groups = "drop") |>
          group_by(mz_label) |>
          arrange(rt) |>
          mutate(
            rt_start = rt,
            rt_end = rt,
            intensity_start = 0,
            intensity_end = intensity
          ) |>
          ungroup()
      }
      
      # Filter to only rows with non-zero intensity (only if we have data)
      if (has_data) {
        sample_chrom <- sample_chrom |>
          filter(intensity > 0)
      }
      
      if ((is.null(sample_chrom) || nrow(sample_chrom) == 0) && !force_plot) next
      
      # Recheck if we still have data after filtering
      has_data <- !is.null(sample_chrom) && nrow(sample_chrom) > 0
      
      if (has_data) {
        # Calculate intensity limits
        y_limit <- max(abs(sample_chrom$intensity), na.rm = TRUE) * 1.05
      } else {
        # No data - use default y limit for empty plot
        y_limit <- 1000
      }
      
      # Add asterisks to marked m/z values (only if we have data)
      if (has_data && !is.na(row$asterisk)) {
        marked_mzs <- strsplit(row$asterisk, ", ")[[1]]
        for (marked_mz in marked_mzs) {
          mz_num <- as.numeric(gsub("mz", "", marked_mz))
          sample_chrom$mz_label <- ifelse(
            sample_chrom$mz_index == mz_num,
            paste0("**", sample_chrom$mz_label, " \\*****"),
            sample_chrom$mz_label
          )
        }
      }
      
      sample_id <- if (!is.na(sample_file)) gsub("\\.mzML$", "", basename(sample_file)) else sprintf("File%d", sample_idx)
      
      if (has_data) {
        # Define fixed color palette for mz indices (matches vp() function palette)
        mz_colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#D4A017", "#0072B2", "#D55E00", "#CC79A7", "#999999")
        names(mz_colors) <- paste0("mz", 0:8)
        
        # Create named vector for actual mz_labels present in the data
        mz_labels_present <- unique(sample_chrom$mz_label)
        mz_indices_present <- unique(sample_chrom$mz_index)
        color_mapping <- setNames(
          mz_colors[paste0("mz", mz_indices_present)],
          mz_labels_present
        )
        
        if (stick) {
          p_rtx <- ggplot(sample_chrom, aes(x = rt, y = intensity, color = mz_label)) +
            geom_segment(aes(xend = rt, yend = 0), linewidth = 0.4)
        } else {
          p_rtx <- ggplot(sample_chrom, aes(x = rt, y = intensity, color = mz_label)) +
            geom_line(linewidth = 0.4)
        }
        
        p_rtx <- p_rtx +
          scale_color_manual(values = color_mapping) +
          scale_y_continuous(
            expand = c(0, 0),
            limits = c(0, y_limit),
            n.breaks = 8,
            labels = scales::label_scientific(digits = 2)
          )
      } else {
        # Create empty plot for no-data case
        p_rtx <- ggplot() +
          annotate("text", x = mean(sample_rt_range), y = y_limit/2, 
                  label = "NO DATA", size = 6, color = "gray50", fontface = "bold") +
          scale_y_continuous(
            expand = c(0, 0),
            limits = c(0, y_limit),
            n.breaks = 8,
            labels = scales::label_scientific(digits = 2)
          )
      }
      
      # Apply x-axis scaling based on whether hard limits are used
      if (use_hard_limits) {
        p_rtx <- p_rtx +
          scale_x_continuous(
            expand = c(0, 0),
            limits = sample_rt_range,
            breaks = function(limits) {
              start <- ceiling(limits[1] * 20) / 20
              end <- floor(limits[2] * 20) / 20
              if (start < end) seq(start, end, by = 0.05) else seq(start, end, by = -0.05)
            },
            minor_breaks = function(limits) {
              start <- ceiling(limits[1] * 40) / 40
              end <- floor(limits[2] * 40) / 40
              if (start < end) seq(start, end, by = 0.025) else seq(start, end, by = -0.025)
            }
          )
      } else {
        p_rtx <- p_rtx +
          scale_x_continuous(
            expand = expansion(mult = c(0.05, 0.05), add = 0),
            breaks = function(limits) {
              start <- ceiling(limits[1] * 20) / 20
              end <- floor(limits[2] * 20) / 20
              if (start < end) seq(start, end, by = 0.05) else seq(start, end, by = -0.05)
            },
            minor_breaks = function(limits) {
              start <- ceiling(limits[1] * 40) / 40
              end <- floor(limits[2] * 40) / 40
              if (start < end) seq(start, end, by = 0.025) else seq(start, end, by = -0.025)
            }
          )
      }
      
      # Create subtitle with RT info
      if (rt_is_fallback) {
        subtitle_text <- sprintf("Sample: %s  |  RT = NA (range: %.2f-%.2f min)", 
                                sample_id, sample_rt_range[1], sample_rt_range[2])
      } else {
        subtitle_text <- sprintf("Sample: %s  |  RT = %.2f min", sample_id, mean(sample_rt_range))
      }
      
      p_rtx <- p_rtx +
        labs(
          title = short_name,
          subtitle = subtitle_text,
          x = "Retention Time (min)",
          y = "Intensity",
          color = NULL
        ) +
        coord_cartesian(clip = "off") +
        theme_classic(base_size = 12) +
        theme(
          panel.grid.major.x = element_line(color = "gray90", linewidth = 0.2),
          panel.grid.minor.x = element_line(color = "gray90", linewidth = 0.2),
          plot.margin = margin(10, 10, 20, 20),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "top",
          legend.justification = "center",
          legend.direction = "horizontal",
          legend.text = ggtext::element_markdown(size = 4),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA),
          legend.key.size = unit(0.25, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.spacing.x = unit(0.02, "cm"),
          legend.box.margin = margin(0, 0, 0, 0),
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
      
      # Store plot (use S1 as standard index for sample-only)
      plot_label <- sprintf("F%d_S1", sample_idx)
      plot_tag <- sprintf("F%d_S1_%s", sample_idx, id_val)
      if (study == "cadaver") {
        plot_tag <- paste0("C_", plot_tag)
      }
      
      compound_result$plots[[plot_label]] <- list(
        plot = p_rtx,
        sample_id = sample_id,
        standard_file = NA,
        plot_tag = plot_tag,
        rt_range = sample_rt_range
      )
      
      # Save RDS if requested (sample-only mode)
      if (save_rds && !is.null(rds_save_folder)) {
        if (exists("config") && !is.null(config$paths$validation_plot_directory)) {
          rds_dir <- file.path(config$paths$validation_plot_directory, rds_save_folder)
        } else {
          rds_dir <- file.path(output_dir, "RDS", rds_save_folder)
        }
        
        dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
        rds_path <- file.path(rds_dir, paste0(plot_tag, ".rds"))
        
        if (!file.exists(rds_path) || overwrite_rds) {
          individual_plot <- list(
            short_name = short_name,
            id = id_val,
            order = order_num,
            plot = p_rtx,
            sample_id = sample_id,
            standard_file = NA,
            plot_tag = plot_tag,
            rt_range = sample_rt_range
          )
          saveRDS(individual_plot, file = rds_path, compress = "gzip")
        }
      }
      
    } else {
      #+ Standard comparison mode
      standards_to_process <- standards
      for (std_idx in seq_along(standards_to_process)) {
        standard_file <- standards_to_process[std_idx]
        
        # Extract chromatograms
        sample_chrom <- extract_chrom_worker(sample_file, target_mzs, sample_rt_range, ppm_tolerance, max_i, worker_cache)
        std_chrom <- extract_chrom_worker(standard_file, target_mzs, sample_rt_range, ppm_tolerance, max_i, worker_cache)
        
        if ((is.null(sample_chrom) || is.null(std_chrom)) && !force_plot) next
        
        # Track whether we have actual data
        has_data <- !is.null(sample_chrom) && !is.null(std_chrom)
        
        if (has_data) {
          # Adjust mz_index if fragment_pare is TRUE to preserve original fragment number
          if (!is.null(original_frag_index)) {
            sample_chrom$mz_index <- original_frag_index
            std_chrom$mz_index <- original_frag_index
          }
          
          # Combine data
          sample_chrom$type <- "Sample"
          std_chrom$type <- "Standard"
          combined_data <- rbind(sample_chrom, std_chrom)
          
          # Create mz labels with actual m/z values (before stick transformation)
          combined_data$mz_label <- sprintf("mz%d: %.4f", combined_data$mz_index, combined_data$mz)
          
          # Apply stick transformation if needed
          if (stick) {
            combined_data <- combined_data |>
              group_by(rt, mz_label, type) |>
              summarize(intensity = max(intensity), .groups = "drop") |>
              group_by(mz_label, type) |>
              arrange(rt) |>
              mutate(
                rt_start = rt,
                rt_end = rt,
                intensity_start = 0,
                intensity_end = intensity
              ) |>
              ungroup()
          }
          
          # Filter to only rows with non-zero intensity (detected fragments)
          combined_data <- combined_data |>
            filter(intensity > 0)
        }
        
        # Recheck if we still have data after filtering
        has_data <- has_data && exists("combined_data") && nrow(combined_data) > 0
      
      if (!has_data && !force_plot) {
        next
      }
      
      if (has_data) {
        # Calculate absolute intensity limits
        max_sample_int <- max(abs(combined_data$intensity[combined_data$type == "Sample"]), na.rm = TRUE)
        max_standard_int <- max(abs(combined_data$intensity[combined_data$type == "Standard"]), na.rm = TRUE)
        y_limit <- max(max_sample_int, max_standard_int) * 1.05
        
        # Set plot intensities (negative for standard, positive for sample)
        combined_data <- combined_data |>
          mutate(plot_intensity = ifelse(type == "Standard", -intensity, intensity))
      } else {
        # No data - use default y limit
        y_limit <- 1000
      }
      
      # Add asterisks to marked m/z values (only if we have data)
      if (has_data && !is.na(row$asterisk)) {
        marked_mzs <- strsplit(row$asterisk, ", ")[[1]]
        for (marked_mz in marked_mzs) {
          mz_num <- as.numeric(gsub("mz", "", marked_mz))
          combined_data$mz_label <- ifelse(
            combined_data$mz_index == mz_num,
            paste0("**", combined_data$mz_label, " \\*****"),
            combined_data$mz_label
          )
        }
      }
      
      sample_id <- if (!is.na(sample_file)) gsub("\\.mzML$", "", basename(sample_file)) else sprintf("File%d", sample_idx)
      
      if (has_data) {
        # Define fixed color palette for mz indices (matches vp() function palette)
        mz_colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#D4A017", "#0072B2", "#D55E00", "#CC79A7", "#999999")
        names(mz_colors) <- paste0("mz", 0:8)
        
        # Create named vector for actual mz_labels present in the data
        mz_labels_present <- unique(combined_data$mz_label)
        mz_indices_present <- unique(combined_data$mz_index)
        color_mapping <- setNames(
          mz_colors[paste0("mz", mz_indices_present)],
          mz_labels_present
        )
        
        if (stick) {
          p_rtx <- ggplot(combined_data, aes(x = rt, y = plot_intensity, color = mz_label, group = interaction(mz_label, type))) +
            geom_segment(aes(xend = rt, yend = 0), linewidth = 0.4) +
            geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4)
        } else {
          p_rtx <- ggplot(combined_data, aes(x = rt, y = plot_intensity, color = mz_label, group = interaction(mz_label, type))) +
            geom_line(linewidth = 0.4) +
            geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4)
        }
        
        p_rtx <- p_rtx +
          scale_color_manual(values = color_mapping) +
          scale_y_continuous(
            expand = c(0, 0),
            limits = c(-y_limit, y_limit),
            labels = function(x) scales::label_scientific(digits = 2)(abs(x)),
            n.breaks = 8
          )
      } else {
        # Create empty plot for no-data case
        p_rtx <- ggplot() +
          annotate("text", x = mean(sample_rt_range), y = 0, 
                  label = "NO DATA", size = 6, color = "gray50", fontface = "bold") +
          geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
          scale_y_continuous(
            expand = c(0, 0),
            limits = c(-y_limit, y_limit),
            labels = function(x) scales::label_scientific(digits = 2)(abs(x)),
            n.breaks = 8
          )
      }
      
      # Apply x-axis scaling based on whether hard limits are used
      if (use_hard_limits) {
        p_rtx <- p_rtx +
          scale_x_continuous(
            expand = c(0, 0),
            limits = sample_rt_range,
            breaks = function(limits) seq(ceiling(limits[1] * 20) / 20, floor(limits[2] * 20) / 20, by = 0.05),
            minor_breaks = function(limits) seq(ceiling(limits[1] * 40) / 40, floor(limits[2] * 40) / 40, by = 0.025)
          )
      } else {
        p_rtx <- p_rtx +
          scale_x_continuous(
            expand = expansion(mult = c(0.05, 0.05), add = 0),
            breaks = function(limits) seq(ceiling(limits[1] * 20) / 20, floor(limits[2] * 20) / 20, by = 0.05),
            minor_breaks = function(limits) seq(ceiling(limits[1] * 40) / 40, floor(limits[2] * 40) / 40, by = 0.025)
          )
      }
      
      # Create subtitle with RT info
      if (rt_is_fallback) {
        subtitle_text <- sprintf("Sample: %s  |  Standard: %s  |  RT = NA (range: %.2f-%.2f min)", 
                                sample_id, standard_file, sample_rt_range[1], sample_rt_range[2])
      } else {
        subtitle_text <- sprintf("Sample: %s  |  Standard: %s  |  RT = %.2f min", 
                                sample_id, standard_file, mean(sample_rt_range))
      }
      
      p_rtx <- p_rtx +
        labs(
          title = short_name,
          subtitle = subtitle_text,
          x = "Retention Time (min)",
          y = sprintf("\u2190 Std  |  %s \u2192", tools::toTitleCase(study)),
          color = NULL
        ) +
        coord_cartesian(clip = "off") +
        theme_classic(base_size = 12) +
        theme(
          panel.grid.major.x = element_line(color = "gray90", linewidth = 0.2),
          panel.grid.minor.x = element_line(color = "gray90", linewidth = 0.2),
          plot.margin = margin(10, 10, 20, 20),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "top",
          legend.justification = "center",
          legend.direction = "horizontal",
          legend.text = ggtext::element_markdown(size = 4),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA),
          legend.key.size = unit(0.25, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.spacing.x = unit(0.02, "cm"),
          legend.box.margin = margin(0, 0, 0, 0),
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
      
      # Store plot
      plot_label <- sprintf("F%d_S%d", sample_idx, std_idx)
      plot_tag <- sprintf("F%d_S%d_%s", sample_idx, std_idx, id_val)
      if (study == "cadaver") {
        plot_tag <- paste0("C_", plot_tag)
      }
      
      compound_result$plots[[plot_label]] <- list(
        plot = p_rtx,
        sample_id = sample_id,
        standard_file = standard_file,
        plot_tag = plot_tag,
        rt_range = sample_rt_range
      )
      
      # Save RDS if requested
      if (save_rds && !is.null(rds_save_folder)) {
        if (exists("config") && !is.null(config$paths$validation_plot_directory)) {
          rds_dir <- file.path(config$paths$validation_plot_directory, rds_save_folder)
        } else {
          rds_dir <- file.path(output_dir, "RDS", rds_save_folder)
        }
        
        dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
        rds_path <- file.path(rds_dir, paste0(plot_tag, ".rds"))
        
        if (!file.exists(rds_path) || overwrite_rds) {
          individual_plot <- list(
            short_name = short_name,
            id = id_val,
            order = order_num,
            plot = p_rtx,
            sample_id = sample_id,
            standard_file = standard_file,
            plot_tag = plot_tag,
            rt_range = sample_rt_range
          )
          saveRDS(individual_plot, file = rds_path, compress = "gzip")
        }
      }
      }  # End else (standard comparison mode)
    }  # End if/else run_standard
  }
  
  # Close all cached mzML files
  for (file_name in ls(worker_cache)) {
    mzml_data <- worker_cache[[file_name]]
    if (!is.null(mzml_data$ms_data)) {
      try(mzR::close(mzml_data$ms_data), silent = TRUE)
    }
  }
  
  return(compound_result)
}

#' Batch RTX Validation: Chromatogram PDF Generator
#'
#' Creates retention time chromatogram plots comparing sample vs standard,
#' compiling them into a single PDF with 6 plots per page (3 rows x 2 columns).
#'
#' @param validation_list Tibble with columns: order, id, short_name, monoisotopic, compound_rt_range, 
#'   mz0-mz3, standards (comma-separated), file1-file6, asterisk, f1_rt_range-f6_rt_range
#' @param study Character string: "tumor" or "cadaver" to determine directory (default: "tumor")
#' @param iterate_through Integer: how many of the top sample files to process (default: 5)
#' @param output_dir Required character string: output directory for final PDF
#' @param pdf_name Character string: name of output PDF file (default: "rtx_validation.pdf")
#' @param ppm_tolerance Numeric mass tolerance in ppm (default: 5)
#' @param rt_lookup Character: "range" uses compound_rt_range (default), "sample" uses file-specific RT ranges, "window" uses file-specific RT ± window/2 as hard limits
#' @param window Numeric: window size in minutes for rt_lookup = "window" mode (e.g., 10/60 for 10 seconds). Only used when rt_lookup = "window".
#' @param stick Logical, whether to plot as vertical sticks (default: FALSE)
#' @param max_i Logical, whether to use maximum intensity (default: FALSE)
#' @param save_rds Logical, whether to save individual plot RDS files (default: TRUE)
#' @param rds_save_folder Character string: subfolder name within validation_plot_directory for saving individual plot RDS files
#' @param overwrite_rds Logical, whether to automatically overwrite existing RDS files without prompting (default: FALSE)
#' @param use_parallel Logical, whether to use parallel processing (default: FALSE)
#' @param n_cores Integer, number of cores to use for parallel processing (default: parallel::detectCores() - 1)
#' @param run_standard Logical, whether to process standard files and create mirror plots (default: TRUE). When FALSE, only sample chromatograms are generated, significantly improving speed.
#' @param fragment_pare Logical, whether to filter to only the top fragment (default: FALSE). When TRUE, looks for top_frag column in validation_list and only plots that specific fragment index.
#'
#' @return Named list of all plots (invisibly)
#'
#' @export
