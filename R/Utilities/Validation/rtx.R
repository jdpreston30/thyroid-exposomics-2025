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
                                   rt_lookup, ppm_tolerance, stick, max_i,
                                   save_rds, rds_save_folder, overwrite_rds, output_dir, study = "tumor") {
  
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
  extract_chrom_worker <- function(file_name, target_mzs, rt_range, ppm_tol, use_max, cache_env) {
    mzml_data <- get_mzml_data_worker(file_name, cache_env)
    if (is.null(mzml_data)) return(NULL)
    
    ms_data <- mzml_data$ms_data
    header_info <- mzml_data$header_info
    
    extract_xic <- function(ms_data, header_info, target_mz, ppm_tol, rt_range, use_max) {
      ms1_scans <- header_info[header_info$msLevel == 1, ]
      ms1_scans <- ms1_scans[ms1_scans$retentionTime >= rt_range[1] * 60 & 
                             ms1_scans$retentionTime <= rt_range[2] * 60, ]
      
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
      
      do.call(rbind, xic_data)
    }
    
    all_xics <- lapply(seq_along(target_mzs), function(i) {
      mz_val <- target_mzs[i]
      xic <- extract_xic(ms_data, header_info, mz_val, ppm_tol, rt_range, use_max)
      if (!is.null(xic)) {
        xic$mz_label <- paste0("mz", i - 1)
        xic
      }
    })
    
    do.call(rbind, all_xics)
  }
  
  # Initialize worker cache
  worker_cache <- new.env(hash = TRUE)
  
  # Extract compound info
  order_num <- row$order
  id_val <- row$id
  short_name <- row$short_name
  
  compound_result <- list(
    short_name = short_name,
    id = id_val,
    order = order_num,
    plots = list()
  )
  
  # Extract m/z values
  target_mzs <- row |>
    dplyr::select(matches("^mz[0-9]+$")) |>
    unlist() |>
    as.numeric() |>
    na.omit()
  
  if (length(target_mzs) == 0) return(compound_result)
  
  # Get RT range
  rt_range_str <- row$compound_rt_range
  if (is.na(rt_range_str)) return(compound_result)
  
  base_rt_range <- eval(parse(text = rt_range_str))
  base_rt_range_expanded <- c(base_rt_range[1] - 0.2, base_rt_range[2] + 0.2)
  
  # Parse standards
  standards <- strsplit(row$standards, ", ")[[1]]
  
  # Get sample files
  all_samples <- c(row$file1, row$file2, row$file3, row$file4, row$file5, row$file6)
  all_samples <- all_samples[!is.na(all_samples)]
  samples_to_process <- all_samples[1:min(iterate_through, length(all_samples))]
  
  # Process each sample
  for (sample_idx in seq_along(samples_to_process)) {
    sample_file <- samples_to_process[sample_idx]
    
    # Determine RT range
    sample_rt_range <- base_rt_range_expanded
    if (rt_lookup == "sample") {
      rt_range_col_name <- paste0("f", sample_idx, "_rt_range")
      if (rt_range_col_name %in% names(row) && !is.na(row[[rt_range_col_name]])) {
        rt_range_str <- row[[rt_range_col_name]]
        sample_rt_range <- eval(parse(text = rt_range_str))
        sample_rt_range <- c(sample_rt_range[1] - 0.2, sample_rt_range[2] + 0.2)
      }
    }
    
    # Process each standard
    for (std_idx in seq_along(standards)) {
      standard_file <- standards[std_idx]
      
      # Extract chromatograms
      sample_chrom <- extract_chrom_worker(sample_file, target_mzs, sample_rt_range, ppm_tolerance, max_i, worker_cache)
      std_chrom <- extract_chrom_worker(standard_file, target_mzs, sample_rt_range, ppm_tolerance, max_i, worker_cache)
      
      if (is.null(sample_chrom) || is.null(std_chrom)) next
      
      # Combine data
      sample_chrom$type <- "Sample"
      std_chrom$type <- "Standard"
      combined_data <- rbind(sample_chrom, std_chrom)
      
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
      
      # Normalize intensities
      max_int <- max(combined_data$intensity, na.rm = TRUE)
      combined_data <- combined_data |>
        mutate(plot_intensity = ifelse(type == "Standard", -intensity / max_int, intensity / max_int))
      
      # Create plot
      mz_labels <- paste0("mz", 0:(length(target_mzs) - 1))
      mz_colors <- c("mz0" = "black", "mz1" = "red", "mz2" = "gold", "mz3" = "blue")[mz_labels]
      
      sample_id <- gsub("\\.mzML$", "", basename(sample_file))
      
      p_rtx <- ggplot(combined_data, aes(x = rt, y = plot_intensity, color = mz_label, group = interaction(mz_label, type))) +
        geom_line(linewidth = 0.4) +
        geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
        scale_color_manual(values = mz_colors) +
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
        scale_x_continuous(
          expand = expansion(mult = c(0.05, 0.05), add = 0),
          breaks = function(limits) seq(ceiling(limits[1] * 20) / 20, floor(limits[2] * 20) / 20, by = 0.05),
          minor_breaks = function(limits) seq(ceiling(limits[1] * 40) / 40, floor(limits[2] * 40) / 40, by = 0.025)
        ) +
        labs(
          title = short_name,
          subtitle = sprintf("Sample: %s  |  Standard: %s  |  RT = %.2f min", sample_id, standard_file, mean(sample_rt_range)),
          x = "Retention Time (min)",
          y = sprintf("\u2190 Std  |  %s \u2192", tools::toTitleCase(study)),
          color = NULL
        ) +
        coord_cartesian(clip = "off") +
        theme_classic(base_size = 12) +
        theme(
          panel.grid.major.x = element_line(color = "gray90", linewidth = 0.2),
          panel.grid.minor.x = element_line(color = "gray90", linewidth = 0.2),
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
      
      # Store plot
      plot_label <- sprintf("F%d_S%d", sample_idx, std_idx)
      plot_tag <- sprintf("F%d_S%d_%s", sample_idx, std_idx, id_val)
      
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
    }
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
#' @param rt_lookup Character: "range" uses compound_rt_range (default), "sample" uses file-specific RT ranges
#' @param stick Logical, whether to plot as vertical sticks (default: FALSE)
#' @param max_i Logical, whether to use maximum intensity (default: FALSE)
#' @param save_rds Logical, whether to save individual plot RDS files (default: TRUE)
#' @param rds_save_folder Character string: subfolder name within validation_plot_directory for saving individual plot RDS files
#' @param overwrite_rds Logical, whether to automatically overwrite existing RDS files without prompting (default: FALSE)
#' @param use_parallel Logical, whether to use parallel processing (default: FALSE)
#' @param n_cores Integer, number of cores to use for parallel processing (default: parallel::detectCores() - 1)
#'
#' @return Named list of all plots (invisibly)
#'
#' @export
rtx <- function(validation_list,
                study = "tumor",
                iterate_through = 5,
                output_dir,
                pdf_name = "rtx_validation.pdf",
                ppm_tolerance = 5,
                rt_lookup = "range",
                stick = FALSE,
                max_i = FALSE,
                save_rds = TRUE,
                rds_save_folder = NULL,
                overwrite_rds = FALSE,
                use_parallel = FALSE,
                n_cores = NULL) {
  
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
  
  # Check for existing RDS files and prompt for overwrite if needed
  if (save_rds && !is.null(rds_save_folder) && !overwrite_rds) {
    # Determine RDS directory
    if (exists("config") && !is.null(config$paths$validation_plot_directory)) {
      rds_dir <- file.path(config$paths$validation_plot_directory, rds_save_folder)
    } else {
      rds_dir <- file.path(output_dir, "RDS", rds_save_folder)
    }
    
    # Check if directory exists and contains .rds files
    if (dir.exists(rds_dir)) {
      existing_rds <- list.files(rds_dir, pattern = "\\.rds$", full.names = FALSE)
      if (length(existing_rds) > 0) {
        cat(sprintf("\n⚠️  Found %d existing RDS files in: %s\n", length(existing_rds), rds_dir))
        cat("Do you want to overwrite them? (Y/N): ")
        response <- toupper(trimws(readline()))
        
        if (response != "Y") {
          cat("Skipping RDS file saving. Set overwrite_rds = TRUE to automatically overwrite.\n")
          save_rds <- FALSE
        } else {
          cat("Proceeding with overwrite...\n")
        }
      }
    }
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
  
  # Initialize mzML file cache
  mzml_cache <- new.env(hash = TRUE)
  
  # Helper function to get cached mzML data or read if not cached
  get_mzml_data <- function(file_name) {
    if (exists(file_name, envir = mzml_cache)) {
      return(mzml_cache[[file_name]])
    }
    
    mzml_file <- paste0(file_name, ".mzML")
    mzml_path <- file.path(mzml_dir, mzml_file)
    
    if (!file.exists(mzml_path)) {
      warning(sprintf("mzML file not found: %s - skipping", mzml_path))
      return(NULL)
    }
    
    ms_data <- mzR::openMSfile(mzml_path)
    header_info <- mzR::header(ms_data)
    
    # Cache the data
    mzml_cache[[file_name]] <- list(
      ms_data = ms_data,
      header_info = header_info,
      path = mzml_path
    )
    
    return(mzml_cache[[file_name]])
  }
  
  # Helper function to extract chromatogram data
  extract_chromatogram <- function(file_name, target_mzs, rt_range, ppm_tol, use_max) {
    mzml_data <- get_mzml_data(file_name)
    
    if (is.null(mzml_data)) {
      return(NULL)
    }
    
    ms_data <- mzml_data$ms_data
    header_info <- mzml_data$header_info
    
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
    # Don't close ms_data here since it's cached and may be reused
    
    return(chromatogram_data)
  }
  
  # Initialize plot storage
  compound_plots <- list()
  
  # Main iteration loop
  cat(sprintf("\nStarting RTX validation for %d compounds...\n", nrow(validation_list)))
  
  # Determine if using parallel processing
  if (use_parallel) {
    if (!requireNamespace("foreach", quietly = TRUE) || 
        !requireNamespace("doParallel", quietly = TRUE)) {
      warning("foreach and doParallel packages required for parallel processing. Falling back to sequential.")
      use_parallel <- FALSE
    }
  }
  
  if (use_parallel) {
    # Setup parallel backend
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
    }
    n_cores <- min(n_cores, parallel::detectCores() - 1, nrow(validation_list))
    
    cat(sprintf("Using parallel processing with %d cores\n", n_cores))
    
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    
    # Export necessary objects and functions to workers
    parallel::clusterExport(cl, c("mzml_dir", "iterate_through", "rt_lookup", 
                                   "stick", "max_i", "ppm_tolerance",
                                   "save_rds", "rds_save_folder", "overwrite_rds",
                                   "output_dir", "study", "config", "process_single_compound"),
                           envir = environment())
    
    # Load required packages on each worker
    parallel::clusterEvalQ(cl, {
      library(mzR)
      library(ggplot2)
      library(dplyr)
      library(ggtext)
    })
    
    # Process compounds in parallel
    compound_results <- foreach::foreach(
      row_idx = 1:nrow(validation_list),
      .packages = c("mzR", "ggplot2", "dplyr", "ggtext"),
      .errorhandling = "pass"
    ) %dopar% {
      row <- validation_list[row_idx, ]
      
      # Process this compound (worker function defined below)
      process_single_compound(
        row = row,
        row_idx = row_idx,
        total_rows = nrow(validation_list),
        mzml_dir = mzml_dir,
        iterate_through = iterate_through,
        rt_lookup = rt_lookup,
        ppm_tolerance = ppm_tolerance,
        stick = stick,
        max_i = max_i,
        save_rds = save_rds,
        rds_save_folder = rds_save_folder,
        overwrite_rds = overwrite_rds,
        output_dir = output_dir,
        study = study
      )
    }
    
    # Stop cluster
    parallel::stopCluster(cl)
    
    # Convert list to named list by compound ID
    compound_plots <- setNames(compound_results, sapply(validation_list$id, as.character))
    
  } else {
    # Sequential processing
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
    
    # Get base RT range
    rt_range_str <- row$compound_rt_range
    if (is.na(rt_range_str)) {
      warning(sprintf("No RT range found for ID '%s' - skipping", id_val))
      next
    }
    
    base_rt_range <- eval(parse(text = rt_range_str))
    base_rt_range_expanded <- c(base_rt_range[1] - 0.2, base_rt_range[2] + 0.2)
    
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
      
      # Determine RT range for this sample
      sample_rt_range <- base_rt_range_expanded
      
      if (rt_lookup == "sample") {
        rt_range_col_name <- paste0("f", sample_idx, "_rt_range")
        if (rt_range_col_name %in% names(row) && !is.na(row[[rt_range_col_name]])) {
          rt_range_str <- row[[rt_range_col_name]]
          sample_rt_range <- eval(parse(text = rt_range_str))
          cat(sprintf("    Using file-specific RT range: [%.2f, %.2f]\n", 
                      sample_rt_range[1], sample_rt_range[2]))
        }
      }
      
      # Extract chromatogram for sample
      sample_chrom <- extract_chromatogram(sample_file, target_mzs, sample_rt_range, ppm_tolerance, max_i)
      
      if (is.null(sample_chrom)) {
        warning(sprintf("Skipping sample %s - data extraction failed", sample_file))
        next
      }
      
      sample_chrom$type <- "Sample"
      sample_chrom$plot_intensity <- sample_chrom$intensity
      
      # Iterate through standards
      for (std_idx in seq_along(standards)) {
        standard_file <- standards[std_idx]
        
        cat(sprintf("    Standard %d/%d: %s\n", std_idx, length(standards), standard_file))
        
        # Extract chromatogram for standard
        standard_chrom <- extract_chromatogram(standard_file, target_mzs, sample_rt_range, ppm_tolerance, max_i)
        
        if (is.null(standard_chrom)) {
          warning(sprintf("Skipping standard %s - data extraction failed", standard_file))
          next
        }
        
        standard_chrom$type <- "Standard"
        standard_chrom$plot_intensity <- -standard_chrom$intensity
        
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
        
        sample_id <- sample_file
        
        # Define fixed color palette for mz indices
        mz_colors <- c("mz0" = "black", "mz1" = "red", "mz2" = "gold", "mz3" = "blue")
        
        # Create named vector for actual mz_labels present in the data
        mz_labels_present <- unique(combined_chrom$mz_label)
        mz_indices_present <- unique(combined_chrom$mz_index)
        color_mapping <- setNames(
          mz_colors[paste0("mz", mz_indices_present)],
          mz_labels_present
        )
        
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
          scale_color_manual(values = color_mapping) +
          scale_x_continuous(limits = sample_rt_range, expand = expansion(mult = c(0.05, 0.05), add = 0),
                           breaks = function(limits) seq(ceiling(limits[1] * 20) / 20, floor(limits[2] * 20) / 20, by = 0.05),
                           minor_breaks = function(limits) seq(ceiling(limits[1] * 40) / 40, floor(limits[2] * 40) / 40, by = 0.025)) +
          scale_y_continuous(
            expand = c(0, 0),
            limits = c(-y_limit_chrom, y_limit_chrom),
            labels = function(x) abs(x),
            n.breaks = 8
          ) +
          labs(
            title = short_name,
            subtitle = sprintf("Sample: %s  |  Standard: %s  |  RT = %.2f min", sample_id, standard_file, mean(sample_rt_range)),
            x = "Retention Time (min)",
            y = sprintf("\u2190 Std  |  %s \u2192", tools::toTitleCase(study)),
            color = NULL
          ) +
          coord_cartesian(clip = "off") +
          theme_classic(base_size = 12) +
          theme(
            panel.grid.major.x = element_line(color = "gray90", linewidth = 0.2),
            panel.grid.minor.x = element_line(color = "gray90", linewidth = 0.2),
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
        
        # Store plot
        plot_label <- sprintf("F%d_S%d", sample_idx, std_idx)
        plot_tag <- sprintf("F%d_S%d_%s", sample_idx, std_idx, id_val)
        
        compound_plots[[id_val]]$plots[[plot_label]] <- list(
          plot = p_rtx,
          sample_id = sample_id,
          standard_file = standard_file,
          plot_tag = plot_tag,
          rt_range = sample_rt_range
        )
        
        cat(sprintf("      Created plot: %s\n", plot_label))
        
        # Save individual plot RDS if requested
        if (save_rds && !is.null(rds_save_folder)) {
          # Get validation_plot_directory from config
          if (exists("config") && !is.null(config$paths$validation_plot_directory)) {
            rds_dir <- file.path(config$paths$validation_plot_directory, rds_save_folder)
          } else {
            # Fallback to local directory if config not available
            rds_dir <- file.path(output_dir, "RDS", rds_save_folder)
          }
          
          dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
          
          # Check if file already exists (skip if overwrite_rds is FALSE)
          rds_path <- file.path(rds_dir, paste0(plot_tag, ".rds"))
          
          if (file.exists(rds_path) && !overwrite_rds) {
            cat(sprintf("      RDS file already exists, skipping: %s\n", plot_tag))
          } else {
            # Create individual plot object with full metadata
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
            
            # Save with plot_tag as filename
            saveRDS(individual_plot, file = rds_path, compress = "gzip")
            cat(sprintf("      Saved RDS: %s\n", plot_tag))
          }
        }
      }
    }
  }
  } # End of sequential/parallel if-else
  
  # Compile into PDF
  cat(sprintf("\nCompiling %d compounds into PDF...\n", length(compound_plots)))
  
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
    
    # Modify plots for PDF: add red tags to subtitles
    pdf_plots <- lapply(plot_labels, function(label) {
      plot_info <- plots[[label]]
      p <- plot_info$plot
      sample_id <- plot_info$sample_id
      standard_file <- plot_info$standard_file
      plot_tag <- plot_info$plot_tag
      
      # Add plot tag to subtitle
      new_subtitle <- sprintf("Sample: %s  |  Standard: %s  |  RT = %.3f min  |  %s",
                              sample_id, standard_file, mean(compound$plots[[label]]$rt_range), plot_tag)
      
      p <- p +
        labs(subtitle = new_subtitle) +
        theme(
          plot.subtitle = ggtext::element_markdown(hjust = 0.5, face = "italic", size = 6,
                                                   color = "black", lineheight = 1.2, margin = margin(0, 0, 3, 0))
        )
      
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
  
  # Copy PDF to validation_plots subdirectory (same location as RDS files)
  if (save_rds && !is.null(rds_save_folder) && pdf_closed) {
    pdf_dest_dir <- if (exists("config") && !is.null(config$paths$validation_plot_directory)) {
      file.path(config$paths$validation_plot_directory, rds_save_folder)
    } else {
      file.path(output_dir, "RDS", rds_save_folder)
    }
    pdf_dest_path <- file.path(pdf_dest_dir, basename(pdf_path))
    tryCatch({
      file.copy(pdf_path, pdf_dest_path, overwrite = TRUE)
      cat(sprintf("PDF copied to: %s\n", pdf_dest_path))
    }, error = function(e) {
      cat(sprintf("Warning: Could not copy PDF to subdirectory: %s\n", e$message))
    })
  }
  
  # Close all cached mzML files
  cat(sprintf("\nClosing %d cached mzML files...\n", length(ls(envir = mzml_cache))))
  for (file_name in ls(envir = mzml_cache)) {
    tryCatch({
      mzR::close(mzml_cache[[file_name]]$ms_data)
    }, error = function(e) {
      # Silently ignore if already closed
    })
  }
  
  # Report RDS saving summary
  if (save_rds && !is.null(rds_save_folder)) {
    total_plots <- sum(sapply(compound_plots, function(x) length(x$plots)))
    cat(sprintf("\nSaved %d individual plot RDS files to: %s\n", 
                total_plots, 
                if (exists("config") && !is.null(config$paths$validation_plot_directory)) {
                  file.path(config$paths$validation_plot_directory, rds_save_folder)
                } else {
                  file.path(output_dir, "RDS", rds_save_folder)
                }))
  }
  
  invisible(compound_plots)
}
