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
        xic$mz_index <- i - 1
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
      
      # Filter to only rows with non-zero intensity (detected fragments)
      combined_data <- combined_data |>
        filter(intensity > 0)
      
      if (nrow(combined_data) == 0) {
        next
      }
      
      # Calculate absolute intensity limits
      max_sample_int <- max(abs(combined_data$intensity[combined_data$type == "Sample"]), na.rm = TRUE)
      max_standard_int <- max(abs(combined_data$intensity[combined_data$type == "Standard"]), na.rm = TRUE)
      y_limit <- max(max_sample_int, max_standard_int) * 1.05
      
      # Set plot intensities (negative for standard, positive for sample)
      combined_data <- combined_data |>
        mutate(plot_intensity = ifelse(type == "Standard", -intensity, intensity))
      
      # Create mz labels with actual m/z values
      combined_data$mz_label <- sprintf("mz%d: %.4f", combined_data$mz_index, combined_data$mz)
      
      # Add asterisks to marked m/z values
      if (!is.na(row$asterisk)) {
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
      
      sample_id <- gsub("\\.mzML$", "", basename(sample_file))
      
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
        scale_color_viridis_d(option = "turbo", end = 0.9) +
        scale_y_continuous(
          expand = c(0, 0),
          limits = c(-y_limit, y_limit),
          labels = function(x) abs(x),
          n.breaks = 8
        ) +
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
#' @param save_compiled_rds Logical, whether to save the compiled compound_plots object to OneDrive (default: FALSE)
#' @param use_parallel Logical, whether to use parallel processing (default: FALSE)
#' @param n_cores Integer, number of cores to use for parallel processing (default: parallel::detectCores() - 1)
#'
#' @return Named list of all plots (invisibly)
#'
#' @export
rtx <- function(validation_list,
                study = "tumor",
                iterate_through = 5,
                output_dir = NULL,
                ppm_tolerance = 5,
                rt_lookup = "range",
                stick = FALSE,
                max_i = FALSE,
                save_rds = TRUE,
                rds_save_folder = NULL,
                overwrite_rds = FALSE,
                save_compiled_rds = FALSE,
                use_parallel = FALSE,
                n_cores = NULL,
                skip_if_disabled = TRUE) {
  
  # Check output_dir only if needed (when config is not available)
  if (is.null(output_dir) && save_rds && !is.null(rds_save_folder)) {
    if (!exists("config", envir = .GlobalEnv) || is.null(get("config", envir = .GlobalEnv)$paths$validation_plot_directory)) {
      stop("output_dir is required when config$paths$validation_plot_directory is not set")
    }
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
  
  # Track total runtime for entire function
  function_start_time <- Sys.time()
  
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
    
    # Calculate total steps for progress tracking
    total_steps <- 0
    for (i in 1:nrow(validation_list)) {
      row <- validation_list[i, ]
      all_samples <- c(row$file1, row$file2, row$file3, row$file4, row$file5, row$file6)
      all_samples <- all_samples[!is.na(all_samples)]
      samples_count <- min(iterate_through, length(all_samples))
      standards_count <- length(strsplit(row$standards, ", ")[[1]])
      total_steps <- total_steps + (samples_count * standards_count)
    }
    
    cat(sprintf("\n▶ Starting parallel processing: %d compounds → %d total plots (using %d cores)\n\n", 
                nrow(validation_list), total_steps, n_cores))
    
    start_time <- Sys.time()
    
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
    
    # Process compounds in parallel with progress tracking
    compounds_completed <- 0
    compound_results <- foreach::foreach(
      row_idx = 1:nrow(validation_list),
      .packages = c("mzR", "ggplot2", "dplyr", "ggtext"),
      .errorhandling = "pass",
      .combine = function(a, b) {
        compounds_completed <<- compounds_completed + 1
        pct <- compounds_completed / nrow(validation_list) * 100
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        eta_secs <- (elapsed / compounds_completed) * (nrow(validation_list) - compounds_completed)
        cat(sprintf("\r[Parallel] %d/%d compounds (%.0f%%) | Time: %02d:%02d:%02d | ETA: %02d:%02d:%02d     ",
                    compounds_completed, nrow(validation_list), pct,
                    floor(elapsed / 3600), floor((elapsed %% 3600) / 60), floor(elapsed %% 60),
                    floor(eta_secs / 3600), floor((eta_secs %% 3600) / 60), floor(eta_secs %% 60)))
        flush.console()
        c(a, list(b))
      },
      .init = list()
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
    
    # Final progress message
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("\r[Parallel] ✓ Completed: %d compounds, %d plots | Total time: %02d:%02d:%02d\n\n",
                nrow(validation_list), total_steps,
                floor(elapsed / 3600), floor((elapsed %% 3600) / 60), floor(elapsed %% 60)))
    
    # Convert list to named list by compound ID
    compound_plots <- setNames(compound_results, sapply(validation_list$id, as.character))
    
  } else {
    # Sequential processing with progress bar
    total_compounds <- nrow(validation_list)
    start_time <- Sys.time()
    
    # Calculate total steps (compounds × samples × standards)
    total_steps <- 0
    for (i in 1:nrow(validation_list)) {
      row <- validation_list[i, ]
      all_samples <- c(row$file1, row$file2, row$file3, row$file4, row$file5, row$file6)
      all_samples <- all_samples[!is.na(all_samples)]
      samples_count <- min(iterate_through, length(all_samples))
      standards_count <- length(strsplit(row$standards, ", ")[[1]])
      total_steps <- total_steps + (samples_count * standards_count)
    }
    
    # Use environment to hold mutable counter
    progress_env <- new.env()
    progress_env$current_step <- 0
    
    cat(sprintf("\n▶ Starting sequential processing: %d compounds → %d total plots\n\n", 
                total_compounds, total_steps))
    
    # Function to update progress display
    update_progress <- function(force_refresh = FALSE) {
      current_step <- progress_env$current_step
      current_time <- Sys.time()
      elapsed <- as.numeric(difftime(current_time, start_time, units = "secs"))
      
      if (current_step > 0) {
        avg_time_per_step <- elapsed / current_step
        remaining_steps <- total_steps - current_step
        eta_secs <- avg_time_per_step * remaining_steps
        eta_str <- sprintf("%02d:%02d:%02d", 
                          floor(eta_secs / 3600), 
                          floor((eta_secs %% 3600) / 60), 
                          floor(eta_secs %% 60))
      } else {
        eta_str <- "calculating..."
      }
      
      elapsed_str <- sprintf("%02d:%02d:%02d", 
                            floor(elapsed / 3600), 
                            floor((elapsed %% 3600) / 60), 
                            floor(elapsed %% 60))
      
      # Progress bar (50 blocks = ~15 steps per block for 756 total)
      pct <- current_step / total_steps
      bar_width <- 50
      filled <- floor(pct * bar_width)
      bar <- paste0(paste0(rep("█", filled), collapse = ""), 
                    paste0(rep("░", bar_width - filled), collapse = ""))
      
      # Get current compound info
      current_compound_info <- if (exists("current_id") && exists("current_name")) {
        sprintf("%s (%s)", current_name, current_id)
      } else {
        "Initializing..."
      }
      
      # Single line progress update
      cat(sprintf("\r[%s] %3.0f%% | Step %d/%d | %s | Time: %s | ETA: %s          ", 
                  bar, pct * 100, current_step, total_steps,
                  current_compound_info, elapsed_str, eta_str))
      
      flush.console()
    }
    
    for (row_idx in 1:nrow(validation_list)) {
    row <- validation_list[row_idx, ]
    
    order_num <- row$order
    id_val <- row$id
    short_name <- row$short_name
    
    # Set current compound for display
    current_id <<- id_val
    current_name <<- short_name
    
    # Calculate plots for this compound
    all_samples <- c(row$file1, row$file2, row$file3, row$file4, row$file5, row$file6)
    all_samples <- all_samples[!is.na(all_samples)]
    
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
    
    # Iterate through sample files (silent processing)
    for (sample_idx in seq_along(samples_to_process)) {
      sample_file <- samples_to_process[sample_idx]
      
      # Determine RT range for this sample
      sample_rt_range <- base_rt_range_expanded
      
      if (rt_lookup == "sample") {
        rt_range_col_name <- paste0("f", sample_idx, "_rt_range")
        if (rt_range_col_name %in% names(row) && !is.na(row[[rt_range_col_name]])) {
          rt_range_str <- row[[rt_range_col_name]]
          sample_rt_range <- eval(parse(text = rt_range_str))
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
      
      # Iterate through standards (silent processing)
      for (std_idx in seq_along(standards)) {
        standard_file <- standards[std_idx]
        
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
              paste0("**", combined_chrom$mz_label, " \\*****"),
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
            y = sprintf("← Standard  |  %s →     ", tools::toTitleCase(study)),
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
        
        compound_plots[[id_val]]$plots[[plot_label]] <- list(
          plot = p_rtx,
          sample_id = sample_id,
          standard_file = standard_file,
          plot_tag = plot_tag,
          rt_range = sample_rt_range
        )
        
        # Increment and update progress after each plot
        progress_env$current_step <- progress_env$current_step + 1
        update_progress()
        
        # Save individual plot RDS if requested (silent)
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
            
            # Save with plot_tag as filename (silent)
            saveRDS(individual_plot, file = rds_path, compress = "gzip")
          }
        }
      }
    }
  }
  } # End of sequential/parallel if-else
  
  # Final progress update for sequential mode
  if (!use_parallel) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    elapsed_str <- sprintf("%02d:%02d:%02d", 
                          floor(elapsed / 3600), 
                          floor((elapsed %% 3600) / 60), 
                          floor(elapsed %% 60))
    cat(sprintf("\r[%s] 100%% | Completed: %d compounds, %d plots | Total time: %s\n\n", 
                paste0(rep("█", 30), collapse = ""), 
                total_compounds, progress_env$current_step, elapsed_str))
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
    local_dir <- if (exists("config") && !is.null(config$paths$validation_plot_directory)) {
      file.path(config$paths$validation_plot_directory, rds_save_folder)
    } else {
      file.path(output_dir, "RDS", rds_save_folder)
    }
    cat(sprintf("\nSaved %d individual plot RDS files to: %s\n", total_plots, local_dir))
    
    # Bulk transfer to OneDrive if configured
    if (exists("config") && !is.null(config$paths$validation_plot_directory_onedrive)) {
      onedrive_dir <- file.path(config$paths$validation_plot_directory_onedrive, rds_save_folder)
      cat(sprintf("\n📤 Transferring RDS files to OneDrive backup...\n"))
      
      dir.create(onedrive_dir, recursive = TRUE, showWarnings = FALSE)
      
      local_files <- list.files(local_dir, pattern = "\\.rds$", full.names = TRUE)
      
      success_count <- 0
      fail_count <- 0
      failed_files <- character(0)
      
      # First pass: transfer all files
      for (local_file in local_files) {
        dest_file <- file.path(onedrive_dir, basename(local_file))
        
        # Attempt copy (may timeout but still succeed)
        suppressWarnings({
          tryCatch({
            file.copy(local_file, dest_file, overwrite = TRUE)
          }, error = function(e) {
            # Ignore errors, will verify existence below
          })
        })
        
        # Verify file actually exists at destination (more reliable than file.copy return value)
        if (file.exists(dest_file)) {
          success_count <- success_count + 1
        } else {
          fail_count <- fail_count + 1
          failed_files <- c(failed_files, local_file)
        }
      }
      
      cat(sprintf("✓ Transferred %d/%d files to OneDrive\n", success_count, length(local_files)))
      
      # Retry failed transfers with exponential backoff
      if (fail_count > 0) {
        cat(sprintf("⚠️  %d files failed to transfer, retrying...\n", fail_count))
        max_retries <- 3
        retry_count <- 0
        
        while (length(failed_files) > 0 && retry_count < max_retries) {
          retry_count <- retry_count + 1
          wait_time <- 2^retry_count  # 2, 4, 8 seconds
          cat(sprintf("   Retry %d/%d (waiting %d seconds)...\n", retry_count, max_retries, wait_time))
          Sys.sleep(wait_time)
          
          still_failed <- character(0)
          for (local_file in failed_files) {
            dest_file <- file.path(onedrive_dir, basename(local_file))
            
            # Attempt copy
            suppressWarnings({
              tryCatch({
                file.copy(local_file, dest_file, overwrite = TRUE)
              }, error = function(e) {
                # Ignore errors, will verify existence below
              })
            })
            
            # Verify file actually exists at destination
            if (file.exists(dest_file)) {
              success_count <- success_count + 1
              fail_count <- fail_count - 1
              cat(sprintf("   ✓ Retry succeeded: %s\n", basename(local_file)))
            } else {
              still_failed <- c(still_failed, local_file)
            }
          }
          failed_files <- still_failed
        }
        
        if (fail_count > 0) {
          cat(sprintf("⚠️  %d files still failed after %d retries\n", fail_count, max_retries))
          cat(sprintf("   Failed files remain in: %s\n", local_dir))
        } else {
          cat(sprintf("✓ All retries succeeded! %d/%d files transferred\n", success_count, length(local_files)))
        }
      }
      

      
      # Clean up local temp files after successful transfer
      if (fail_count == 0) {
        cat(sprintf("\n🗑️  Cleaning up local temp files...\n"))
        deleted_count <- 0
        for (local_file in local_files) {
          if (file.remove(local_file)) {
            deleted_count <- deleted_count + 1
          }
        }
        cat(sprintf("✓ Removed %d temp RDS files from local directory\n", deleted_count))
        

        
        # Try to remove the run-specific folder if empty (ignore .DS_Store and system files)
        remaining_files <- list.files(local_dir, all.files = TRUE, no.. = TRUE)
        # Filter out .DS_Store and other system files
        remaining_files <- remaining_files[!remaining_files %in% c(".DS_Store", "Thumbs.db")]
        if (length(remaining_files) == 0) {
          # Remove .DS_Store if it exists before removing directory
          ds_store <- file.path(local_dir, ".DS_Store")
          if (file.exists(ds_store)) file.remove(ds_store)
          
          unlink(local_dir, recursive = TRUE)
          cat(sprintf("✓ Removed empty temp directory: %s\n", basename(local_dir)))
          
          # Try to remove parent temp_RDS folder if completely empty
          parent_temp_dir <- dirname(local_dir)
          if (basename(parent_temp_dir) == "temp_RDS") {
            parent_remaining <- list.files(parent_temp_dir, all.files = TRUE, no.. = TRUE)
            parent_remaining <- parent_remaining[!parent_remaining %in% c(".DS_Store", "Thumbs.db")]
            if (length(parent_remaining) == 0) {
              # Remove .DS_Store if it exists
              ds_store <- file.path(parent_temp_dir, ".DS_Store")
              if (file.exists(ds_store)) file.remove(ds_store)
              
              unlink(parent_temp_dir, recursive = TRUE)
              cat(sprintf("✓ Removed empty parent temp directory: temp_RDS\n"))
            }
          }
        }
      }
    }
  }
  
  # Save compiled compound_plots object if requested
  if (save_compiled_rds && !is.null(rds_save_folder)) {
    if (exists("config") && !is.null(config$paths$validation_plot_directory_onedrive)) {
      compiled_rds_dir <- config$paths$validation_plot_directory_onedrive
      dir.create(compiled_rds_dir, recursive = TRUE, showWarnings = FALSE)
      compiled_rds_path <- file.path(compiled_rds_dir, paste0(rds_save_folder, "_compiled.rds"))
      
      cat(sprintf("\n💾 Saving compiled compound_plots object...\n"))
      saveRDS(compound_plots, compiled_rds_path)
      cat(sprintf("✓ Saved compiled object to: %s\n", compiled_rds_path))
    } else {
      warning("save_compiled_rds = TRUE but config$paths$validation_plot_directory_onedrive is not set. Skipping compiled RDS save.")
    }
  }
  # Print final total time for entire workflow
  total_elapsed <- as.numeric(difftime(Sys.time(), function_start_time, units = "secs"))
  cat(sprintf("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"))
  cat(sprintf("✓ RTX Plot Generation Complete\n"))
  cat(sprintf("  Total Runtime: %02d:%02d:%02d (plot generation + RDS saving + OneDrive backup)\n",
              floor(total_elapsed / 3600), 
              floor((total_elapsed %% 3600) / 60), 
              floor(total_elapsed %% 60)))
  cat(sprintf("  Generated %d compounds with %d total plots\n",
              length(compound_plots),
              sum(sapply(compound_plots, function(x) length(x$plots)))))
  cat(sprintf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"))
  
  return(compound_plots)
}
