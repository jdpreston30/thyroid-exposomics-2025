#' Compile Validation Plots into PDF
#'
#' Takes a compound_plots object (output from rtx()) and compiles it into a PDF
#'
#' @param compound_plots List of compound plots from rtx() function
#' @param output_dir Directory to save PDF
#' @param pdf_name Name of PDF file
#' @param add_plot_tags Logical, whether to add plot tags to subtitles (default TRUE)
#' @param external_subfolder Optional subfolder name within OneDrive validation_plots directory
#'
#' @return Path to created PDF file
#' @export
compile_validation_pdf <- function(compound_plots, 
                                   output_dir, 
                                   pdf_name,
                                   add_plot_tags = TRUE,
                                   skip_if_disabled = TRUE,
                                   external_subfolder = NULL) {
  
  if (is.null(compound_plots) || length(compound_plots) == 0) {
    stop("No plots provided. compound_plots must be a non-empty list.")
  }
  
  cat(sprintf("\n📄 Compiling %d compounds into PDF...\n", length(compound_plots)))
  
  pdf_path <- file.path(output_dir, pdf_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Ensure no lingering graphics devices before opening PDF
#! BOUNDED, and on.exit-guarded. An unbounded `while (dev.cur() > 1) dev.off()` spins forever when a device refuses to close -- which is precisely what a failed PDF write produces ("internal read error in PDF_endpage"). On 2026-09-04 that hung the run, and the device it left behind later segfaulted script 09 inside ggplotGrob's text measurement. The on.exit guarantees the sweep runs even if this function errors out, so a failed diagnostic PDF can never poison a downstream script. Nothing reads this PDF; it must never be able to stop the pipeline.
  for (.i in seq_len(20)) { if (dev.cur() <= 1) break; try(dev.off(), silent = TRUE) }
  on.exit({ for (.i in seq_len(20)) { if (dev.cur() <= 1) break; try(dev.off(), silent = TRUE) } }, add = TRUE)
  
#! OPEN DEFECT (punch list 133): some compiles fail here with one `write failed`, then every later page cascades `internal read error in PDF_endpage` off the broken device. It is TRANSIENT, not content-dependent: on 2026-09-09 pt1 failed on its first compound's page 2 while pt2 (same grob content, same markdown, primes and asterisks) compiled all 20 compounds clean, and which variant block fails has varied across runs. RULED OUT: disk space (reproduced at 26% capacity, 35 GB free; the 09-08 run's 92% was coincidence); plot count (pt1 dies ~13 plots in, IARC cadaver does 198); any specific compound (CP3017 renders in one block and fails in another); grob content (pt2 proves it); the `output_dir` argument (rtx() only uses it for a fallback RDS path); a leaked cluster (rtx.R:291 stops it). Best remaining lead: the failing compile is the one that follows the LONGEST rtx() block (pt1, 47 min, 240 plots) and starts the instant the 240-file RDS transfer to the exFAT drive returns -- so check memory pressure and drive I/O at the moment pdf() opens, and capture the real error with `tryCatch(..., error = function(e) print(e))` around the first `grid.draw`, since `write failed` is R's generic message for any failed fwrite on the device. Harmless today only because of the sweep above plus the try() wrapper at every call site -- fix the trigger, do not remove the containment.
  pdf(pdf_path, width = 8.5, height = 11, family = "Helvetica", onefile = TRUE)
  
  first_page <- TRUE
  .diag_done <- FALSE
  
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
#! [["plot"]] not $plot: the entry no longer has a `plot` element once rtx() has dropped it, and `$` partial-matches to `plot_tag` and returns the tag STRING -- which passed is.null() and reached arrangeGrob() as a character. Caught by the unit test for the 134 refactor.
      p <- plot_info[["plot"]]
#! rtx() no longer returns ggplot objects (see process_single_compound.R) -- each entry carries rds_name/rds_folder instead, and the grob is read back here one at a time so the parent never holds a whole block in memory. By the time compile runs, rtx() has already transferred the RDS to validation_plot_directory_onedrive and deleted the local copy, so look there first and fall back to the local staging dir only if the transfer was skipped.
      if (is.null(p) && !is.null(plot_info$rds_name)) {
        .cands <- c(
          if (exists("config") && !is.null(config$paths$validation_plot_directory_onedrive))
            file.path(config$paths$validation_plot_directory_onedrive, plot_info$rds_folder, plot_info$rds_name),
          if (exists("config") && !is.null(config$paths$validation_plot_directory))
            file.path(config$paths$validation_plot_directory, plot_info$rds_folder, plot_info$rds_name),
          file.path(output_dir, "RDS", plot_info$rds_folder, plot_info$rds_name)
        )
        .hit <- .cands[file.exists(.cands)][1]
        if (is.na(.hit)) stop(sprintf("compile_validation_pdf: grob %s not found in any of: %s",
                                      plot_info$rds_name, paste(.cands, collapse = " | ")))
        .obj <- readRDS(.hit)
        p <- .obj$plot
        if (is.null(plot_info$subtitle_rt)) plot_info$subtitle_rt <- .obj$subtitle_rt
      }
      
      if (add_plot_tags) {
        sample_id <- plot_info$sample_id
        standard_file <- plot_info$standard_file
        plot_tag <- plot_info$plot_tag
        rt_range <- plot_info$rt_range
        
        # Add plot tag to subtitle
#! Same defect as remove_standard()/zoom_fragment() had: mean(rt_range) is the midpoint of two %.2f-rounded endpoints, not the measured RT. Diagnostic PDFs only -- nothing outside script 08 reads initial_compile/ -- but there is no reason for it to print a different number than the supplement.
        new_subtitle <- sprintf("Sample: %s  |  Standard: %s  |  RT = %.3f min  |  %s",
                               sample_id, standard_file, subtitle_rt_of(plot_info), plot_tag)
        
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
      
      # Pad with blank plots to maintain 3x2 grid layout
      n_plots_on_page <- length(page_plots)
      if (n_plots_on_page < plots_per_page) {
        blank_plot <- ggplot2::ggplot() + ggplot2::theme_void()
        n_blanks_needed <- plots_per_page - n_plots_on_page
        for (i in 1:n_blanks_needed) {
          page_plots[[n_plots_on_page + i]] <- blank_plot
        }
      }
      
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
#! Diagnostics for punch list 133. `write failed` is the pdf() device's generic message for any short fwrite, so on the FIRST failure of a compile dump everything that could explain it: the full condition and call, the device stack, how big the partial PDF is, free disk under it, and system memory pressure. Later pages in the same compile always fail with PDF_endpage because the device is already broken, so only the first one carries information.
        if (!isTRUE(.diag_done)) {
          .diag_done <<- TRUE
          cat("  --- compile_validation_pdf first-failure diagnostics ---\n")
          cat(sprintf("  condition : %s\n", paste(class(e), collapse = "/")))
          cat(sprintf("  call      : %s\n", paste(deparse(conditionCall(e)), collapse = " ")))
          cat(sprintf("  devices   : cur=%d  open=%s\n", grDevices::dev.cur(),
                      paste(names(grDevices::dev.list()), collapse = ",")))
          cat(sprintf("  pdf so far: %s bytes at %s\n",
                      if (file.exists(pdf_path)) file.size(pdf_path) else "absent", pdf_path))
          .df <- tryCatch(system2("df", c("-h", shQuote(dirname(pdf_path))), stdout = TRUE), error = function(e2) "df unavailable")
          cat(sprintf("  disk      : %s\n", paste(.df, collapse = " || ")))
          .mp <- tryCatch(system2("memory_pressure", stdout = TRUE, stderr = TRUE), error = function(e2) "memory_pressure unavailable")
          cat(sprintf("  memory    : %s\n", paste(tail(.mp, 1), collapse = " ")))
          .w <- warnings()
          if (length(.w)) cat(sprintf("  warnings  : %s\n", paste(head(names(.w), 5), collapse = " | ")))
          cat("  --------------------------------------------------------\n")
        }
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
    # Force close any remaining devices (bounded -- see the note at the top)
    for (.i in seq_len(20)) { if (dev.cur() <= 1) break; try(dev.off(), silent = TRUE) }
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
      # Build OneDrive path with optional subfolder
      if (!is.null(external_subfolder)) {
        onedrive_dir <- file.path(config$paths$validation_plot_directory_onedrive, external_subfolder)
        dir.create(onedrive_dir, recursive = TRUE, showWarnings = FALSE)
        onedrive_pdf_path <- file.path(onedrive_dir, pdf_name)
      } else {
        onedrive_pdf_path <- file.path(config$paths$validation_plot_directory_onedrive, pdf_name)
      }
      
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