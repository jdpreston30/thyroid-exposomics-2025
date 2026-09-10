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
#! Two hard-won constraints shape this function (2026-09-10; punch list 133 and the compile memory finding):
#! 1. Nothing is written inside the repo while a device is open. The repo lives in ~/Desktop, an iCloud Drive file-provider domain; iCloud swapped the target's inode while pdf() held it open and the next page write failed. All rendering goes to R's tempdir and the finished PDF is copied in with one file.copy().
#! 2. Every compound renders in a FORKED CHILD. Under ggplot2 4.0 a ggplot deserialised from RDS retains ~165 MB in the parent's heap every time ggplotGrob() builds it (reachable by gc, not released by rm/gc/dev.off; a freshly constructed plot does not do this). 176 grobs built in one process put 51 GB into swap and made one compile take 51 minutes. A child builds its compound's pages, writes one part file and exits; the parent's heap stays flat (measured: +1 MB over 6 builds). The parts are merged with pdfunite.
  parts_dir <- tempfile("compile_parts_")
  dir.create(parts_dir)
  on.exit(unlink(parts_dir, recursive = TRUE), add = TRUE)

#! Reads one grob back from disk. rtx() no longer returns ggplot objects (punch list 134): each entry carries rds_name/rds_folder, and by the time compile runs the RDS has been transferred to validation_plot_directory_onedrive and the local copy deleted, so look there first.
  .load_plot <- function(plot_info) {
    p <- plot_info[["plot"]]
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
#! Same defect as remove_standard()/zoom_fragment() had: mean(rt_range) is the midpoint of two %.2f-rounded endpoints, not the measured RT. Diagnostic PDFs only, but there is no reason to print a different number than the supplement.
      new_subtitle <- sprintf("Sample: %s  |  Standard: %s  |  RT = %.3f min  |  %s",
                              plot_info$sample_id, plot_info$standard_file, subtitle_rt_of(plot_info), plot_info$plot_tag)
      p <- p +
        ggplot2::labs(subtitle = new_subtitle) +
        ggplot2::theme(
          plot.subtitle = ggtext::element_markdown(hjust = 0.5, face = "italic", size = 6,
                                                  color = "black", lineheight = 1.2, margin = ggplot2::margin(0, 0, 3, 0))
        )
    }
    p
  }

#! Renders one compound (all its pages, 3x2 per page) to part_path. Runs inside the fork; returns the page count, or the error message as a string.
  .render_compound <- function(compound, part_path) {
    plots <- compound$plots
    plot_labels <- names(plots)
    plot_labels <- plot_labels[order(
      as.numeric(gsub("F([0-9]+)_S([0-9]+)", "\\1", plot_labels)),
      as.numeric(gsub("F([0-9]+)_S([0-9]+)", "\\2", plot_labels))
    )]
    pdf_plots <- lapply(plot_labels, function(label) .load_plot(plots[[label]]))
    plots_per_page <- 6
    n_pages <- ceiling(length(pdf_plots) / plots_per_page)
    pdf(part_path, width = 8.5, height = 11, family = "Helvetica", onefile = TRUE)
    on.exit(try(dev.off(), silent = TRUE), add = TRUE)
    for (page_num in seq_len(n_pages)) {
      start_idx <- (page_num - 1) * plots_per_page + 1
      end_idx <- min(page_num * plots_per_page, length(pdf_plots))
      page_plots <- pdf_plots[start_idx:end_idx]
      if (length(page_plots) < plots_per_page) {
        blank_plot <- ggplot2::ggplot() + ggplot2::theme_void()
        for (i in seq_len(plots_per_page - length(page_plots))) page_plots[[length(page_plots) + 1]] <- blank_plot
      }
      title_text <- if (n_pages == 1) sprintf("%s (%s)", compound$short_name, compound$id)
                    else sprintf("%s (%s) - page %d/%d", compound$short_name, compound$id, page_num, n_pages)
      title_grob <- grid::textGrob(title_text, gp = grid::gpar(fontsize = 14, fontface = "bold"), x = 0.5, y = 0.95, just = "top")
      grid_plot <- gridExtra::arrangeGrob(grobs = page_plots, ncol = 2, top = title_grob)
      if (page_num > 1) grid::grid.newpage()
      grid::grid.draw(grid_plot)
    }
    n_pages
  }

  part_files <- character()
  n_failed <- 0L
  for (compound_id in names(compound_plots)) {
    compound <- compound_plots[[compound_id]]
    if (length(compound$plots) == 0) next
    cat(sprintf("  Adding %s (%s): %d plots\n", compound$short_name, compound$id, length(compound$plots)))
    part_path <- file.path(parts_dir, sprintf("%03d_%s.pdf", length(part_files) + 1L, compound$id))
    job <- parallel::mcparallel(tryCatch(.render_compound(compound, part_path), error = function(e) conditionMessage(e)))
    res <- parallel::mccollect(job)[[1]]
    if (is.numeric(res) && file.exists(part_path) && file.size(part_path) > 0) {
      part_files <- c(part_files, part_path)
    } else {
      n_failed <- n_failed + 1L
      cat(sprintf("  Error rendering compound %s: %s\n", compound$id, if (is.character(res)) res else "child returned no PDF"))
    }
  }
  if (length(part_files) == 0) {
    warning("PDF may not have been created successfully: every compound failed to render")
    return(invisible(NULL))
  }

  pdf_tmp <- tempfile("compile_", fileext = ".pdf")
  merged <- if (length(part_files) == 1) file.copy(part_files, pdf_tmp) else
    system2("pdfunite", c(shQuote(part_files), shQuote(pdf_tmp)), stdout = FALSE, stderr = FALSE) == 0
  if (!isTRUE(merged) || !file.exists(pdf_tmp)) {
    warning("PDF may not have been created successfully: pdfunite failed to merge the compound parts")
    return(invisible(NULL))
  }
  if (!file.copy(pdf_tmp, pdf_path, overwrite = TRUE)) {
    warning(sprintf("PDF may not have been created successfully: could not copy %s to %s", pdf_tmp, pdf_path))
    return(invisible(NULL))
  }
  unlink(pdf_tmp)
  cat(sprintf("\n✓ PDF saved to: %s (%d compounds%s)\n", pdf_path, length(part_files),
              if (n_failed > 0) sprintf(", %d FAILED to render", n_failed) else ""))

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
