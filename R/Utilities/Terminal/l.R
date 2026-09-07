#' Quick Load Validation Plot(s)
#'
#' Loads one or more validation plot RDS files and assigns to global environment
#' with their plot names. Convenience function for manual inspection.
#'
#' @param ... Plot tags (unquoted), e.g., F1_S1_CP2486 or F1_S1_CP2486, F6_S1_CP3113
#'
#' @return Invisibly returns the plot object(s)
#' @export
#'
#' @examples
#' l(F1_S1_CP2486)  # Loads single plot
#' l(F1_S1_CP2486, F6_S1_CP3113, F1_S1_CP2107)  # Loads multiple plots
l <- function(...) {
  # Capture all unquoted names
  plot_names <- as.list(substitute(list(...)))[-1]
  plot_tags <- sapply(plot_names, deparse)
  
  if (length(plot_tags) == 0) {
    stop("No plot names provided. Usage: l(plot_name) or l(plot1, plot2, ...)")
  }
  
  base_path <- config$paths$validation_plot_directory_onedrive
  loaded_plots <- list()
  
#! curated/original no longer exists in either store, so this searched one dead path and returned an empty list; vp() then hard-stops at its "Invalid plot object" check. Script 09 line 107 is the only l() call in the pipeline and it aborted the run there. Search vp()'s own directories as a fallback -- l()'s purpose is to force a fresh load from disk rather than reuse a modified in-memory copy, which the variant_rtx original serves.
  search_dirs <- c(file.path("curated", "original"), "variant_rtx", "iarc_tumor_rtx", "iarc_cadaver_rtx")
  for (plot_tag in plot_tags) {
    # Construct path to RDS file
    rds_path <- NULL
    for (dir in search_dirs) {
      test_path <- file.path(base_path, dir, paste0(plot_tag, ".rds"))
      if (file.exists(test_path)) { rds_path <- test_path; break }
    }
    # Check if file exists
    if (is.null(rds_path)) {
      warning(sprintf("Plot file not found in any directory: %s", plot_tag))
      next
    }

    # Load the plot
    plot_obj <- readRDS(rds_path)
    
    # Assign to global environment with the plot tag name
    assign(plot_tag, plot_obj, envir = .GlobalEnv)
    loaded_plots[[plot_tag]] <- plot_obj
    
    cat(sprintf("✓ Loaded %s\n", plot_tag))
  }
  
  if (length(loaded_plots) == 1) {
    invisible(loaded_plots[[1]])
  } else {
    invisible(loaded_plots)
  }
}
