#' Global Validation Plot Adjustments
#'
#' Applies a series of global adjustments to validation plots in a single function:
#' - X-axis range adjustments (if validation_curated provided)
#' - Color standardization
#' - Y-axis title customization
#' - Grid line removal
#'
#' @param validation_plots Named list of plot objects
#' @param validation_curated Data frame with plot names and rt_adjustment column (optional)
#' @param y_title Character string for y-axis title (default: "← Standard | Tumor →")
#' @param remove_vertical_grid Logical, remove vertical grid lines (default: TRUE)
#' @param remove_horizontal_grid Logical, remove horizontal grid lines (default: FALSE)
#'
#' @return Named list of adjusted plots
#' @export
global_VP_adjust <- function(validation_plots,
                              validation_curated = NULL,
                              y_title = "← Standard | Tumor →",
                              remove_vertical_grid = TRUE,
                              remove_horizontal_grid = FALSE) {
  
  if (!is.list(validation_plots) || length(validation_plots) == 0) {
    stop("validation_plots must be a non-empty list")
  }
  
  cat(sprintf("\n🔧 Starting global adjustments for %d plots...\n", length(validation_plots)))
  
  adjusted_plots <- validation_plots
  
  # Step 1: Adjust x-axis ranges (if validation_curated provided)
  if (!is.null(validation_curated)) {
    cat("\n[1/4] Adjusting x-axis ranges...\n")
    adjusted_plots <- adjust_VP_ranges(
      validation_plots = adjusted_plots,
      validation_curated = validation_curated
    )
  } else {
    cat("\n[1/4] Skipping x-axis range adjustment (no validation_curated provided)\n")
  }
  
  # Step 2: Adjust colors
  cat("\n[2/4] Adjusting colors...\n")
  adjusted_plots <- adjust_VP_colors(
    validation_plots = adjusted_plots
  )
  
  # Step 3: Adjust y-axis title
  cat("\n[3/4] Adjusting y-axis title...\n")
  adjusted_plots <- adjust_VP_yaxis_title(
    validation_plots = adjusted_plots,
    y_title = y_title
  )
  
  # Step 4: Remove grid lines
  cat("\n[4/4] Adjusting grid lines...\n")
  adjusted_plots <- adjust_VP_grid(
    validation_plots = adjusted_plots,
    remove_vertical = remove_vertical_grid,
    remove_horizontal = remove_horizontal_grid
  )
  
  # Step 5: Force scientific notation on y-axis
  cat("\n[5/5] Forcing scientific notation on y-axis...\n")
  adjusted_plots <- lapply(adjusted_plots, function(p) {
    p + ggplot2::scale_y_continuous(labels = scales::scientific)
  })
  
  cat(sprintf("\n✅ Global adjustments complete for %d plots\n\n", length(adjusted_plots)))
  
  return(adjusted_plots)
}
