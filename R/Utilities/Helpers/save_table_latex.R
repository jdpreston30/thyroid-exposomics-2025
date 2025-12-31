#' Save GT Table as LaTeX
#'
#' Converts a gt table to LaTeX and saves to the specified path.
#' Applies post-processing to ensure proper header fill.
#'
#' @param gt_obj A gt table object
#' @param file_path Full path where the .tex file should be saved
#' @param fix_header_fill Whether to apply LaTeX header fill fix (default: TRUE)
#'
#' @return Invisibly returns the file path
#'
save_table_latex <- function(gt_obj, file_path, fix_header_fill = TRUE) {
  
  # Ensure directory exists
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  
  # Convert to LaTeX
  latex_code <- gt::as_latex(gt_obj) |> as.character()
  
  # Apply header fill fix if requested
  if (fix_header_fill) {
    latex_code <- fix_latex_header_fill(latex_code)
  }
  
  # Save
  writeLines(latex_code, file_path)
  
  message("Table saved to: ", file_path)
  invisible(file_path)
}
