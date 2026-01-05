#' Print Plot to SVG with Auto-Refresh for macOS Preview
#'
#' Saves ggplot objects as scalable vector graphics (SVG) files with automatic 
#' Preview.app integration on macOS. The function handles file path management, 
#' automatic directory creation, and provides seamless workflow integration with 
#' Preview's auto-refresh functionality for iterative plot development.
#'
#' @param plot ggplot object or other plot object compatible with ggplot2::ggsave()
#' @param filename Character string for the SVG filename. The ".svg" extension
#'   will be automatically added if not present
#' @param width Numeric value for plot width in inches (default: 8.5 for letter size)
#' @param height Numeric value for plot height in inches (default: 11 for letter size)  
#' @param output_dir Character string specifying the output directory. If NULL (default),
#'   uses \code{config$paths$figures} or falls back to a default path
#' @param auto_open Logical indicating whether to automatically open the SVG in
#'   Preview.app on first save (default: TRUE). Subsequent saves will auto-refresh
#'   if Preview is already open
#'
#' @return Invisibly returns the full file path to the created SVG file
#'
#' @details
#' This function provides a streamlined workflow for saving and viewing plots during
#' analysis. Key features include:
#' 
#' - **Vector format**: Creates resolution-independent SVG files perfect for PDF output
#' - **Auto-refresh workflow**: Opens SVG in Preview on first save, subsequent saves
#'   automatically refresh the Preview window
#' - **Path management**: Automatically handles file extensions and directory creation
#' - **Config integration**: Uses dynamic configuration paths when available
#' - **High-quality output**: Vector graphics maintain quality at any scale
#' - **Background handling**: Ensures white background for clean output
#'
#' SVG files are ideal for:
#' - Scientific publications (vector graphics scale perfectly)
#' - PDF generation (maintains text as editable text, not pixels)
#' - Presentations (crisp at any zoom level)
#' - Further editing in Illustrator/Inkscape
#'
#' @examples
#' \dontrun{
#'   # Basic usage
#'   p <- ggplot(data, aes(x, y)) + geom_point()
#'   print_to_svg(p, "scatter_plot")
#'   
#'   # Custom dimensions for wide figure
#'   print_to_svg(p, "wide_plot", width = 12, height = 6)
#'   
#'   # Specify custom output directory
#'   print_to_svg(p, "analysis_figure", output_dir = "~/Desktop/figures")
#' }
#'
#' @importFrom ggplot2 ggsave
#' @export
print_to_svg <- function(plot, filename, width = 8.5, height = 11,
                         output_dir = NULL, auto_open = TRUE) {
  
  # Use config path if output_dir not specified
  if (is.null(output_dir)) {
    if (exists("config") && !is.null(config$paths$figures)) {
      output_dir <- config$paths$figures
    } else {
      # Fallback if config not loaded
      output_dir <- "/Users/jdp2019/Desktop/PGD_figures"
      warning("Config not found, using fallback path: ", output_dir)
    }
  }
  
  # Ensure filename has .svg extension
  if (!grepl("\\.svg$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".svg")
  }

  # Create full path
  filepath <- file.path(output_dir, filename)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Check if file already exists (for auto-open logic)
  file_exists <- file.exists(filepath)

  # Save the plot as SVG
  ggplot2::ggsave(
    filename = filepath,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = "svg",
    bg = "white"
  )

  # Auto-open in Preview only on first run (or if specified)
  if (auto_open && !file_exists) {
    system(paste("open", shQuote(filepath)))
    cat("SVG saved and opened in Preview:", filepath, "\n")
    cat("Preview will auto-refresh when you re-run this function!\n")
  } else {
    cat("SVG updated:", filepath, "\n")
  }

  # Return path invisibly
  invisible(filepath)
}
