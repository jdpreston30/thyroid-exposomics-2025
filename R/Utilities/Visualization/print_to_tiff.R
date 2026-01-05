#' Print Plot to TIFF with High Resolution
#'
#' Saves ggplot objects as high-resolution TIFF files. TIFF format is ideal for
#' publication-quality figures that require lossless compression and archival stability.
#'
#' @param plot ggplot object or other plot object compatible with ggplot2::ggsave()
#' @param filename Character string for the TIFF filename. The ".tiff" extension
#'   will be automatically added if not present
#' @param width Numeric value for plot width in inches (default: 8.5 for letter size)
#' @param height Numeric value for plot height in inches (default: 11 for letter size)  
#' @param dpi Numeric value for resolution in dots per inch (default: 500 for high-quality output)
#' @param output_dir Character string specifying the output directory. If NULL (default),
#'   uses "Outputs/Figures/TIFF"
#' @param compression Character string specifying compression type. Options: "none", 
#'   "rle", "lzw", "jpeg", "zip" (default: "lzw" for good compression without quality loss)
#'
#' @return Invisibly returns the full file path to the created TIFF file
#'
#' @details
#' TIFF format is ideal for figures destined for publication or long-term archival.
#' Key features include:
#' 
#' - **Lossless compression**: LZW compression preserves all image data
#' - **High resolution**: Default 500 DPI suitable for print publication
#' - **Archival quality**: TIFF format is stable and widely supported
#' - **Path management**: Automatically handles file extensions and directory creation
#'
#' @examples
#' \dontrun{
#'   # Basic usage
#'   p <- ggplot(data, aes(x, y)) + geom_point()
#'   print_to_tiff(p, "figure1")
#'   
#'   # Custom resolution for higher quality
#'   print_to_tiff(p, "figure2", dpi = 600)
#'   
#'   # Specify custom output directory
#'   print_to_tiff(p, "figure3", output_dir = "Outputs/Figures/TIFF")
#' }
#'
#' @importFrom ggplot2 ggsave
#' @export
print_to_tiff <- function(plot, filename, width = 8.5, height = 11, dpi = 500,
                          output_dir = "Outputs/Figures/TIFF", compression = "lzw") {
  
  # Ensure filename doesn't have extension
  filename_base <- sub("\\.(tiff?|TIFF?)$", "", filename, ignore.case = TRUE)
  
  # Create TIFF filename
  tiff_filename <- paste0(filename_base, ".tiff")
  
  # Handle output directory
  if (is.null(output_dir)) {
    output_dir <- "Outputs/Figures/TIFF"
  }
  
  # Create directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Full path
  full_path <- file.path(output_dir, tiff_filename)
  
  # Save TIFF with specified parameters
  ggplot2::ggsave(
    filename = full_path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    device = "tiff",
    bg = "white",
    compression = compression
  )
  
  # Return full path invisibly
  invisible(full_path)
}
