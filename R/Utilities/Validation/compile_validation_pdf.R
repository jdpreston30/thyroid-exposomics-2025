#' Compile validation plots into supplementary figure subsection PDFs
#'
#' Creates multi-page PDFs for supplementary figure subsections with validation plots
#' arranged two per page (top and bottom). Each sf_sub value (e.g., "2.1", "2.2") gets
#' its own PDF file. Fragment plots (_F variants) are paired with their parent plots on
#' the same page. Blank placeholders are inserted when needed to maintain layout.
#'
#' @param metadata Data frame with validation plot metadata containing columns:
#'   \itemize{
#'     \item sf_sub: Subsection identifier (e.g., "2.1", "3.2")
#'     \item panel: Position on page ("top" or "bottom")
#'     \item grob: Plot object (ggplot grob)
#'   }
#' @param sf_sub_value Character sf_sub identifier to compile (e.g., "2.1" for SF2.1)
#' @param output_dir Character path to output directory for PDFs.
#'   Default: here::here("Outputs", "Figures")
#'
#' @return Invisibly returns NULL. Creates PDF file as side effect and prints
#'   confirmation message with output path.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Filters metadata to specified sf_sub value
#'   \item Separates plots by panel position ("top" or "bottom")
#'   \item Creates pages pairing top and bottom plots in order
#'   \item Inserts blank plots when a pair has only one plot
#'   \item Saves to PDF with dimensions 8.5 x 11 inches (US letter)
#'   \item Names output file as SF{sf_sub}.pdf (e.g., SF2.1.pdf, SF3.2.pdf)
#' }
#'
#' @examples
#' \dontrun{
#' compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "2.1")
#' compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "3.2")
#' }
#'
#' @export
compile_sf_sub_pdf <- function(metadata, sf_sub_value, output_dir = here::here("Outputs", "Figures")) {
  # Create blank plot for padding
  blank_plot <- ggplot() + theme_void()
  
  sf_data <- metadata %>% filter(sf_sub == sf_sub_value)
  
  # Get all plot pairs (rows with panel 'top' and 'bottom')
  top_plots <- sf_data %>% filter(panel == "top") %>% pull(grob)
  bottom_plots <- sf_data %>% filter(panel == "bottom") %>% pull(grob)
  
  # Determine number of pages needed
  n_pages <- max(length(top_plots), length(bottom_plots))
  
  # Create pages
  pages <- map(seq_len(n_pages), function(i) {
    # Get top and bottom plots for this page
    top_plot <- if (i <= length(top_plots)) top_plots[[i]] else blank_plot
    bottom_plot <- if (i <= length(bottom_plots)) bottom_plots[[i]] else blank_plot
    
    # Create page with top and bottom plots
    cowplot::plot_grid(
      top_plot,
      bottom_plot,
      ncol = 1,
      nrow = 2
    )
  })
  
  # Save to PDF
  output_path <- file.path(output_dir, paste0("SF", sf_sub_value, ".pdf"))
  pdf(output_path, width = 8.5, height = 11)
  walk(pages, print)
  dev.off()
  
  message("Created: ", output_path)
}
