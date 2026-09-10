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
compile_sf_sub_pdf <- function(metadata, sf_sub_value, output_dir = here::here("Supplementary", "Components", "Figures", "PDF")) {
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
  output_path <- file.path(output_dir, paste0("S", sf_sub_value, ".pdf"))
  
  # Clean up any existing graphics devices before opening new PDF
  while (!is.null(dev.list())) {
    dev.off()
  }
  
#! cairo_pdf, not pdf(). The base pdf() device ignores the font faces gridtext sets, so every element_markdown title and subtitle came out plain: the italic locants (*o*-Toluidine, Benz[*a*]anthracene), the bold titles and the italic subtitles never rendered in any supplement build before 2026-09-10, although the PNGs (ragg) showed them. Its Type 1 Helvetica also cannot encode U+2032, so 2,4′-Methoxychlor printed with an acute accent, and its hyphens came out as U+2212. cairo_pdf embeds TrueType Helvetica with all four faces and renders the grobs exactly as the PNGs do.
#! family = "Arial": the figure convention (all figures Arial), and the sans face cairo picks for the grobs' default family lacks the U+2190/U+2192 arrows of the "← Standard | Tumor →" axis title (they rendered as boxes in a test); Arial carries them and the prime, and embeds all four faces.
  cairo_pdf(output_path, width = 8.5, height = 11, onefile = TRUE, family = "Arial")
  walk(pages, print)
  
  # Robust device cleanup
  while (!is.null(dev.list())) {
    dev.off()
  }
  
  message("Created: ", output_path)
}
