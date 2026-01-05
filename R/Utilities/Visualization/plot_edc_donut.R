#' Plot EDC Classification Donut Chart
#'
#' Creates a donut chart showing the distribution of chemicals across EDC
#' classification categories (Non-EDC and Potential EDC). The donut displays 
#' counts within each segment with legend to the right.
#'
#' @param edc_data A data frame with the following columns:
#'   \describe{
#'     \item{Potential_EDC}{Factor with levels: Non-EDC, Potential EDC}
#'     \item{n}{Integer count of chemicals in each category}
#'   }
#'
#' @return A ggplot object showing a donut chart with EDC colors,
#'   counts within segments, and legend to right.
#'
#' @details
#' The plot includes:
#' \itemize{
#'   \item Title: "EDC Classification" (bold, 12pt) centered above donut
#'   \item Donut segments colored by EDC_colors palette
#'   \item Count labels ("n=X") within each segment
#'   \item Legend to right of donut
#'   \item Arial font throughout
#' }
#'
#' @examples
#' \dontrun{
#' p2C <- plot_edc_donut(EDC)
#' ggsave("edc_donut.pdf", p2C, width = 8, height = 6)
#' }
#'
#' @seealso \code{\link{EDC_colors}} for color palette
#'
#' @export
plot_edc_donut <- function(edc_data) {
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  
  # Prepare plot data
  plot_data <- edc_data |>
    mutate(
      fraction = n / sum(n),
      ymax = cumsum(fraction),
      ymin = c(0, head(ymax, n = -1)),
      labelPosition = (ymax + ymin) / 2,
      label = paste0("n = ", n)
    )
  
  # Calculate total
  total_n <- sum(edc_data$n)
  
  # Create donut chart (inner=2.5, outer=6 for much thicker ring)
  donut <- ggplot(plot_data, aes(ymax = ymax, ymin = ymin, xmax = 6, xmin = 2.5, fill = Potential_EDC)) +
    geom_rect(color = "black", linewidth = 0.3) +
    geom_text(
      aes(x = 4.25, y = labelPosition, label = label),
      size = 8 / .pt,
      fontface = "bold.italic",
      family = "Arial"
    ) +
    scale_fill_manual(
      values = EDC_colors,
      breaks = c("Potential EDC", "Non-EDC"),
      labels = c("Potential EDC", "Non-EDC"),
      drop = FALSE
    ) +
    coord_polar(theta = "y") +
    xlim(c(0, 6)) +
    labs(fill = NULL) +
    theme_void(base_family = "Arial") +
    theme(
      legend.position = c(1.0, 0.5),
      legend.justification = c(0, 0.5),
      legend.text = element_text(size = 8, face = "plain", family = "Arial"),
      legend.key = element_rect(color = "black", linewidth = 0.2),
      legend.key.size = unit(0.4, "cm"),
      legend.key.spacing.y = unit(0.08, "cm"),
      legend.background = element_rect(fill = "transparent", color = NA),
      plot.margin = margin(t = 0, r = 10, b = 10, l = 10)
    ) +
    guides(fill = guide_legend(override.aes = list(linewidth = 0.2)))
  
  # Use donut with built-in legend (no plot_grid complications)
  combined <- donut
  
  # Add title using ggplot instead of draw_label
  title <- ggplot() +
    annotate(
      "text",
      x = 0.4,
      y = 0.5,
      label = "Endocrine Disrupting Chemicals",
      fontface = "bold",
      size = 10 / .pt,
      family = "Arial",
      hjust = 0.5,
      vjust = 0.5
    ) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(t = 5, r = 0, b = 5, l = 0))
  
  # Final plot with title
  final_plot <- plot_grid(
    title,
    combined,
    ncol = 1,
    rel_heights = c(0.08, 1)
  )
  
  return(final_plot)
}
