#' Plot Chemical Detection Scatter Plot
#'
#' Creates a scatter plot showing the number of detected chemicals per sample,
#' grouped by thyroid cancer variant. Points are colored by variant type with
#' a horizontal legend positioned at the top of the plot area.
#'
#' @param detection_data A data frame with the following columns:
#'   \describe{
#'     \item{variant}{Character indicating the cancer variant (Follicular, FV-PTC, Papillary)}
#'     \item{total_detected}{Numeric count of total chemicals detected per sample}
#'   }
#'
#' @return A ggplot object showing individual sample detection counts as points
#'   colored by variant. Legend displays horizontally at the top (y = 339) in order:
#'   Follicular, FV-PTC, Papillary.
#'
#' @details
#' The plot includes:
#' \itemize{
#'   \item Y-axis range: 260-340 with automatic breaks
#'   \item X-axis: Variant groups
#'   \item Variant-specific colors from the global variant_colors palette
#'   \item Horizontal legend centered at top of plot
#'   \item Arial font, size 14 bold for axis titles, size 12 bold for axis text
#' }
#'
#' @examples
#' \dontrun{
#' p1B <- plot_detection_scatter(detection_no_endog)
#' ggsave("detection_scatter.pdf", p1B, width = 6, height = 6)
#' }
#'
#' @seealso \code{\link{variant_colors}} for color palette
#'
#' @export
plot_detection_scatter <- function(detection_data) {
  library(ggplot2)
  library(dplyr)
  
  # Ensure variant order
  plot_data <- detection_data %>%
    mutate(variant = factor(variant, levels = c("Follicular", "FV-PTC", "Papillary")))
  
  # Create scatter plot
  p <- ggplot(plot_data, aes(x = variant, y = total_detected, color = variant)) +
    geom_point(size = 3, alpha = 0.7, position = position_jitter(width = 0.2, seed = 42)) +
    scale_color_manual(
      values = variant_colors,
      breaks = c("Follicular", "FV-PTC", "Papillary")
    ) +
    scale_y_continuous(
      limits = c(260, 340),
      expand = expansion(mult = c(0, 0.05), add = 0)
    ) +
    labs(
      x = NULL,
      y = "# of Chemicals Detected",
      color = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.position = c(0.5, 1),
      legend.justification = c(0.5, 0),
      legend.direction = "horizontal",
      legend.text = element_text(size = 12, face = "plain", family = "Arial"),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.spacing.x = unit(0.2, "cm"),
      axis.text.x = element_text(face = "bold", color = "black", size = 12),
      axis.text.y = element_text(face = "bold", color = "black", size = 12),
      axis.title.y = element_text(face = "bold", color = "black", size = 14, margin = margin(r = 10)),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      plot.margin = margin(t = 30, r = 10, b = 10, l = 10)
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 4, alpha = 1)
    ))
  
  return(p)
}
