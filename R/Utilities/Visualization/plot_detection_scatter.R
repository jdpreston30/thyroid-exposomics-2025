#' Plot Chemical Detection Scatter Plot
#'
#' Creates a scatter plot showing the number of detected chemicals per sample,
#' grouped by thyroid cancer variant. Points are colored by variant type with
#' optional statistical test result displayed.
#'
#' @param detection_data A data frame with the following columns:
#'   \describe{
#'     \item{variant}{Character indicating the cancer variant (Follicular, FV-PTC, Papillary)}
#'     \item{total_detected}{Numeric count of total chemicals detected per sample}
#'   }
#' @param p_value Optional numeric p-value from statistical test (e.g., Kruskal-Wallis).
#'   If provided, will be displayed in the top-left corner of the plot.
#'
#' @return A ggplot object showing individual sample detection counts as points
#'   colored by variant with median lines. If p_value provided, displays formatted
#'   p-value text in top-left corner.
#'
#' @details
#' The plot includes:
#' \itemize{
#'   \item Y-axis range: 260-340 with automatic breaks
#'   \item X-axis: Variant groups
#'   \item Variant-specific colors from the global variant_colors palette
#'   \item Optional p-value annotation in top-left corner
#'   \item Arial font, size 12 bold for axis titles, size 10 bold for axis text
#' }
#'
#' @examples
#' \dontrun{
#' kw_p <- kruskal.test(total_detected ~ variant, data = detection_no_endog)$p.value
#' p1B <- plot_detection_scatter(detection_no_endog, p_value = kw_p)
#' ggsave("detection_scatter.pdf", p1B, width = 6, height = 6)
#' }
#'
#' @seealso \code{\link{variant_colors}} for color palette
#'
#' @export
plot_detection_scatter <- function(detection_data, p_value = NULL) {
  library(ggplot2)
  library(dplyr)
  
  # Ensure variant order
  plot_data <- detection_data |>
    mutate(variant = factor(variant, levels = c("Follicular", "FV-PTC", "Papillary")))
  
  # Calculate medians for each variant
  medians <- plot_data |>
    group_by(variant) |>
    summarise(median_value = median(total_detected, na.rm = TRUE), .groups = "drop")
  
  # Create scatter plot
  p <- ggplot(plot_data, aes(x = variant, y = total_detected, color = variant)) +
    geom_point(size = 2.5, alpha = 1, shape = 16, position = position_jitter(width = 0.35, seed = 42)) +
    geom_errorbar(
      data = medians,
      aes(x = variant, y = median_value, ymin = median_value, ymax = median_value, color = variant),
      width = 0.6,
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
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
      legend.position = "none",
      axis.text.x = element_text(face = "bold", color = "black", size = 10),
      axis.text.y = element_text(face = "bold", color = "black", size = 10),
      axis.title.y = element_text(face = "bold", color = "black", size = 12, margin = margin(r = 10)),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  # Add p-value annotation if provided
  if (!is.null(p_value)) {
    p_text <- if (p_value < 0.001) {
      "Kruskal-Wallis p < 0.001"
    } else {
      sprintf("Kruskal-Wallis p = %.3f", p_value)
    }
    
    p <- p + annotate(
      "text",
      x = -Inf,
      y = 339.5,
      label = p_text,
      hjust = -0.1,
      vjust = 1.5,
      size = 8 / .pt,
      family = "Arial",
      fontface = "plain"
    )
  }
  
  return(p)
}
