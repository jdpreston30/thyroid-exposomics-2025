#' Plot Chemical Detection Frequency Distribution
#'
#' Creates a stacked bar plot showing the distribution of detected chemicals
#' across thyroid cancer variants. Bars are stacked by variant type with
#' standardized colors and formatting for publication.
#'
#' @param freq_dist_data A data frame with the following columns:
#'   \describe{
#'     \item{Bin}{Factor indicating detection count bins (e.g., "[260,270)")}
#'     \item{Follicular}{Integer count of Follicular variant samples}
#'     \item{FV-PTC}{Integer count of FV-PTC variant samples}
#'     \item{Papillary}{Integer count of Papillary variant samples}
#'   }
#'
#' @return A ggplot object with stacked bars showing sample counts by detection
#'   bin and variant type. Legend displays variants in order: Follicular, FV-PTC,
#'   Papillary. X-axis labels are rotated 45 degrees.
#'
#' @details
#' The plot includes:
#' \itemize{
#'   \item Y-axis range: 0-25 with breaks every 5 samples
#'   \item Variant-specific colors from the global variant_colors palette
#'   \item Black borders on bars and legend keys
#'   \item Arial font, size 14 bold for axis titles, size 12 bold for axis text
#' }
#'
#' @examples
#' \dontrun{
#' p2A <- plot_detection_distribution(freq_dist_bins)
#' ggsave("detection_distribution.pdf", p2A, width = 8, height = 6)
#' }
#'
#' @seealso \code{\link{variant_colors}} for color palette
#'
#' @export
plot_detection_distribution <- function(freq_dist_data) {
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  # Convert to long format for stacking
  plot_data <- freq_dist_data |>
    pivot_longer(
      cols = c(Follicular, `FV-PTC`, Papillary),
      names_to = "variant",
      values_to = "count"
    ) |>
    mutate(variant = factor(variant, levels = c("Follicular", "FV-PTC", "Papillary")))
  
  # Create stacked bar plot
  p <- ggplot(plot_data, aes(x = Bin, y = count, fill = variant)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.4, width = 0.6) +
    scale_fill_manual(
      values = variant_colors,
      breaks = c("Follicular", "FV-PTC", "Papillary")
    ) +
    guides(fill = guide_legend(
      keywidth = unit(0.25, "cm"),
      keyheight = unit(0.25, "cm"),
      byrow = FALSE,
      label.position = "right",
      ncol = 1
    )) +
    scale_y_continuous(
      limits = c(0, 25),
      breaks = seq(0, 25, by = 5),
      expand = expansion(mult = c(0, 0.05), add = 0)
    ) +
    labs(
      x = "Chemicals Detected",
      y = "Number of Samples",
      fill = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.position = c(0.02, 1),
      legend.justification = c(0, 1),
      legend.direction = "horizontal",
      legend.text = element_text(size = 8, face = "plain", family = "Arial"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "white", color = "black", linewidth = 0.25),
      legend.key.size = unit(0.25, "cm"),
      legend.spacing.x = unit(0.05, "cm"),
      legend.spacing.y = unit(0.02, "cm"),
      legend.box.spacing = unit(0, "cm"),
      axis.text.x = element_text(face = "bold", color = "black", size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(face = "bold", color = "black", size = 10),
      axis.title.x = element_text(face = "bold", color = "black", size = 12),
      axis.title.y = element_text(face = "bold", color = "black", size = 12, margin = margin(r = 10)),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8)
    )
  
  return(p)
}
