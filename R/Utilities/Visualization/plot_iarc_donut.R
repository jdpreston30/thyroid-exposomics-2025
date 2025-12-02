#' Plot IARC Carcinogen Class Donut Chart
#'
#' Creates a donut chart showing the distribution of chemicals across IARC
#' carcinogen classification groups. The donut displays counts within each
#' segment, with "Not Classified" count shown in the center hole. Legend
#' appears to the right with total count below.
#'
#' @param iarc_data A data frame with the following columns:
#'   \describe{
#'     \item{IARC_Group}{Factor with levels: Group 1, Group 2A, Group 2B, Group 3, Not Classified}
#'     \item{n}{Integer count of chemicals in each group}
#'   }
#'
#' @return A ggplot object showing a donut chart with IARC group colors,
#'   counts within segments, Not Classified count in center, legend to right,
#'   and total count below legend.
#'
#' @details
#' The plot includes:
#' \itemize{
#'   \item Title: "IARC Carcinogen Class" (bold, 14pt) centered above donut
#'   \item Donut segments colored by IARC_colors palette
#'   \item Count labels ("n=X") within each segment
#'   \item Not Classified count in center hole (italic, not bold)
#'   \item Legend to right of donut
#'   \item Total count below legend ("Total Chemicals:\nn=X")
#'   \item Arial font throughout
#' }
#'
#' @examples
#' \dontrun{
#' p2B <- plot_iarc_donut(IARC)
#' ggsave("iarc_donut.pdf", p2B, width = 8, height = 6)
#' }
#'
#' @seealso \code{\link{IARC_colors}} for color palette
#'
#' @export
plot_iarc_donut <- function(iarc_data) {
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  
  # Separate Not Classified from other groups
  not_classified <- iarc_data |>
    filter(IARC_Group == "Not Classified") |>
    pull(n)
  
  # Prepare all data including Not Classified (but with minimal fraction)
  # This ensures all groups are treated equally in the legend
  plot_data <- iarc_data |>
    mutate(
      # Use actual fractions for Groups 1-3, tiny fraction for Not Classified
      fraction = ifelse(IARC_Group == "Not Classified", 0.00000000001, 
                       n / sum(iarc_data$n[iarc_data$IARC_Group != "Not Classified"])),
      ymax = cumsum(fraction),
      ymin = c(0, head(ymax, n = -1)),
      labelPosition = (ymax + ymin) / 2,
      label = paste0("n = ", n)
    )
  
  # Get only classified groups for text labels
  classified_data <- plot_data |>
    filter(IARC_Group != "Not Classified")
  
  # Calculate total
  total_n <- sum(iarc_data$n)
  
  # Create donut chart (inner=2.5, outer=6 for much thicker ring)
  # Plot all groups so Not Classified is treated as a real legend item
  donut <- ggplot(plot_data, aes(ymax = ymax, ymin = ymin, xmax = 6, xmin = 2.5, fill = IARC_Group)) +
    geom_rect(color = "black", linewidth = 0.3) +
    geom_text(
      data = classified_data,
      aes(x = 4.25, y = labelPosition, label = label),
      size = 8 / .pt,
      fontface = "bold.italic",
      family = "Arial"
    ) +
    # Add Not Classified count in center
    annotate(
      "text",
      x = 0,
      y = 0,
      label = paste0("n = ", not_classified),
      size = 8 / .pt,
      fontface = "italic",
      family = "Arial"
    ) +
    scale_fill_manual(
      values = IARC_colors,
      breaks = c("Group 1", "Group 2A", "Group 2B", "Group 3", "Not Classified"),
      labels = c("Group 1", "Group 2A", "Group 2B", "Group 3", "Not Classified"),
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
      label = "IARC Carcinogen Class",
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
