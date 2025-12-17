#' Plot Chemical Class Distribution
#'
#' Creates a horizontal bar plot showing the number of chemicals annotated
#' in each use class. Classes are ordered from top to bottom by frequency,
#' with all bars in black.
#'
#' @param class_data A data frame with the following columns:
#'   \describe{
#'     \item{Graph_Class}{Character indicating the chemical use class}
#'     \item{n}{Integer count of chemicals in that class}
#'   }
#'
#' @return A ggplot object showing horizontal bars with chemical classes
#'   on the y-axis and counts on the x-axis. All bars are black with no fill color.
#'
#' @details
#' The plot includes:
#' \itemize{
#'   \item X-axis range: 0-100 with label "Number of Chemicals annotated"
#'   \item Y-axis: Chemical use classes ordered by frequency
#'   \item Black bars with no color differentiation
#'   \item Superscript symbols (†, ‡, §) in class labels
#'   \item Arial font, size 12 bold for axis titles, size 10 bold for axis text
#' }
#'
#' @examples
#' \dontrun{
#' p2A <- plot_class_distribution(all)
#' ggsave("class_distribution.pdf", p2A, width = 8, height = 6)
#' }
#'
#' @export
plot_class_distribution <- function(class_data, x_max = 100, sup = FALSE) {
  library(ggplot2)
  library(dplyr)
  
  # Process the data and order by frequency (reverse for bottom-to-top)
  plot_data <- class_data |>
    mutate(Graph_Class = factor(Graph_Class, levels = rev(Graph_Class)))
  
  # Calculate appropriate breaks based on x_max
  break_interval <- if (x_max <= 100) 20 else if (x_max <= 200) 50 else 50
  
  # Adjust font sizes for supplementary figures
  axis_text_x_size <- 10
  axis_text_y_size <- if (sup) 7 else 10
  axis_title_size <- if (sup) 10 else 12
  
  # Create horizontal bar plot
  p <- ggplot(plot_data, aes(x = n, y = Graph_Class)) +
    geom_col(fill = "black", color = "black", width = 0.6) +
    scale_x_continuous(
      limits = c(0, x_max),
      breaks = seq(0, x_max, by = break_interval),
      expand = expansion(mult = c(0, 0.05), add = 0)
    ) +
    labs(
      x = "Number of Chemicals Annotated",
      y = NULL
    ) +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      axis.text.x = element_text(face = "bold", color = "black", size = axis_text_x_size),
      axis.text.y = element_text(face = "bold", color = "black", size = axis_text_y_size, hjust = 1),
      axis.title.x = element_text(face = "bold", color = "black", size = axis_title_size, margin = margin(t = 10)),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  return(p)
}
