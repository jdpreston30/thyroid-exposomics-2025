#' Plot Qualitative Features Heatmap
#'
#' Creates a heatmap showing percent detection of qualitative features across
#' thyroid cancer variants. Uses hierarchical clustering to order compounds
#' by similarity.
#'
#' @param data A data frame in long format with columns:
#'   \describe{
#'     \item{variant}{Factor with levels: Follicular, FV-PTC, Papillary}
#'     \item{short_name}{Factor, chemical compound names (ordered by clustering)}
#'     \item{pct_detection}{Numeric, percent detection (0-100)}
#'   }
#'
#' @return A ggplot2 object showing a heatmap with:
#'   \itemize{
#'     \item Blue-to-red gradient (0% = #1A5FA5, 50% = #DDC9C2, 100% = #BF303B)
#'     \item Black borders between tiles (0.3 linewidth)
#'     \item Horizontal legend at top with "Percent Detection" title
#'     \item X-axis variant labels rotated 45 degrees
#'     \item Bold Arial text throughout
#'   }
#'
#' @examples
#' \dontrun{
#' p3B <- plot_qualitative_heatmap(qualitative_heatmap_data)
#' ggsave("qualitative_heatmap.pdf", p3B, width = 4, height = 5)
#' }
#'
#' @export
plot_qualitative_heatmap <- function(data) {
  library(ggplot2)
  
  p <- ggplot(data, aes(x = variant, y = short_name, fill = pct_detection)) +
    geom_tile(color = "black", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "#0C5EA5", mid = "#DDC9C2", high = "#BE4E4D", 
      midpoint = 50, limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      name = "Percent Detection",
      guide = guide_colorbar(
        title.position = "left",
        title.hjust = 0.5,
        barwidth = 0.8,
        barheight = 7,
        ticks.colour = "black",
        ticks.linewidth = 0.3,
        frame.colour = "black",
        frame.linewidth = 0.3,
        direction = "vertical"
      )
    ) +
    labs(x = NULL, y = NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), position = "right") +
    # coord_fixed(ratio = 0.5) +
    theme_minimal(base_size = 10, base_family = "Arial") +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.97, face = "bold", size = 10, color = "black", family = "Arial"),
      axis.text.y.right = element_text(hjust = 0, size = 9, face = "bold", color = "black", family = "Arial"),
      axis.text.y.left = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.1, "cm"),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10, color = "black", family = "Arial", angle = 90),
      legend.text = element_text(face = "bold", size = 8, color = "black", family = "Arial"),
      legend.background = element_rect(fill = "transparent", color = NA),
      plot.margin = margin(t = 5, r = 10, b = 10, l = 10)
    )
  
  return(p)
}
