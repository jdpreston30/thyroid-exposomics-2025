#' Plot Qualitative Features Heatmap
#'
#' Creates a heatmap showing percent detection of qualitative features across
#' thyroid cancer variants. Uses hierarchical clustering to order compounds
#' by similarity.
#'
#' @param data A data frame in long format with columns:
#'   \describe{
#'     \item{variant}{Factor with levels: Follicular, IEFVPTC, Papillary}
#'     \item{short_name}{Factor, chemical compound names (ordered by clustering)}
#'     \item{pct_detection}{Numeric, percent detection (0-100)}
#'   }
#'
#' All five weight defaults below are the values approved for the published Figure 3, set by measuring
#' the render at 300 dpi rather than by computing them. Do not "correct" them with the px formula: it
#' predicts the nominal stroke, but \code{panel.border} and the colourbar frame are both clipped at their
#' element's edge, so part of the stroke is lost and the drawn result comes out lighter than the formula
#' says. That is why the outer border needs 1.55 to read as the same weight as the 0.8 frames elsewhere
#' in the figure. Change one of these only against a measured render.
#'
#' @param border_linewidth Outer panel border weight. Default 1.55, tuned so the drawn border matches the
#'   frames/axis lines of panels 3A, 3D, 3E and 3F.
#' @param tile_linewidth Weight of the lines between heatmap cells. Default 0.39, deliberately lighter
#'   than the outer border -- same heavy-outer / light-inner relationship 3A has.
#' @param tick_linewidth Axis tick weight. Default 0.8, matching the ticks in 3A/3D/3E/3F. Ticks are not
#'   clipped, so here the px formula does hold: 0.8 renders 7.1 px at 300 dpi.
#' @param tick_length_px Tick length in pixels at 300 dpi. Default 14. Converted to inches internally so
#'   the physical length holds at the 800 dpi export.
#' @param legend_frame_linewidth Weight of the box around the colourbar. Default 0.8. Separate from
#'   \code{border_linewidth} because the colourbar frame and the panel border do not render a given value
#'   at the same drawn thickness. The colourbar has no tick marks -- see the note at \code{legend.ticks}.
#'
#' Nominal conversion (ticks only): \code{linewidth = px / 300 * 96 / 2.845276}. So 7.1 px = 0.8.
#'
#' @return A ggplot2 object showing a heatmap with:
#'   \itemize{
#'     \item Blue-to-red gradient (0% = #1A5FA5, 50% = #DDC9C2, 100% = #BF303B)
#'     \item Black borders between tiles (\code{tile_linewidth}, default 0.3)
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
plot_qualitative_heatmap <- function(data,
                                     border_linewidth = 1.55,
                                     tile_linewidth   = 0.39,
                                     tick_linewidth   = 0.8,
                                     tick_length_px   = 14,
                                     legend_frame_linewidth = 0.8) {
  library(ggplot2)
  
  p <- ggplot(data, aes(x = variant, y = short_name, fill = pct_detection)) +
    geom_tile(color = "black", linewidth = tile_linewidth) +
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
#! No colourbar ticks: they are drawn inward from both edges, so the 0 and 100 marks land on the bar's
#! ends and get clipped by the frame. There is no expand option -- the bar spans the scale limits by
#! construction -- and widening `limits` to make room would change the value-to-colour mapping for every
#! value except the midpoint, heatmap tiles included (0 shifts #0C5EA5 -> #3667A8). Dropping the marks
#! keeps all five labels and the mapping. Must be element_blank(): `ticks = FALSE` is deprecated as of
#! ggplot2 3.5 and 4.0.1 silently discards it into `...`, leaving the ticks at their DEFAULT WHITE.
        theme = theme(legend.ticks = element_blank()),
#! frame.colour/frame.linewidth are also deprecated but DO still apply in 4.0.1 (verified by diffing
#! renders). The frame draws ~2x thicker than the same value elsewhere, hence its own parameter. If they
#! go defunct, fold into the theme above: legend.frame = element_rect(linewidth = legend_frame_linewidth)
        frame.colour = "black",
        frame.linewidth = legend_frame_linewidth,
        direction = "vertical"
      )
    ) +
    labs(x = NULL, y = NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    # coord_fixed(ratio = 0.5) +
    theme_minimal(base_size = 10, base_family = "Arial") +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.97, face = "bold", size = 10, color = "black", family = "Arial"),
      axis.text.y = element_text(hjust = 1, size = 9, face = "bold", color = "black", family = "Arial"),
#! Defaults match 3A/3D/3E/3F at 7.1 px. Thin lines antialias wider than nominal, so a 0.3 tick
#! computes to 2.7 px but measures nearer 4.
      axis.ticks = element_line(color = "black", linewidth = tick_linewidth),
#! Inches, so the length is exact at 300 dpi and holds physically at the 800 dpi export. Shorter than
#! the 0.15 cm (17.7 px) used by 3D/3E/3F -- this panel is narrow and long ticks crowd the 45-deg labels.
      axis.ticks.length = unit(tick_length_px/300, "in"),
      axis.line = element_blank(),
      panel.grid = element_blank(),
#! Outer border, independent of the tile dividers -- same heavy-outer / light-inner relationship 3A has.
      panel.border = element_rect(color = "black", fill = NA, linewidth = border_linewidth),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10, color = "black", family = "Arial", angle = 90),
      legend.text = element_text(face = "bold", size = 8, color = "black", family = "Arial"),
      legend.background = element_rect(fill = "transparent", color = NA),
      plot.margin = margin(t = 5, r = 10, b = 10, l = 10)
    )
  
  return(p)
}
