#' Draw a single frame around the whole panel block, independent of panel.border
#'
#' \code{panel.border} is per-panel and applies to all four edges, so in a faceted plot with
#' \code{panel.spacing = 0} the internal dividers are the same element as the outer edge and cannot be
#' weighted separately. This adds one rectangle spanning every panel in the gtable, letting
#' \code{panel.border} carry the divider weight while the outer frame gets its own.
#'
#' Extents come from the gtable layout rather than hardcoded coordinates, so the frame stays flush with
#' the panel block regardless of strip text, axis label or legend size.
#'
#' @param p A ggplot object.
#' @param linewidth Frame weight in ggplot \code{linewidth} units (converted to grid \code{lwd}).
#' @param color Frame colour.
#' @return A gtable, which \code{cowplot::draw_plot()} and \code{ggplot2::ggsave()} both accept.
#' @export
add_panel_frame <- function(p, linewidth = 0.8, color = "black") {
  g <- ggplot2::ggplotGrob(p)
#+ 1: Locate every panel cell in the layout
  pn <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
  if (nrow(pn) == 0) return(g)
#+ 2: One rect spanning the full panel block
#! linewidth -> lwd is x2.845276 (ggplot2's .pt); grid lwd 1 = 1/96 inch.
  frame <- grid::rectGrob(gp = grid::gpar(col = color, fill = NA, lwd = linewidth * 2.845276))
  gtable::gtable_add_grob(
    g, frame,
    t = min(pn$t), b = max(pn$b), l = min(pn$l), r = max(pn$r),
    name = "outer-panel-frame", clip = "off"
  )
}
