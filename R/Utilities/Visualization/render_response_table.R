#' Render a Small Table to PNG for Pasting into a Reviewer Response
#'
#' Writes a clean, self-sizing PNG of a small data frame for inclusion in a
#' reviewer-response document. Intended for exhibits that support the response
#' but are not part of the published supplement, so they never enter the
#' manuscript's table numbering.
#'
#' Uses \code{gridExtra::tableGrob} rendered through \code{ragg::agg_png} rather
#' than \code{gt::gtsave()}, which requires a headless browser (webshot2 /
#' chromote) that is not installed in this project.
#'
#' Column headers may contain \code{\\n} to wrap onto two lines. The default font is
#' Times New Roman so the exhibit matches response-letter body text; this is
#' deliberately unlike the project's figures, which use Arial, because these tables
#' are pasted into a prose document rather than the manuscript.
#'
#' @param data Data frame to render. All columns are drawn as-is, so format
#'   numbers to strings beforehand.
#' @param path Output PNG path. Parent directory is created if needed.
#' @param title Optional bold title drawn above the table.
#' @param footnote Optional footnote drawn below the table, in smaller italic
#'   type. May contain \code{\\n} for multiple lines.
#' @param font_family Font used throughout (default "Times New Roman", matching
#'   response-letter body text).
#' @param dpi Resolution (default 300, suitable for print).
#' @param base_size Base font size in points (default 9).
#'
#' @return Invisibly, the path written.
#'
#' @examples
#' \dontrun{
#' render_response_table(
#'   covariate_screen_table,
#'   "Outputs/Revisions/covariate_screen_table.png",
#'   title = "Association of covariates with chemical features",
#'   footnote = "Approximately 5% of features are expected to reach P < 0.05 by chance."
#' )
#' }
#'
#' @export
render_response_table <- function(data, path, title = NULL, footnote = NULL,
                                  font_family = "Times New Roman", dpi = 300, base_size = 9) {
  stopifnot(is.data.frame(data), nrow(data) > 0)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  df <- as.data.frame(lapply(data, as.character), check.names = FALSE, stringsAsFactors = FALSE)
  names(df) <- names(data)
  tt <- gridExtra::ttheme_minimal(
    base_size = base_size,
    core = list(
      fg_params = list(fontfamily = font_family, hjust = 0, x = 0.04),
      bg_params = list(fill = c("white", "grey97"), col = NA)
    ),
    colhead = list(
      fg_params = list(fontfamily = font_family, fontface = "bold", hjust = 0, x = 0.04),
      bg_params = list(fill = "white", col = NA)
    )
  )
  tg <- gridExtra::tableGrob(df, rows = NULL, theme = tt)
  # Rules above and below the header, and one closing the table (booktabs feel)
  hline <- function(y, lwd) {
    grid::segmentsGrob(x0 = grid::unit(0, "npc"), x1 = grid::unit(1, "npc"),
                       y0 = grid::unit(y, "npc"), y1 = grid::unit(y, "npc"),
                       gp = grid::gpar(lwd = lwd))
  }
  tg <- gtable::gtable_add_grob(tg, hline(1, 1.4), t = 1, b = 1, l = 1, r = ncol(tg))
  tg <- gtable::gtable_add_grob(tg, hline(0, 0.8), t = 1, b = 1, l = 1, r = ncol(tg))
  tg <- gtable::gtable_add_grob(tg, hline(0, 1.4), t = nrow(tg), b = nrow(tg), l = 1, r = ncol(tg))
  #! Heights are accumulated per-part and summed, so the canvas is derived from the actual
  #! grobs rather than padding constants; a multi-line footnote would otherwise be clipped.
  parts <- list()
  hts <- list()
  if (!is.null(title)) {
    tgrob <- grid::textGrob(title, x = 0.02, hjust = 0,
      gp = grid::gpar(fontfamily = font_family, fontface = "bold", fontsize = base_size + 1))
    parts <- c(parts, list(tgrob))
    hts <- c(hts, list(grid::grobHeight(tgrob) + grid::unit(5, "mm")))
  }
  parts <- c(parts, list(tg))
  hts <- c(hts, list(sum(tg$heights)))
  w_tbl <- grid::convertWidth(sum(tg$widths), "in", valueOnly = TRUE)
  if (!is.null(footnote)) {
    #! Wrap over-long lines to the table width; existing \n breaks are preserved, so a
    #! footnote wider than the table neither clips nor forces a banner-wide canvas.
    fgp <- grid::gpar(fontfamily = font_family, fontface = "italic", fontsize = base_size - 1.5)
    footnote <- paste(unlist(lapply(strsplit(footnote, "\n", fixed = TRUE)[[1]], function(ln) {
      if (!nzchar(ln)) return(ln)
      w <- grid::convertWidth(grid::grobWidth(grid::textGrob(ln, gp = fgp)), "in", valueOnly = TRUE)
      if (w <= w_tbl) return(ln)
      strwrap(ln, width = max(20L, floor(nchar(ln) / w * w_tbl)))
    })), collapse = "\n")
    fgrob <- grid::textGrob(footnote, x = 0.02, hjust = 0, vjust = 1, gp = fgp)
    parts <- c(parts, list(fgrob))
    hts <- c(hts, list(grid::grobHeight(fgrob) + grid::unit(6, "mm")))
  }
  heights <- do.call(grid::unit.c, hts)
  combined <- gridExtra::arrangeGrob(grobs = parts, ncol = 1, heights = heights)
  # Canvas takes the widest element, not just the table, so nothing can be clipped
  w_parts <- vapply(parts, function(g) {
    if (inherits(g, "gtable")) grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE)
    else grid::convertWidth(grid::grobWidth(g), "in", valueOnly = TRUE)
  }, numeric(1))
  w_in <- max(w_parts) + 0.3
  h_in <- grid::convertHeight(sum(heights), "in", valueOnly = TRUE) + 0.2
  ragg::agg_png(path, width = w_in, height = h_in, units = "in", res = dpi, background = "white")
  grid::grid.newpage()
  grid::grid.draw(combined)
  grDevices::dev.off()
  message("Created: ", path)
  invisible(path)
}
