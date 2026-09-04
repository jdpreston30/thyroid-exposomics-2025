#' Render chemical locants italic in figure labels
#'
#' Converts a chemical name into a plotmath string (or expression) with its
#' alphabetic locant italicised, matching the manuscript convention.
#'
#' Plotmath rather than `ggtext::element_markdown()`: gridtext is a separate
#' layout engine that re-measures text, which shifts panel geometry and collapses
#' interior spaces ("Hexanoic acid" -> "Hexanoicacid"). Plotmath is grid's native
#' typesetting, so metrics are unchanged -- verified pixel-identical panel bounds.
#' This is the same mechanism already used for the `z`-score axis title.
#'
#' `face = "bold"` does NOT reach a plotmath label, so weight must be explicit;
#' every branch below emits `bold()` or `bolditalic()`.
#'
#' @param x Character vector of chemical names.
#'
#' @return `locant_plotmath()` returns plotmath source strings, for
#'   `ggplot2::as_labeller(..., default = label_parsed)`. `locant_expr()` returns
#'   an expression vector, for `scale_*_discrete(labels = )`.
#'
#' @name locant_plotmath
NULL

#' @rdname locant_plotmath
locant_plotmath <- function(x) {
#! Quoting is what keeps spaces, the U+2032 prime and the U+1D43 superscript intact -- plotmath only reinterprets unquoted tokens.
  .q <- function(s) paste0('"', gsub('"', '', s, fixed = TRUE), '"')
  vapply(x, function(one) {
    if (is.na(one) || !nzchar(one)) return('bold("")')
#! Leading locant: o-Cresol, N-MeFOSAA.
    m <- regmatches(one, regexec("^([NOSomnpd]|sec|tert|cis|trans|meta)-(.*)$", one))[[1]]
    if (length(m) == 3) return(sprintf("bolditalic(%s)*bold(%s)", .q(m[2]), .q(paste0("-", m[3]))))
#! Internal locant: Di-n-octyl phthalate, Tri-o-cresyl phosphate.
    m2 <- regmatches(one, regexec("^(.+)-([NOSomnpd]|sec|tert|cis|trans|meta)-(.*)$", one))[[1]]
    if (length(m2) == 4) return(sprintf("bold(%s)*bolditalic(%s)*bold(%s)",
                                        .q(paste0(m2[2], "-")), .q(m2[3]), .q(paste0("-", m2[4]))))
    sprintf("bold(%s)", .q(one))
  }, character(1), USE.NAMES = FALSE)
}

#' @rdname locant_plotmath
locant_expr <- function(x) {
  as.expression(lapply(locant_plotmath(x), function(s) parse(text = s)[[1]]))
}

#' @rdname locant_plotmath
#' @description `locant_title()` additionally pads the expression with symmetric phantom descenders.
#'   Plotmath reserves an expression's true bounding box, so a title without descenders
#'   ("o-Toluidine") reserves ~10 px less height at 300 dpi than one with them ("4-Aminobiphenyl"),
#'   which misaligns panel tops and title baselines between otherwise identical panels. Plain text
#'   never had this problem because it reserves a fixed line height regardless of glyphs. The
#'   phantoms are symmetric so horizontal centring is unaffected.
locant_title <- function(x) {
  parse(text = sprintf('paste(phantom("g"), %s, phantom("g"))', locant_plotmath(x)))[[1]]
}
