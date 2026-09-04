#' Mark chemical locants for LaTeX italicisation
#'
#' `mark_locants()` wraps alphabetic locants in private-use sentinels;
#' `unmark_locants()` converts those sentinels into a LaTeX command once the text
#' is final. The two steps are separate because `gt::as_latex()` escapes
#' backslashes and braces, so raw LaTeX cannot be injected before `gt` runs --
#' the same reason the dagger markers are substituted post hoc in
#' `fix_latex_header_fill()`. Sentinels are ordinary characters and pass through
#' `gt` untouched.
#'
#' Only the letters are marked: hyphens, brackets and the commas separating
#' fusion locants stay upright, matching the manuscript.
#'
#' @param x Character vector of chemical names.
#' @param latex Character vector of LaTeX text containing sentinels.
#' @param cmd LaTeX command to apply, without the backslash. Use `"textit"` in
#'   upright surroundings and `"textup"` inside an italic run (ST6 row labels),
#'   where the convention reverses so the locant stays distinguishable.
#'
#' @return Character vector of the same length.
#'
#' @name mark_locants
NULL

#! U+E000/U+E001 are private-use: no font renders them, so an unconverted sentinel is a visible bug rather than a silent one.
.LOCANT_OPEN  <- "\uE000"
.LOCANT_CLOSE <- "\uE001"

#' @rdname mark_locants
mark_locants <- function(x) {
  o <- .LOCANT_OPEN; c_ <- .LOCANT_CLOSE
#! Fusion locants: each letter marked separately so the separating comma stays upright (dibenz[a,h]anthracene).
  x <- stringr::str_replace_all(x, "(?<=[A-Za-z])\\[([a-z](?:,[a-z0-9]+)*)\\](?=[a-z])", function(m) {
    vapply(m, function(one) {
      parts <- strsplit(stringr::str_sub(one, 2L, -2L), ",", fixed = TRUE)[[1]]
      parts <- ifelse(grepl("^[a-z]$", parts), paste0(o, parts, c_), parts)
      paste0("[", paste(parts, collapse = ","), "]")
    }, character(1), USE.NAMES = FALSE)
  })
#! N,N- and O,O-: both letters marked, the comma between them left upright.
  x <- stringr::str_replace_all(x, "(?<![A-Za-z])([NOS]),\\1-(?=[A-Za-z])",
                                paste0(o, "\\1", c_, ",", o, "\\1", c_, "-"))
#! Single locants and structural prefixes, at the start of a name or after a hyphen, space or parenthesis. Numeric locants (4-, 2,4-) and Greek ones (gamma-) are untouched: the alternation is letters only.
  x <- stringr::str_replace_all(
    x, "(?<=^|[-\\s(])(sec|tert|cis|trans|ortho|meta|para|[NOSomnpd])-(?=[A-Za-z])",
    paste0(o, "\\1", c_, "-"))
  x
}

#' @rdname mark_locants
unmark_locants <- function(latex, cmd = "textit") {
  stringr::str_replace_all(
    latex,
    paste0(.LOCANT_OPEN, "([^", .LOCANT_CLOSE, "]*)", .LOCANT_CLOSE),
    paste0("\\\\", cmd, "{\\1}"))
}
