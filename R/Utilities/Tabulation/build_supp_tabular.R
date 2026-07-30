#' Build a supplementary table tabular (shared emitter for ST5 and ST6)
#'
#' Emits a self-contained LaTeX tabular (no \\begin{table} wrapper). Used by both
#' ST5 (covariate-chemical association screen) and ST6 (covariate-adjusted
#' tumor-type effects from \code{\link{adjust_type_effect}}); the two tables are
#' structurally identical, so they share one emitter rather than duplicating it.
#' Visual style matches \code{build_ST3()} / \code{build_ST4()}: section
#' headers (e.g. "Quantitative features") are plain underlined text at the base
#' indent, and individual chemicals are italic and indented
#' (\code{\\hspace*\{0.2cm\}}). Font (Times New Roman) is inherited from the
#' supplement at compile time.
#'
#' Two-line column headers are built with \code{\\shortstack} and vertically
#' centered with \code{\\raisebox}, matching ST1-ST4. Statistical symbols
#' (\emph{P}, \emph{n}) are italicized.
#'
#' @param data A data frame with a leading label column (\code{Chemical}) plus
#'   pre-formatted character value columns, and a logical \code{.section} column
#'   flagging mode-header rows. Value columns are emitted in the order given;
#'   supply them already rounded/formatted (this function does no numeric
#'   formatting, matching the contract of \code{build_ST4()}). Two-line headers
#'   are given as "line1|line2".
#'
#' @return A character scalar of LaTeX code ready to \\input.
#'
build_supp_tabular <- function(data) {
  df <- as.data.frame(data, check.names = FALSE, stringsAsFactors = FALSE)
  is_sec <- if (".section" %in% names(df)) as.logical(df$.section) else rep(FALSE, nrow(df))
  df$.section <- NULL
  label_col <- names(df)[1]
  value_cols <- setdiff(names(df), label_col)
  # LaTeX-escape cell content (% and & appear in detection/effect cells)
  esc <- function(x) {
    x <- ifelse(is.na(x), "", as.character(x))
    x <- gsub("%", "\\\\%", x)
    x <- gsub("&", "\\\\&", x)
    x
  }
  vmat <- as.matrix(df[value_cols])
  # Header cells: 10pt bold, vertically centered. Two-line labels are given as
  # "line1|line2" and become \shortstack; italicize the statistical symbols.
  ital <- function(s) {
    s <- gsub("\\bP\\b", "\\\\textit{P}", s)
    gsub("\\bn\\b", "\\\\textit{n}", s)
  }
  hdr_cell <- function(h) {
    parts <- strsplit(h, "|", fixed = TRUE)[[1]]
    inner <- if (length(parts) == 2L) {
      paste0("\\shortstack{", ital(parts[1]), "\\\\", ital(parts[2]), "}")
    } else {
      ital(h)
    }
    paste0("\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{", inner, "}}")
  }
  header_cells <- c(hdr_cell(label_col), vapply(value_cols, hdr_cell, character(1), USE.NAMES = FALSE))
  header_line <- paste0("\\rule{0pt}{12pt}", paste(header_cells, collapse = " & "), " \\\\")
  # Body: mode headers underlined at base indent; chemicals italic + indented
  body <- vapply(seq_len(nrow(df)), function(i) {
    lab <- esc(df[[label_col]][i])
    vals <- esc(vmat[i, ])
    if (is_sec[i]) {
      paste0("\\underline{", lab, "} & ", paste(vals, collapse = " & "))
    } else {
      paste0("\\hspace*{0.2cm}\\textit{", lab, "} & ", paste(vals, collapse = " & "))
    }
  }, character(1))
  colspec <- paste0("l", strrep("c", length(value_cols)))
  paste(c(
    paste0("\\begin{tabular}{", colspec, "}"),
    "\\toprule[0.5pt]",
    header_line,
    "\\midrule[0.5pt]",
    paste0(body, " \\\\"),
    "\\bottomrule[0.5pt]",
    "\\end{tabular}"
  ), collapse = "\n")
}
