#' Build Supplementary Table 4 (ST4) — Candidate confounder assessment
#'
#' Emits a self-contained LaTeX tabular (no \\begin{table} wrapper) for the
#' covariate-confounding results (reviewer #2). Visual style mirrors
#' \code{build_ST3()}: section-header rows (e.g. "Covariate vs. tumor type") are
#' plain, underlined text at the base indent — like ST3's variable rows — and the
#' individual test rows are italic and indented (\\hspace*{0.2cm}), like ST3's
#' stratum rows. Font (Times New Roman) is inherited from the supplement at
#' compile time. The table is small and single-page, so it does not use the
#' longtable machinery.
#'
#' @param data A data frame with a leading label column (\code{Comparison}) plus
#'   any number of value columns (e.g. \code{Test}, \code{P}), all character and
#'   pre-formatted, and a logical \code{.section} column flagging section-header
#'   rows. On section rows the value columns are typically blank.
#'
#' @return A character scalar of LaTeX code ready to \\input.
#'
build_ST4 <- function(data) {
  df <- as.data.frame(data, check.names = FALSE, stringsAsFactors = FALSE)
  is_sec <- if (".section" %in% names(df)) as.logical(df$.section) else rep(FALSE, nrow(df))
  df$.section <- NULL
  label_col <- names(df)[1]
  value_cols <- setdiff(names(df), label_col)
  # LaTeX-escape cell content (% and & appear in test/label cells)
  esc <- function(x) {
    x <- ifelse(is.na(x), "", as.character(x))
    x <- gsub("%", "\\\\%", x)
    x <- gsub("&", "\\\\&", x)
    x
  }
  vmat <- as.matrix(df[value_cols])
  # Header cells: 10pt bold, vertically centered via \raisebox (matches ST3). The
  # statistical symbol P is italicized; all other headers are plain bold.
  hdr_cell <- function(h) {
    inner <- if (h == "P") "\\textit{P}" else h
    paste0("\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{", inner, "}}")
  }
  header_cells <- c(hdr_cell(label_col), vapply(value_cols, hdr_cell, character(1), USE.NAMES = FALSE))
  header_line <- paste0("\\rule{0pt}{12pt}", paste(header_cells, collapse = " & "), " \\\\")
  # Body: section rows underlined-plain at base indent; entry rows italic + indented
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
