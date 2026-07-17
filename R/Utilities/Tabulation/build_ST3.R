#' Build Supplementary Table 3 (ST3) — Cohort demographics
#'
#' Emits a self-contained LaTeX tabular (no \\begin{table} wrapper) for the
#' DTC-vs-cadaver demographic comparison produced by \code{TernTables::ternG()}.
#' Portrait, 0.5pt booktabs rules, headers vertically centered via
#' \code{\\raisebox} (matching ST1/ST2). Font (Times New Roman) is inherited from
#' the supplement at compile time. This table is small and single-page, so it
#' does not use the longtable machinery in \code{fix_latex_header_fill()}.
#'
#' Column-1 styling follows a shifted TernTables hierarchy (no top-level
#' categories in this table):
#'   - Variable rows (Sex, Age, Sample Collection Timing): plain text at the
#'     GROUP indent level (flush left); underlined only if the variable has
#'     sub-strata (i.e. categorical).
#'   - Strata rows (Female/Male, collection bins): italic at the Class indent
#'     level (\\hspace*{0.2cm}), matching ST2's Class indent.
#' Statistical symbols n, N, and P are italicized in the header.
#'
#' @param data The ternG display tibble; must include the \code{.indent} column
#'   (call \code{ternG()} with \code{indent_info_column = TRUE}).
#'
#' @return A character scalar of LaTeX code ready to \\input.
#'
build_ST3 <- function(data) {
  df <- as.data.frame(data, check.names = FALSE, stringsAsFactors = FALSE)
  # ternG .indent is absolute (e.g. 2 = variable row, 6 = stratum); key off the
  # relative depth so any scheme works. base = variable (GROUP) level.
  ind <- if (".indent" %in% names(df)) df$.indent else rep(0L, nrow(df))
  df$.indent <- NULL
  if ("test" %in% names(df)) df$test <- NULL  # drop ternG's per-row test column (show_test=TRUE); methods live in the caption
  base_ind <- min(ind, na.rm = TRUE)
  # Clean embedded newlines in headers (e.g. "Total\n(N = 68)")
  names(df) <- gsub("\\s*\n\\s*", " ", names(df))
  value_cols <- setdiff(names(df), "Variable")
  # P-value formatting: <0.001 collapses scientific notation; else 3 dp; keep "-"/""
  if ("P" %in% names(df)) {
    num <- suppressWarnings(as.numeric(df$P))
    df$P <- ifelse(!is.na(num) & num < 0.001, "< 0.001",
             ifelse(!is.na(num), formatC(num, format = "f", digits = 3), df$P))
  }
  # LaTeX-escape cell content (% and & appear in n (%) cells; ± in mean ± SD)
  esc <- function(x) {
    x <- ifelse(is.na(x), "", as.character(x))
    x <- gsub("%", "\\\\%", x)
    x <- gsub("&", "\\\\&", x)
    x <- gsub("±", "$\\\\pm$", x)
    x
  }
  vmat <- as.matrix(df[value_cols])
  # Header cells match ST2 exactly: 10pt bold, vertically centered via \raisebox.
  # Group columns become two-line \shortstack (label / count) — that height is what
  # the \raisebox centers against; single-line cells (Variable, P) then center too.
  # Statistical symbols n / N / P are italicized.
  hdr_cell <- function(h) {
    m <- regmatches(h, regexec("^(.*) \\((n|N) = ([0-9]+)\\)$", h))[[1]]
    inner <- if (length(m) == 4L) {
      paste0("\\shortstack{", m[2], "\\\\\\\\(\\textit{", m[3], "} = ", m[4], ")}")
    } else if (h == "P") {
      "\\textit{P}"
    } else {
      h
    }
    paste0("\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{", inner, "}}")
  }
  header_cells <- c(hdr_cell("Variable"), vapply(value_cols, hdr_cell, character(1), USE.NAMES = FALSE))
  header_line  <- paste0("\\rule{0pt}{12pt}", paste(header_cells, collapse = " & "), " \\\\")
  # Body: variable rows plain (underline iff sub-strata follow); strata italic + indented
  body <- vapply(seq_len(nrow(df)), function(i) {
    var  <- df$Variable[i]
    vals <- esc(vmat[i, ])
    if (ind[i] == base_ind) {
      has_sub <- i < nrow(df) && ind[i + 1] > base_ind
      lab <- if (has_sub) paste0("\\underline{", var, "}") else var
      paste0(lab, " & ", paste(vals, collapse = " & "))
    } else {
      paste0("\\hspace*{0.2cm}\\textit{", var, "} & ", paste(vals, collapse = " & "))
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
