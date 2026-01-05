#' Fix ST3 LaTeX Table Formatting
#'
#' Post-processes gt's LaTeX output for ST3:
#' - Convert to tabular (single page, not longtable)
#' - Set proper column alignment and widths for landscape
#' - Add line breaks to headers
#' - Vertical center headers
#' - Set header row height to match ST1/ST2
#' - Bold Usage Class header rows
#' - Superscript letters (a, b, c, e) in column headers
#' - Match border style from ST1/ST2
#'
#' @param latex_text Character string of LaTeX code
#'
#' @return Modified LaTeX text
#'
fix_ST3_latex <- function(latex_text) {
  
  # Fix font size to 8pt for table content
  latex_text <- gsub("\\\\fontsize\\{6\\.0pt\\}\\{7\\.0pt\\}", "\\\\fontsize{8.0pt}{9.6pt}", latex_text)
  latex_text <- gsub("\\\\fontsize\\{7\\.0pt\\}\\{8\\.4pt\\}", "\\\\fontsize{8.0pt}{9.6pt}", latex_text)
  latex_text <- gsub("\\\\fontsize\\{10\\.0pt\\}\\{12\\.0pt\\}", "\\\\fontsize{8.0pt}{9.6pt}", latex_text)
  
  # Equalize borders to 0.5pt
  latex_text <- gsub("\\\\toprule(?!\\[)", "\\\\toprule[0.5pt]", latex_text, perl = TRUE)
  latex_text <- gsub("\\\\midrule(?!\\[)", "\\\\midrule[0.5pt]", latex_text, perl = TRUE)
  
  # Convert tabular* to tabular with fixed column widths for landscape
  # All columns center except Name (left). Total: 20.5cm (fits in 25.4cm landscape)
  # Widths: Name(4.5cm), CAS(2cm), IARC(1cm), Control(2.75cm), Tumor(2.75cm), Range(2.25cm), AT(1.75cm), Urine(1.75cm), Plasma(1.75cm)
  latex_text <- gsub(
    "\\\\begin\\{tabular\\*\\}\\{\\\\linewidth\\}\\{@\\{\\\\extracolsep\\{\\\\fill\\}\\}[rlc]+\\}",
    "\\\\centering\\\\begin{tabular}{>{\\\\raggedright\\\\arraybackslash}p{4.5cm}>{\\\\centering\\\\arraybackslash}p{2cm}>{\\\\centering\\\\arraybackslash}p{1cm}>{\\\\centering\\\\arraybackslash}p{2.75cm}>{\\\\centering\\\\arraybackslash}p{2.75cm}>{\\\\centering\\\\arraybackslash}p{2.25cm}>{\\\\centering\\\\arraybackslash}p{1.75cm}>{\\\\centering\\\\arraybackslash}p{1.75cm}>{\\\\centering\\\\arraybackslash}p{1.75cm}}",
    latex_text
  )
  latex_text <- gsub("\\\\end\\{tabular\\*\\}", "\\\\end{tabular}", latex_text)
  
  # Replace NA with empty strings
  latex_text <- gsub("& NA &", "& &", latex_text)
  latex_text <- gsub("& NA \\\\\\\\", "& \\\\\\\\", latex_text)
  

  # Fix en-dash
  latex_text <- gsub("–", "--", latex_text)
  

  # Fix escaped \textsuperscript from gt (gt escapes backslashes)
  latex_text <- gsub("\\\\textbackslash\\{\\}textsuperscript\\\\\\{([^}]+)\\\\\\}", "\\\\textsuperscript{\\1}", latex_text)
  
  # Split into lines for processing
  lines <- strsplit(latex_text, "\n")[[1]]
  
  # Process data rows to apply formatting
  for (i in seq_along(lines)) {
    line <- lines[i]
    # Skip header and structural lines
    if (grepl("^\\\\(toprule|midrule|bottomrule|begin|end|rule)", line)) next
    
    # Process lines that contain data
    if (grepl(" & ", line)) {
      parts <- strsplit(line, " & ", fixed = TRUE)[[1]]
      if (length(parts) >= 1) {
        first_col <- parts[1]
        rest <- if(length(parts) > 1) paste(parts[-1], collapse = " & ") else ""
        
        # Detect row type based on first column content
        # Check if it's a Usage Class header: starts with letter (no leading spaces),
        # and all other columns are empty (header rows have empty data)
        other_cols_empty <- length(parts) > 1 && all(sapply(parts[-1], function(x) trimws(x) == "" || grepl("^\\s*\\\\\\\\\\s*$", x)))
        is_header_row <- grepl("^[A-Za-z]", first_col) && 
                         !grepl("^  ", first_col) && 
                         !grepl("^\\\\hspace", first_col) &&
                         other_cols_empty
        
        if (is_header_row) {
          # Usage Class header: bold
          # Examples: "Combustion Byproducts (PAH)", "Dye Intermediates", "Polychlorinated Biphenyls (PCBs)"
          first_col <- paste0("\\textbf{", first_col, "}")
        } else if (grepl("^  [A-Za-z0-9]", first_col)) {
          # Chemical name (2-space indent): apply indent but keep plain text
          # Match letters OR numbers (for chemicals like 2-Naphthylamine, 4-Aminobiphenyl)
          # Use 0.2cm to match ST2 first-level indent
          text <- sub("^  ", "", first_col)
          first_col <- paste0("\\hspace*{0.2cm}", text)
        }
        # Empty/spacer rows stay as-is
        
        # Rebuild line
        if (length(parts) > 1) {
          lines[i] <- paste(first_col, rest, sep = " & ")
        } else {
          lines[i] <- first_col
        }
      }
    }
  }
  
  # Find and replace header line
  header_idx <- which(grepl("^Name & CAS", lines) | grepl("Name & CAS &", lines))
  
  if (length(header_idx) > 0) {
    header_idx <- header_idx[1]
    
    # Build new headers with proper formatting
    # All bold, vertically centered using \raisebox with \shortstack
    # † and ‡ need to be superscripted in headers
    # Using \raisebox{-0.5\height} to vertically center like ST1/ST2
    new_headers <- c(
      "\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{\\shortstack[c]{Name}}}",
      "\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{\\shortstack[c]{CAS}}}",
      "\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{\\shortstack[c]{IARC\\\\Group}}}",
      "\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{\\shortstack[c]{Mean Non-Cancer\\\\Thyroid Conc.\\\\(PPB)}}}",
      "\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{\\shortstack[c]{Mean Tumor\\\\Conc.\\\\(PPB)}}}",
      "\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{\\shortstack[c]{Range\\\\(PPB)$^{\\text{a}}$}}}",
      "\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{\\shortstack[c]{Adipose\\\\Tissue\\\\(PPB)$^{\\text{b}}$}}}",
      "\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{\\shortstack[c]{Urine\\\\(PPB)$^{\\text{c}}$}}}",
      "\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{\\shortstack[c]{Serum/\\\\Plasma\\\\(PPB)$^{\\text{c}}$}}}"
    )
    
    # Rebuild header line with 18pt height for 3-line headers (6pt per line)
    lines[header_idx] <- paste0("\\rule{0pt}{18pt}", paste(new_headers, collapse = " & "), " \\\\")
  }
  
  # Single bottomrule to match ST1/ST2
  bottomrule_idx <- which(grepl("\\\\bottomrule", lines))
  if (length(bottomrule_idx) > 0) {
    lines[bottomrule_idx[1]] <- "\\bottomrule"
  }
  
  # Rejoin lines
  latex_text <- paste(lines, collapse = "\n")
  
  return(latex_text)
}
