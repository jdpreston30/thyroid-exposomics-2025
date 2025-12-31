#' Fix ST2 LaTeX Table Formatting
#'
#' Post-processes gt's LaTeX output for ST2:
#' - Convert to longtable for multi-page support
#' - Set proper column alignment (left, center, center, center, left)
#' - Add line breaks to headers
#' - Vertical center headers
#' - Set header row height to match ST1 (increased for 4-line first column)
#' - Replace NA with empty strings
#' - Fix checkmark encoding
#'
#' @param latex_text Character string of LaTeX code
#'
#' @return Modified LaTeX text
#'
fix_ST2_latex <- function(latex_text) {
  
  # Fix font size to 8pt for table content
  latex_text <- gsub("\\\\fontsize\\{6\\.0pt\\}\\{7\\.0pt\\}", "\\\\fontsize{8.0pt}{9.6pt}", latex_text)
  latex_text <- gsub("\\\\fontsize\\{7\\.0pt\\}\\{8\\.4pt\\}", "\\\\fontsize{8.0pt}{9.6pt}", latex_text)
  latex_text <- gsub("\\\\fontsize\\{10\\.0pt\\}\\{12\\.0pt\\}", "\\\\fontsize{8.0pt}{9.6pt}", latex_text)
  
  # Equalize borders to 0.5pt
  latex_text <- gsub("\\\\toprule(?!\\[)", "\\\\toprule[0.5pt]", latex_text, perl = TRUE)
  latex_text <- gsub("\\\\midrule(?!\\[)", "\\\\midrule[0.5pt]", latex_text, perl = TRUE)
  
  # Convert tabular* to longtable with fixed column widths
  # Column 1: 8cm left, Column 2: 2cm center (CAS), Columns 3-4: 1.75cm each center (EDC/IARC), Column 5: 9cm left
  latex_text <- gsub("\\\\begin\\{tabular\\*\\}\\{\\\\linewidth\\}\\{@\\{\\\\extracolsep\\{\\\\fill\\}\\}[rlc]+\\}",
                     "\\\\begin{longtable}{>{\\\\raggedright\\\\arraybackslash}p{8cm}>\\\\centering p{2cm}>\\\\centering p{1.75cm}>\\\\centering p{1.75cm}>{\\\\raggedright\\\\arraybackslash}p{9cm}}",
                     latex_text)
  latex_text <- gsub("\\\\end\\{tabular\\*\\}",
                     "\\\\end{longtable}",
                     latex_text)
  
  # Replace NA with empty strings
  latex_text <- gsub("& NA &", "& &", latex_text)
  latex_text <- gsub("& NA \\\\\\\\", "& \\\\\\\\", latex_text)
  
  # Fix dagger escaping
  latex_text <- gsub("\\$\\^\\\\dagger\\$", "$^\\\\dagger$", latex_text)
  
  # Fix checkmark encoding
  latex_text <- gsub("✓", "$\\\\checkmark$", latex_text)
  
  # Replace NA: with -: in Superclass:Class column
  latex_text <- gsub("NA:", "-:", latex_text)
  
  # Split into lines for processing
  lines <- strsplit(latex_text, "\n")[[1]]
  
  # Process data rows to apply formatting and preserve indentation using phantom spaces
  for (i in seq_along(lines)) {
    line <- lines[i]
    # Skip header and structural lines
    if (grepl("^\\\\(toprule|midrule|bottomrule|endfirsthead|endhead|endfoot|addlinespace|rule)", line)) next
    
    # Process lines that contain data
    if (grepl(" & ", line)) {
      parts <- strsplit(line, " & ", fixed = TRUE)[[1]]
      if (length(parts) >= 1) {
        first_col <- parts[1]
        rest <- if(length(parts) > 1) paste(parts[-1], collapse = " & ") else ""
        
        # Detect and apply formatting + indentation using \hspace*{} for spacing
        if (first_col == "PAGEBREAK") {
          # Page break marker - replace entire line with newpage command
          lines[i] <- "\\newpage"
        } else if (grepl("^[A-Z][A-Z ]+$", first_col)) {
          # GROUP (no spaces, all caps): bold, 10pt font, no indent
          first_col <- paste0("{\\fontsize{10pt}{12pt}\\selectfont\\textbf{", first_col, "}}")
          if (length(parts) > 1) {
            lines[i] <- paste(first_col, rest, sep = " & ")
          } else {
            lines[i] <- first_col
          }
        } else if (grepl("^  [A-Z]", first_col)) {
          # Class (2 spaces): bold, with indent
          text <- sub("^  ", "", first_col)
          first_col <- paste0("\\hspace*{0.2cm}\\textbf{", text, "}")
        } else if (grepl("^    [A-Za-z]", first_col)) {
          # Subclass (4 spaces): italic + underline, with indent
          text <- sub("^    ", "", first_col)
          first_col <- paste0("\\hspace*{0.4cm}\\textit{\\underline{", text, "}}")
        } else if (grepl("^      ", first_col)) {
          # Chemical (6 spaces, always): plain text, with indent
          text <- sub("^      ", "", first_col)
          first_col <- paste0("\\hspace*{0.6cm}", text)
        }
        
        # Update line if not already done (for non-GROUP/PageBreak rows)
        if (!grepl("^[A-Z][A-Z ]+$", parts[1]) && parts[1] != "PAGEBREAK") {
          if (length(parts) > 1) {
            lines[i] <- paste(first_col, rest, sep = " & ")
          } else {
            lines[i] <- first_col
          }
        }
      }
    }
  }
  
  # Find and replace header line (look for line with "GROUP Class Subclass Name")
  header_idx <- which(grepl("GROUP Class Subclass Name", lines, fixed = TRUE))
  
  if (length(header_idx) > 0) {
    header_idx <- header_idx[1]
    
    # Extract original headers
    header_line <- lines[header_idx]
    if (grepl(" & ", header_line)) {
      # Build first column header with colon, selective formatting, matching ST1 font size
      first_col_header <- "\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{GROUP: Class:} \\textit{\\underline{Subclass}}: Name}"
      
      # Other headers - all bold, vertically centered, with line breaks where needed
      new_headers <- c(
        first_col_header,
        "\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{CAS}}",
        "\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{\\shortstack{Potential\\\\\\\\EDC}}}",
        "\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{\\shortstack{IARC\\\\\\\\Group}}}",
        "\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{Superclass: Class}}"
      )
      
      # Rebuild header line with ST1 matching height
      lines[header_idx] <- paste(new_headers, collapse = " & ") %>%
        paste0("\\rule{0pt}{12pt}", ., " \\\\")
      
      # Find midrule after header
      midrule_idx <- which(grepl("\\\\midrule\\[0\\.5pt\\]", lines))[1]
      
      if (!is.na(midrule_idx) && midrule_idx > header_idx) {
        continuation_lines <- c(
          "\\addlinespace[2.5pt]",
          "\\endfirsthead",
          "\\toprule[0.5pt]",
          lines[header_idx],
          "\\midrule[0.5pt]",
          "\\addlinespace[2.5pt]",
          "\\endhead",
          "\\bottomrule",
          "\\endfoot"
        )
        
        lines <- c(
          lines[1:midrule_idx],
          continuation_lines,
          lines[(midrule_idx + 1):length(lines)]
        )
      }
    }
  }
  
  # Rejoin lines
  latex_text <- paste(lines, collapse = "\n")
  
  return(latex_text)
}
