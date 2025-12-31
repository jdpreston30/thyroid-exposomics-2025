#' Fix LaTeX Table Formatting
#'
#' Post-processes gt's LaTeX output to fix font size, borders, and make headers bold.
#'
#' @param latex_text Character string of LaTeX code
#' @param fill_color Unused (kept for compatibility)
#'
#' @return Modified LaTeX text
#'
fix_latex_header_fill <- function(latex_text, fill_color = NULL) {
  
  # Fix font size to 8pt for table content
  latex_text <- gsub("\\\\fontsize\\{6\\.0pt\\}\\{7\\.0pt\\}", "\\\\fontsize{8.0pt}{9.6pt}", latex_text)
  latex_text <- gsub("\\\\fontsize\\{7\\.0pt\\}\\{8\\.4pt\\}", "\\\\fontsize{8.0pt}{9.6pt}", latex_text)
  latex_text <- gsub("\\\\fontsize\\{10\\.0pt\\}\\{12\\.0pt\\}", "\\\\fontsize{8.0pt}{9.6pt}", latex_text)
  latex_text <- gsub("\\\\fontsize\\{12\\.0pt\\}\\{14\\.0pt\\}", "\\\\fontsize{8.0pt}{9.6pt}", latex_text)
  
  # Equalize borders to 0.5pt
  latex_text <- gsub("\\\\toprule(?!\\[)", "\\\\toprule[0.5pt]", latex_text, perl = TRUE)
  latex_text <- gsub("\\\\midrule(?!\\[)", "\\\\midrule[0.5pt]", latex_text, perl = TRUE)
  
  # Fix column alignment: first 3 cols left (l,l,l), rest center (c) - 9 columns total
  # Pattern matches: various permutations of l/r/c
  latex_text <- gsub("\\{@\\{\\\\extracolsep\\{\\\\fill\\}\\}[rlc]+\\}", 
                     "{@{\\\\extracolsep{\\\\fill}}lllcccccc}", 
                     latex_text)
  
  # Fix dagger escaping - convert back from gt's over-escaping
  latex_text <- gsub("\\\\\\$\\\\textasciicircum\\{\\}\\\\textbackslash\\{\\}dagger\\\\\\$", 
                     "$^\\\\dagger$", 
                     latex_text)
  
  # Split into lines for easier processing
  lines <- strsplit(latex_text, "\n")[[1]]
  
  # Find toprule line
  toprule_idx <- which(grepl("\\\\toprule\\[0\\.5pt\\]", lines))
  
  if (length(toprule_idx) > 0) {
    # Header should be next line
    header_idx <- toprule_idx + 1
    if (header_idx <= length(lines)) {
      header_line <- lines[header_idx]
      
      # Remove trailing backslash(es) and whitespace
      header_line <- gsub("\\\\+\\s*$", "", header_line)
      
      # Split by & to get columns
      cols <- strsplit(header_line, " & ", fixed = TRUE)[[1]]
      
      # Process each column - make bold and italicize mz columns
      cols_fixed <- sapply(cols, function(col) {
        col <- trimws(col)
        
        # Check if this is an mz column (mz0, mz1, etc.)
        if (grepl("^mz\\d+$", col)) {
          # Make 10pt, italic AND bold
          col <- paste0("{\\fontsize{10pt}{12pt}\\selectfont\\textbf{\\textit{", col, "}}}")
        } else {
          # Just make 10pt bold
          # Escape # if not already escaped
          if (grepl("#", col, fixed = TRUE) && !grepl("\\\\#", col)) {
            col <- gsub("#", "\\\\#", col, fixed = TRUE)
          }
          col <- paste0("{\\fontsize{10pt}{12pt}\\selectfont\\textbf{", col, "}}")
        }
        
        col
      }, USE.NAMES = FALSE)
      
      # Rebuild header line with proper ending
      header_line <- paste(cols_fixed, collapse = " & ")
      header_line <- paste0(header_line, " \\\\ ")
      
      lines[header_idx] <- header_line
    }
  }
  
  # Rejoin lines
  latex_text <- paste(lines, collapse = "\n")
  
  return(latex_text)
}
