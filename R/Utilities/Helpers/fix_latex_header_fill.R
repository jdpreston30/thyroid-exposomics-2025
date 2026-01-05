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
  
  # Convert tabular* to longtable for multi-page support
  latex_text <- gsub("\\\\begin\\{tabular\\*\\}\\{\\\\linewidth\\}\\{@\\{\\\\extracolsep\\{\\\\fill\\}\\}[rlc]+\\}",
                     "\\\\begin{longtable}{clllcccccc}",
                     latex_text)
  latex_text <- gsub("\\\\end\\{tabular\\*\\}",
                     "\\\\end{longtable}",
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
      
      # Process each column - make all headers bold (no italic)
      cols_fixed <- sapply(cols, function(col) {
        col <- trimws(col)
        
        # Escape # if not already escaped
        if (grepl("#", col, fixed = TRUE) && !grepl("\\\\#", col)) {
          col <- gsub("#", "\\\\#", col, fixed = TRUE)
        }
        
        # Add line breaks for long column names using shortstack
        if (col == "Monoisotopic Mass") {
          col <- "\\shortstack{Monoisotopic \\\\\\\\ Mass}"
        } else if (col == "Target RT (min.)") {
          col <- "\\shortstack{Target \\\\\\\\ RT (min.)}"
        } else {
          # Single-line headers also get shortstack
          col <- paste0("\\shortstack{", col, "}")
        }
        
        # Make 10pt bold and vertically center with raisebox
        col <- paste0("\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{", col, "}}")
        
        col
      }, USE.NAMES = FALSE)
      
      # Rebuild header line with proper ending
      # Use raisebox for vertical centering in 12pt height header
      header_line <- paste(cols_fixed, collapse = " & ")
      header_line <- paste0("\\rule{0pt}{12pt}", header_line, " \\\\ ")
      
      lines[header_idx] <- header_line
    }
  }
  
  # Rejoin lines
  latex_text <- paste(lines, collapse = "\n")
  
  # For longtable: add commands to define first page header and continuation page headers
  # Find the header and midrule lines in the rejoined text
  if (length(toprule_idx) > 0 && (toprule_idx + 1) <= length(lines)) {
    # Get the header line that was just modified
    header_line_full <- lines[toprule_idx + 1]
    
    # Split text into lines again to insert longtable commands
    text_lines <- strsplit(latex_text, "\n")[[1]]
    
    # Find the midrule line after the header
    midrule_idx <- which(grepl("\\\\midrule\\[0\\.5pt\\]\\\\addlinespace\\[2\\.5pt\\]", text_lines))
    
    if (length(midrule_idx) > 0) {
      midrule_idx <- midrule_idx[1]  # Use first occurrence
      
      # Build the longtable commands
      longtable_commands <- c(
        "\\endfirsthead",
        "\\toprule[0.5pt]",
        header_line_full,
        "\\midrule[0.5pt]\\addlinespace[2.5pt]",
        "\\endhead",
        "\\bottomrule",
        "\\endfoot"
      )
      
      # Insert commands after the first midrule
      text_lines <- c(
        text_lines[1:midrule_idx],
        longtable_commands,
        text_lines[(midrule_idx + 1):length(text_lines)]
      )
      
      latex_text <- paste(text_lines, collapse = "\n")
    }
  }
  
  return(latex_text)
}
