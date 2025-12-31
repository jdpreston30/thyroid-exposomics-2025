#' Fix ST3 LaTeX Table Formatting
#'
#' Post-processes gt's LaTeX output for ST3 (landscape):
#' - Convert to longtable for multi-page support  
#' - Fix column specifications for 9 columns
#' - Match ST1/ST2 header formatting
#'
#' @param latex_text Character string of LaTeX code
#'
#' @return Modified LaTeX text
#'
fix_ST3_latex <- function(latex_text) {
  
  # Fix font sizes to 8pt
  latex_text <- gsub("\\\\fontsize\\{6\\.0pt\\}\\{7\\.0pt\\}", "\\\\fontsize{8.0pt}{9.6pt}", latex_text)
  latex_text <- gsub("\\\\fontsize\\{7\\.0pt\\}\\{8\\.4pt\\}", "\\\\fontsize{8.0pt}{9.6pt}", latex_text)
  latex_text <- gsub("\\\\fontsize\\{10\\.0pt\\}\\{12\\.0pt\\}", "\\\\fontsize{8.0pt}{9.6pt}", latex_text)
  latex_text <- gsub("\\\\fontsize\\{12\\.0pt\\}\\{14\\.0pt\\}", "\\\\fontsize{8.0pt}{9.6pt}", latex_text)
  
  # Equalize borders to 0.5pt
  latex_text <- gsub("\\\\toprule(?!\\[)", "\\\\toprule[0.5pt]", latex_text, perl = TRUE)
  latex_text <- gsub("\\\\midrule(?!\\[)", "\\\\midrule[0.5pt]", latex_text, perl = TRUE)
  latex_text <- gsub("\\\\bottomrule", "\\\\bottomrule[0.5pt]", latex_text)
  
  # Convert tabular* to longtable with ST3's 8-column spec
  # Name (left), CAS (center), IARC (center), 5 numeric columns (right)
  latex_text <- gsub("\\\\begin\\{tabular\\*\\}\\{\\\\linewidth\\}\\{@\\{\\\\extracolsep\\{\\\\fill\\}\\}[lrcLRC]+\\}",
                     "\\\\begin{longtable}{lccrrrrr}",
                     latex_text)
  latex_text <- gsub("\\\\begin\\{longtable\\}\\{[clrCLR ]+\\}",
                     "\\\\begin{longtable}{lccrrrrr}",
                     latex_text)
  latex_text <- gsub("\\\\end\\{tabular\\*\\}",
                     "\\\\end{longtable}",
                     latex_text)
  
  # Fix dagger/special characters
  latex_text <- gsub("§", "\\\\textsuperscript{§}", latex_text)
  
  # Replace NA with empty
  latex_text <- gsub("& NA &", "& &", latex_text)
  latex_text <- gsub("& NA \\\\\\\\", "& \\\\\\\\", latex_text)
  
  # Fix headers - add \rule and wrap in raisebox/textbf/shortstack
  # This matches ST1/ST2 formatting exactly
  lines <- strsplit(latex_text, "\\n")[[1]]
  toprule_idx <- grep("\\\\toprule\\[0\\.5pt\\]", lines)
  
  if (length(toprule_idx) > 0) {
    # Find header row (may span multiple lines due to newlines in column names)
    header_start <- toprule_idx + 1
    header_lines <- c()
    current_idx <- header_start
    
    # Collect lines until we hit \\
    while (current_idx <= length(lines) && !grepl("\\\\\\\\\\s*$", lines[current_idx])) {
      header_lines <- c(header_lines, lines[current_idx])
      current_idx <- current_idx + 1
    }
    if (current_idx <= length(lines)) {
      header_lines <- c(header_lines, lines[current_idx])
    }
    header_end <- current_idx
    
    # Combine and clean
    header_content <- paste(header_lines, collapse = " ")
    header_content <- gsub("\\s+", " ", header_content)  # Normalize whitespace
    header_content <- sub("\\s*\\\\\\\\\\s*$", "", header_content)
    
    # Split by &
    headers <- strsplit(header_content, "\\s*&\\s*")[[1]]
    
    # Format each header
    formatted_headers <- sapply(headers, function(h) {
      h <- trimws(h)
      # Add line breaks for specific multi-word headers
      if (h == "IARC Group") {
        h <- "IARC \\\\\\\\ Group"
      } else if (grepl("Mean Tumor", h)) {
        h <- "Mean Tumor \\\\\\\\ Concentration (PPB)"
      } else if (grepl("Range", h)) {
        h <- "Range \\\\\\\\ (PPB)*"
      } else if (grepl("Adipose", h)) {
        h <- "Adipose Tissue \\\\\\\\ (PPB)†"
      } else if (grepl("Urine", h)) {
        h <- "Urine \\\\\\\\ (PPB)‡"
      } else if (grepl("Serum", h)) {
        h <- "Serum/Plasma \\\\\\\\ (PPB)‡"
      }
      paste0("\\raisebox{-0.5\\height}{\\fontsize{10pt}{12pt}\\selectfont\\textbf{\\shortstack{", h, "}}}")
    })
    
    # Create new header line
    new_header <- paste0(
      "\\rule{0pt}{12pt}",
      paste(formatted_headers, collapse = " & "),
      " \\\\"
    )
    
    # Replace in lines
    lines <- c(
      lines[1:toprule_idx],
      new_header,
      lines[(header_end + 1):length(lines)]
    )
    
    latex_text <- paste(lines, collapse = "\n")
  }
  
  return(latex_text)
}
