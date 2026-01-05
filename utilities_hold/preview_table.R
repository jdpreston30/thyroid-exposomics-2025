#' Preview GT Table as PDF
#'
#' Wraps a gt table object in a minimal LaTeX document, compiles to PDF,
#' and automatically opens for quick inspection.
#'
#' @param gt_obj A gt table object
#' @param filename Name for the preview PDF (default: "table_preview.pdf")
#' @param landscape Whether to use landscape orientation (default: TRUE)
#' @param open_pdf Whether to automatically open the PDF (default: TRUE)
#'
#' @return Invisibly returns the path to the generated PDF
#'
preview_table <- function(gt_obj, 
                         filename = "table_preview.pdf", 
                         landscape = TRUE,
                         open_pdf = TRUE) {
  
  # Convert to LaTeX
  latex_code <- gt::as_latex(gt_obj) |> as.character()
  
  # Apply header fill fix
  latex_code <- fix_latex_header_fill(latex_code)
  
  # Build minimal document
  preamble <- c(
    "\\documentclass[10pt]{article}",
    "\\usepackage{booktabs}",
    "\\usepackage{longtable}",
    "\\usepackage{array}",
    "\\usepackage{multirow}",
    "\\usepackage{xcolor}",
    "\\usepackage{colortbl}",
    if (landscape) "\\usepackage{pdflscape}",
    if (landscape) {
      "\\usepackage[landscape,margin=0.3in,paperwidth=11in,paperheight=8.5in]{geometry}"
    } else {
      "\\usepackage[margin=0.3in]{geometry}"
    },
    "\\usepackage{fontspec}",
    "\\setmainfont{Times New Roman}",
    "\\begin{document}",
    latex_code,
    "\\end{document}"
  )
  
  full_doc <- paste(preamble, collapse = "\n")
  
  # Write temp file
  temp_tex <- tempfile(fileext = ".tex")
  writeLines(full_doc, temp_tex)
  
  # Compile with xelatex for Times New Roman support
  output_pdf <- file.path(dirname(temp_tex), filename)
  
  result <- system2(
    "xelatex",
    args = c(
      "-interaction=nonstopmode",
      "-output-directory", dirname(temp_tex),
      paste0("-jobname=", tools::file_path_sans_ext(filename)),
      temp_tex
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  
  # Check if PDF was created
  if (!file.exists(output_pdf)) {
    warning("PDF compilation failed. Check LaTeX errors.")
    message("LaTeX output:\n", paste(result, collapse = "\n"))
    return(invisible(NULL))
  }
  
  message("Preview PDF created: ", output_pdf)
  
  # Open if requested
  if (open_pdf) {
    system2("open", output_pdf)
  }
  
  invisible(output_pdf)
}
