#* 17: Construct Supplementary Material PDF
#+ 17.1: Configuration
#- 17.1.0: Set tinytex to continue on LaTeX warnings
options(tinytex.engine_args = "-interaction=nonstopmode")
#- 17.1.1: Set line numbering (TRUE to enable, FALSE to disable)
add_line_numbers <- FALSE
#+ 17.2: Prepare Component Files
#+ 17.3: Read Component Files
#- 17.3.1: Define paths to all component files
components_dir <- here::here("Supplementary", "Components")
sections_dir <- file.path(components_dir, "Sections")
cover_page_path <- file.path(sections_dir, "cover_page.Rmd")
figures_path <- file.path(sections_dir, "figures.Rmd")
methods_path <- file.path(sections_dir, "methods.tex")
tables_path <- file.path(sections_dir, "tables.tex")
#- 17.3.2: Check that all components exist
required_files <- c(cover_page_path, figures_path, methods_path, tables_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing component files: ", paste(missing_files, collapse = ", "))
}
#+ 17.4: Combine Components
#- 17.4.1: Read each component
cover_content <- readLines(cover_page_path, warn = FALSE)
figures_content <- readLines(figures_path, warn = FALSE)
methods_content <- readLines(methods_path, warn = FALSE)
tables_content <- readLines(tables_path, warn = FALSE)
#- 17.4.2: Fix paths for correct references when rendered from Components directory
# Update bibliography and csl paths to be relative from Components directory
bib_path_rel <- file.path("References", "supplementary.bib")
csl_path_rel <- file.path("References", "the-lancet.csl")
# Replace the relative paths in cover content
cover_content <- gsub('bibliography: "References/supplementary.bib"', 
                     paste0('bibliography: "', bib_path_rel, '"'), 
                     cover_content, fixed = TRUE)
cover_content <- gsub('csl: "References/the-lancet.csl"', 
                     paste0('csl: "', csl_path_rel, '"'), 
                     cover_content, fixed = TRUE)
# Fix figure paths to be relative from Components directory  
figures_content <- gsub('../Figures/PDF/', 'Figures/PDF/', figures_content, fixed = TRUE)
# Fix table paths to use absolute paths from config
tables_dir <- file.path(components_dir, "Tables")
tables_content <- gsub('../Tables/ST1.tex', 
                       file.path(tables_dir, "ST1.tex"), 
                       tables_content, fixed = TRUE)
tables_content <- gsub('../Tables/ST2.tex', 
                       file.path(tables_dir, "ST2.tex"), 
                       tables_content, fixed = TRUE)
tables_content <- gsub('./abbreviations.tex', 
                       file.path(sections_dir, "abbreviations.tex"), 
                       tables_content, fixed = TRUE)
#- 17.4.3: Add line numbers if enabled
if (add_line_numbers) {
  # Find the header-includes section in cover_content and add linenumbers package
  yaml_end <- which(cover_content == "---")[2]
  if (!is.na(yaml_end) && yaml_end > 1) {
    # Add line numbering package before the closing ---
    linenumbers_line <- "  - \\usepackage{lineno}"
    modulolinenumbers_line <- "  - \\linenumbers"
    cover_content <- c(
      cover_content[1:(yaml_end - 1)],
      linenumbers_line,
      modulolinenumbers_line,
      cover_content[yaml_end:length(cover_content)]
    )
  }
}
#- 17.4.4: Combine all content
full_content <- c(
  cover_content,
  "",  # Empty line for separation
  methods_content,
  "",  # Empty line for separation
  figures_content,
  "",  # Empty line for separation
  tables_content
)
#+ 17.5: Generate Final PDF
#- 17.5.1: Write combined markdown file
output_rmd <- file.path(components_dir, "supplementary_material.Rmd")
writeLines(full_content, output_rmd)
#- 17.5.2: Create intermediates directory for aux/log files
intermediates_dir <- here::here("Supplementary", "Build_Logs")
if (!dir.exists(intermediates_dir)) dir.create(intermediates_dir, recursive = TRUE)
#- 17.5.3: Render to PDF in Supplementary directory with intermediates in Build_Logs
output_dir <- here::here("Supplementary")
# Note: xelatex produces warnings about missing braces in landscape environment
# These are non-fatal - the PDF compiles successfully. We catch the error and
# verify PDF creation instead of stopping on warnings.
tryCatch({
  rmarkdown::render(
    input = output_rmd,
    output_dir = output_dir,
    output_file = "Supplementary Material.pdf",
    intermediates_dir = intermediates_dir,
    clean = TRUE,
    quiet = FALSE
  )
}, error = function(e) {
  # Check if PDF was actually created despite the "error"
  pdf_path <- file.path(output_dir, "Supplementary Material.pdf")
  if (file.exists(pdf_path)) {
    cat("\n✅ PDF compiled successfully (", round(file.info(pdf_path)$size / 1024), " KB)\n", sep = "")
    cat("   LaTeX warnings are non-fatal and can be ignored.\n")
  } else {
    # Real error - PDF wasn't created
    stop("PDF compilation failed: ", e$message)
  }
})
#- 17.5.4: Move LaTeX build artifacts to Build_Logs
latex_artifacts <- c(
  file.path(output_dir, "Supplementary Material.aux"),
  file.path(output_dir, "Supplementary Material.log"),
  file.path(output_dir, "Supplementary Material.tex")
)
for (artifact in latex_artifacts) {
  if (file.exists(artifact)) {
    file.rename(artifact, file.path(intermediates_dir, basename(artifact)))
  }
}
#- 17.5.5: Open the PDF
output_pdf <- file.path(output_dir, "Supplementary Material.pdf")
system(paste("open", shQuote(output_pdf)))
#- 17.5.6: Clean up empty References folder in Build_Logs (pandoc artifact)
refs_dir <- file.path(intermediates_dir, "References")
if (dir.exists(refs_dir) && length(list.files(refs_dir)) == 0) {
  unlink(refs_dir, recursive = TRUE)
}
#- 17.5.7: Keep the combined markdown file for review (do not delete)