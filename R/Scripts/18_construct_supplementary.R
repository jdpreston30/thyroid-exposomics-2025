#* 18: Construct Supplementary Material PDF
#+ 18.1: Configuration
#- 18.1.1: Set line numbering (TRUE to enable, FALSE to disable)
add_line_numbers <- FALSE
#+ 18.2: Prepare Component Files
#+ 18.3: Read Component Files
#- 18.3.1: Define paths to all component files
components_dir <- here::here("Supplementary", "Components")
sections_dir <- file.path(components_dir, "Sections")
cover_page_path <- file.path(sections_dir, "cover_page.Rmd")
note1_path <- file.path(sections_dir, "note1.tex")
note2_path <- file.path(sections_dir, "note2.tex")
figures_path <- file.path(sections_dir, "figures.Rmd")
tables_path <- file.path(sections_dir, "tables.tex")
#- 18.3.2: Check that all components exist
required_files <- c(cover_page_path, note1_path, note2_path, figures_path, tables_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing component files: ", paste(missing_files, collapse = ", "))
}
#+ 18.4: Combine Components
#- 18.4.1: Read each component
cover_content <- readLines(cover_page_path, warn = FALSE)
note1_content <- readLines(note1_path, warn = FALSE)
note2_content <- readLines(note2_path, warn = FALSE)
figures_content <- readLines(figures_path, warn = FALSE)
tables_content <- readLines(tables_path, warn = FALSE)
#- 18.4.2: Fix paths for correct references when rendered from Components directory
# Fix figure paths to be relative from Components directory  
figures_content <- gsub('../Figures/PDF/', 'Figures/PDF/', figures_content, fixed = TRUE)
#- 18.4.3: Add line numbers if enabled
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
#- 18.4.4: Combine all content
full_content <- c(
  cover_content,
  "",  # Empty line for separation
  note1_content,     # Supplementary Note 1 (before the figures)
  "",  # Empty line for separation
  note2_content,     # Supplementary Note 2 (confounding/robustness methods)
  "",  # Empty line for separation
  figures_content,
  "",  # Empty line for separation
  tables_content
)
#+ 18.5: Generate Final PDF
#- 18.5.1: Write combined markdown file
output_rmd <- file.path(components_dir, "supplementary_material.Rmd")
writeLines(full_content, output_rmd)
#- 18.5.2: Create intermediates directory for aux/log files
intermediates_dir <- here::here("Supplementary", "Build_Logs")
if (!dir.exists(intermediates_dir)) dir.create(intermediates_dir, recursive = TRUE)
#- 18.5.3: Render to PDF in Supplementary directory with intermediates in Build_Logs
output_dir <- here::here("Supplementary")
rmarkdown::render(
  input = output_rmd,
  output_dir = output_dir,
  output_file = "Supplementary Data.pdf",
  intermediates_dir = intermediates_dir,
  clean = TRUE
)
#- 18.5.4: Move LaTeX build artifacts to Build_Logs
latex_artifacts <- c(
  file.path(output_dir, "Supplementary Data.aux"),
  file.path(output_dir, "Supplementary Data.log"),
  file.path(output_dir, "Supplementary Data.tex")
)
for (artifact in latex_artifacts) {
  if (file.exists(artifact)) {
    file.rename(artifact, file.path(intermediates_dir, basename(artifact)))
  }
}
#- 18.5.5: Open the PDF
output_pdf <- file.path(output_dir, "Supplementary Data.pdf")
system(paste("open", shQuote(output_pdf)))
#- 18.5.6: Clean up empty References folder in Build_Logs (pandoc artifact)
refs_dir <- file.path(intermediates_dir, "References")
if (dir.exists(refs_dir) && length(list.files(refs_dir)) == 0) {
  unlink(refs_dir, recursive = TRUE)
}
#- 18.5.7: Keep the combined markdown file for review (do not delete)
