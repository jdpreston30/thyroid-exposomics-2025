#* 8: Generate Supporting Information PDF
#+ 8.1: Read Component Files
#- 8.1.1: Define paths to all component files
components_dir <- here::here("Supporting Information", "Components")
sections_dir <- file.path(components_dir, "Sections")
cover_page_path <- file.path(sections_dir, "cover_page.Rmd")
figures_path <- file.path(sections_dir, "figures.Rmd")
methods_path <- file.path(sections_dir, "methods.tex")
#- 8.1.2: Check that all components exist
required_files <- c(cover_page_path, figures_path, methods_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing component files: ", paste(missing_files, collapse = ", "))
}
#+ 8.2: Combine Components
#- 8.2.1: Read each component
cover_content <- readLines(cover_page_path, warn = FALSE)
figures_content <- readLines(figures_path, warn = FALSE)
methods_content <- readLines(methods_path, warn = FALSE)
#- 8.2.2: Fix paths in cover page for correct references when rendered from Components directory
# Update bibliography and csl paths to be relative from Components directory
references_dir <- file.path(components_dir, "References")
bib_path_rel <- file.path("References", "Supporting_AJT.bib")
csl_path_rel <- file.path("References", "jama.csl")
# Replace the relative paths in cover content
cover_content <- gsub('../References/Supporting_AJT.bib', bib_path_rel, cover_content, fixed = TRUE)
cover_content <- gsub('../References/jama.csl', csl_path_rel, cover_content, fixed = TRUE)
# Fix figure paths to be relative from Components directory  
figures_content <- gsub('../Figures/PDF/', 'Figures/PDF/', figures_content, fixed = TRUE)
#- 8.2.3: Combine all content
full_content <- c(
  cover_content,
  "",  # Empty line for separation
  figures_content,
  "",  # Empty line for separation
  methods_content
)

#+ 8.3: Generate Final PDF
#- 8.3.1: Write combined markdown file
output_rmd <- file.path(components_dir, "supporting_info.Rmd")
writeLines(full_content, output_rmd)

#- 8.3.2: Render to PDF in Supporting Information directory
output_dir <- here::here("Supporting Information")
rmarkdown::render(
  input = output_rmd,
  output_dir = output_dir,
  output_file = "Supporting Information.pdf"
)