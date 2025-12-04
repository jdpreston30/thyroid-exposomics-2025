#* 8: Generate Supplementary Material PDF
#+ 8.1: Read Component Files
#- 8.1.1: Define paths to all component files
components_dir <- here::here("Supplementary", "Components")
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
#- 8.2.2: Fix paths for correct references when rendered from Components directory
# Update bibliography and csl paths to be relative from Components directory
bib_path_rel <- file.path("References", "suplemental.bib")
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
output_rmd <- file.path(components_dir, "supplementary_material.Rmd")
writeLines(full_content, output_rmd)
#- 8.3.2: Render to PDF in Supplementary directory
output_dir <- here::here("Supplementary")
rmarkdown::render(
  input = output_rmd,
  output_dir = output_dir,
  output_file = "Supplementary Material.pdf"
)