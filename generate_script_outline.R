#' Generate Script Outline
#' 
#' Crawls through all R script files in the Scripts directory and extracts
#' all section headers (lines starting with #*, #+, or #-) to create an
#' outline of the entire analysis pipeline.

library(here)

# Get all R scripts in Scripts directory
scripts_dir <- here("R", "Scripts")
script_files <- list.files(scripts_dir, pattern = "\\.R$", full.names = TRUE)
script_files <- sort(script_files)  # Sort alphabetically

# Initialize output
output_lines <- c()
output_lines <- c(output_lines, paste(rep("=", 80), collapse = ""))
output_lines <- c(output_lines, "THYROID EXPOSOMICS ANALYSIS PIPELINE - SCRIPT OUTLINE")
output_lines <- c(output_lines, paste("Generated:", Sys.time()))
output_lines <- c(output_lines, paste(rep("=", 80), collapse = ""))
output_lines <- c(output_lines, "")

# Process each script file
for (script_path in script_files) {
  # Get script name
  script_name <- basename(script_path)
  
  # Read file contents
  lines <- readLines(script_path, warn = FALSE)
  
  # Extract section headers (lines starting with #*, #+, or #-)
  header_lines <- grep("^#[\\*\\+\\-]", lines, value = TRUE)
  
  # Only include if there are headers
  if (length(header_lines) > 0) {
    # Add script name
    output_lines <- c(output_lines, "")
    output_lines <- c(output_lines, paste(rep("-", 80), collapse = ""))
    output_lines <- c(output_lines, toupper(script_name))
    output_lines <- c(output_lines, paste(rep("-", 80), collapse = ""))
    
    # Add all headers
    output_lines <- c(output_lines, header_lines)
  }
}

# Add footer
output_lines <- c(output_lines, "")
output_lines <- c(output_lines, paste(rep("=", 80), collapse = ""))
output_lines <- c(output_lines, paste("Total scripts processed:", length(script_files)))
output_lines <- c(output_lines, paste(rep("=", 80), collapse = ""))

# Write to output file in root directory
output_file <- here("SCRIPT_OUTLINE.txt")
writeLines(output_lines, output_file)

# Print confirmation
cat("Script outline generated successfully!\n")
cat("Output saved to:", output_file, "\n")
cat("Total scripts processed:", length(script_files), "\n")

# Open the file
system(paste("open", shQuote(output_file)))
