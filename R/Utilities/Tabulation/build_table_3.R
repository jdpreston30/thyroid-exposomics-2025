#' Build Table 3 (Significant Chemicals)
#'
#' Takes the significant chemicals table and exports it to a formatted Excel file
#' with proper styling, alignment, and formatting.
#'
#' @param data The table 3 tibble with chemical names, usage classes, and statistics
#' @param export_path Path to export the Excel file
#'
#' @return The formatted tibble
#'
build_table_3 <- function(data, export_path) {
  library(openxlsx)
  
  # Clean up Usage Class (Type) values to match standardized format
  data_formatted <- data |>
    mutate(
      `Usage Class (Type)` = case_when(
        # Plasticizer variations
        str_detect(`Usage Class (Type)`, "Plasticizer \\(Plasticizer Metabolites\\)") ~ "Plasticizer Metabolite",
        str_detect(`Usage Class (Type)`, "Plasticizer \\(Plastic Additives\\)") ~ "Plastic Additive",
        str_detect(`Usage Class (Type)`, "Plasticizer \\(Plasticizers\\)") ~ "Plasticizer",
        
        # Combustion Byproduct PAH
        str_detect(`Usage Class (Type)`, "Combustion Byproduct \\(Polycyclic Aromatic Hydrocarbon\\)") ~ "Combustion Byproduct (PAH)",
        
        # PFAS
        str_detect(`Usage Class (Type)`, "Per- and Polyfluoroalkyl Substances \\(PFAS\\)") ~ "PFAS",
        
        # Insecticide/Pesticide variations
        str_detect(`Usage Class (Type)`, "Insecticide/Pesticide \\(Insect Repellents\\)") ~ "Insecticide/Pesticide (Insect Repellent)",
        
        # Organic UV Filters (plural to singular)
        str_detect(`Usage Class (Type)`, "Organic UV Filters") ~ "Organic UV Filter",
        
        # Dye intermediate - capitalize Intermediate
        `Usage Class (Type)` == "Dye intermediate" ~ "Dye Intermediate",
        
        # Preservative - singular
        str_detect(`Usage Class (Type)`, "Preservative \\(Parabens\\)") ~ "Preservative (Paraben)",
        
        # Keep everything else as-is
        TRUE ~ `Usage Class (Type)`
      )
    )
  
  # Note: Periods kept as regular decimal points (not replaced with middle dots)
  
  # Create workbook
  wb <- createWorkbook()
  addWorksheet(wb, "Table 3")
  
  # Write data starting at row 1
  writeData(wb, sheet = 1, x = data_formatted, startRow = 1, colNames = TRUE)
  
  # Create styles
  # Header style - middle vertical alignment, bold, gray fill
  header_style_left <- createStyle(
    fontSize = 7.5,
    fontName = "Times New Roman",
    textDecoration = "bold",
    halign = "left",
    valign = "center",
    fgFill = "#D9D9D9"
  )
  
  header_style_center <- createStyle(
    fontSize = 7.5,
    fontName = "Times New Roman",
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    fgFill = "#D9D9D9"
  )
  
  # Data styles - bottom vertical alignment
  data_style_left <- createStyle(
    fontSize = 7,
    fontName = "Times New Roman",
    halign = "left",
    valign = "bottom"
  )
  
  data_style_center <- createStyle(
    fontSize = 7,
    fontName = "Times New Roman",
    halign = "center",
    valign = "bottom"
  )
  
  # Apply header styles (row 1)
  # Columns 1-2 (Chemical Name, Usage Class) - left aligned
  addStyle(wb, sheet = 1, style = header_style_left, rows = 1, cols = 1:2, gridExpand = TRUE, stack = FALSE)
  # Columns 3-6 (FTC_let, FV-PTC, PTC, p-value) - center aligned
  addStyle(wb, sheet = 1, style = header_style_center, rows = 1, cols = 3:6, gridExpand = TRUE, stack = FALSE)
  
  # Apply data styles (rows 2 onwards)
  data_rows <- 2:(nrow(data_formatted) + 1)
  
  # Columns 1-2 (Chemical Name, Usage Class) - left aligned
  addStyle(wb, sheet = 1, style = data_style_left, rows = data_rows, cols = 1:2, gridExpand = TRUE, stack = FALSE)
  # Columns 3-6 (FTC_let, FV-PTC, PTC, p-value) - center aligned
  addStyle(wb, sheet = 1, style = data_style_center, rows = data_rows, cols = 3:6, gridExpand = TRUE, stack = FALSE)
  
  # Add blank row after data
  blank_row_num <- nrow(data_formatted) + 2
  writeData(wb, sheet = 1, x = "", startRow = blank_row_num, startCol = 1)
  
  # Add double border at bottom of blank row (across all columns)
  blank_row_border <- createStyle(
    border = "bottom",
    borderColour = "black",
    borderStyle = "double"
  )
  addStyle(wb, sheet = 1, style = blank_row_border, rows = blank_row_num, cols = 1:6, gridExpand = TRUE, stack = TRUE)
  
  # Merge cells in row below blank row
  merge_row_num <- blank_row_num + 1
  mergeCells(wb, sheet = 1, cols = 1:6, rows = merge_row_num)
  
  # Add footnote text to merged cell
  footnote_text <- paste(
    "ᵃ Possible, likely, or known carcinogen",
    "ᵇ Potential endocrine disrupting chemical",
    "ᶜ Indicates level 2 identification",
    "Abbreviations: 5-NOT = 5-Nitro-o-toluidine; DEET = N,N-Diethyl-meta-toluamide; DNOP = Di-n-octyl phthalate; MDA = 4,4'-Diaminodiphenylmethane; MEHP = Mono-2-ethylhexyl phthalate; N-MeFOSAA = N-Methylperfluoro-1-octanesulfonamidoacetic acid (linear); OD-PABA = Octyl-dimethyl-p-aminobenzoic acid; PAH = polycyclic aromatic hydrocarbon; TEEP = Tetraethyl ethylenediphosphonate; TTBNP = Tris(tribromoneopentyl); UV = ultraviolet",
    sep = "\n"
  )
  writeData(wb, sheet = 1, x = footnote_text, startRow = merge_row_num, startCol = 1)
  
  # Add style to merged cell with footnote and double border at bottom
  merge_cell_style <- createStyle(
    fontSize = 8,
    fontName = "Times New Roman",
    textDecoration = "italic",
    halign = "left",
    valign = "center",
    border = "bottom",
    borderColour = "black",
    borderStyle = "double",
    wrapText = TRUE
  )
  addStyle(wb, sheet = 1, style = merge_cell_style, rows = merge_row_num, cols = 1:6, gridExpand = TRUE, stack = TRUE)
  
  # Set row height for merged cell to add padding (dynamic based on content)
  setRowHeights(wb, sheet = 1, rows = merge_row_num, heights = 100)
  
  # Add borders to header row (top and bottom)
  header_border_top <- createStyle(
    border = "top",
    borderColour = "black",
    borderStyle = "thin"
  )
  header_border_bottom <- createStyle(
    border = "bottom",
    borderColour = "black",
    borderStyle = "thin"
  )
  addStyle(wb, sheet = 1, style = header_border_top, rows = 1, cols = 1:6, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet = 1, style = header_border_bottom, rows = 1, cols = 1:6, gridExpand = TRUE, stack = TRUE)
  
  # Set column widths
  setColWidths(wb, sheet = 1, cols = 1, widths = 25)  # Chemical Name
  setColWidths(wb, sheet = 1, cols = 2, widths = 30)  # Usage Class
  setColWidths(wb, sheet = 1, cols = 3:5, widths = 12)  # FTC_let, FV-PTC, PTC
  setColWidths(wb, sheet = 1, cols = 6, widths = 12)  # p-value
  
  # Save
  saveWorkbook(wb, export_path, overwrite = TRUE)
  
  # Return the formatted tibble
  data_formatted
}

