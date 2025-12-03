#' Build Table 3 with n (%) Counts
#'
#' Creates a hierarchical table with proper indentation and percentage calculations
#' at each level (header -> class -> subclass). Orders by count at each level.
#' Optionally exports to formatted Excel file.
#'
#' @param data Dataframe containing the hierarchical data
#' @param header_col Name of the column containing top-level headers
#' @param class_col Name of the column containing class categories
#' @param subclass_col Name of the column containing subclass categories
#' @param export_path Optional path to export formatted Excel file. If NULL, no export.
#'
#' @return A tibble with two columns: Category (with indentation) and n (%)
#'
build_table_2 <- function(data, header_col, class_col, subclass_col, export_path = NULL) {
  
  # Calculate total n
  total_n <- data |> 
    filter(!is.na(.data[[header_col]])) |> 
    nrow()
  
  # Get counts for headers with custom ordering
  header_counts <- data |>
    filter(!is.na(.data[[header_col]])) |>
    count(.data[[header_col]], name = "n") |>
    rename(Header = 1) |>
    mutate(
      pct = n / total_n * 100,
      # Headers always show 1 decimal place
      pct_display = if_else(pct < 5, "< 5%", paste0(round(pct, 1), "%")),
      display = paste0(n, " (", pct_display, ")")
    )
  
  # Define custom header order based on actual unique values in the data
  custom_header_order <- c(
    "Agrochemicals",
    "Other Chemicals",
    "Pollutants and Industrial Chemicals"
  )
  
  # Apply custom ordering only if the header exists in the data
  header_counts <- header_counts |>
    mutate(Header = factor(Header, levels = custom_header_order)) |>
    arrange(Header) |>
    mutate(Header = as.character(Header))
  
  # Get counts for classes within headers
  class_counts <- data |>
    filter(!is.na(.data[[header_col]]), !is.na(.data[[class_col]])) |>
    count(.data[[header_col]], .data[[class_col]], name = "n") |>
    rename(Header = 1, Class = 2) |>
    left_join(
      data |> 
        filter(!is.na(.data[[header_col]])) |> 
        count(.data[[header_col]], name = "header_n") |>
        rename(Header = 1),
      by = "Header"
    ) |>
    mutate(
      pct = n / header_n * 100,
      pct_display = if_else(pct < 5, "< 5%", paste0(floor(pct + 0.5), "%")),
      display = paste0(n, " (", pct_display, ")")
    ) |>
    arrange(Header, desc(n))
  
  # Get counts for subclasses within classes
  subclass_counts <- data |>
    filter(
      !is.na(.data[[header_col]]), 
      !is.na(.data[[class_col]]), 
      !is.na(.data[[subclass_col]]), 
      .data[[subclass_col]] != ""
    ) |>
    count(.data[[header_col]], .data[[class_col]], .data[[subclass_col]], name = "n") |>
    rename(Header = 1, Class = 2, Subclass = 3) |>
    left_join(
      data |> 
        filter(!is.na(.data[[class_col]])) |> 
        count(.data[[header_col]], .data[[class_col]], name = "class_n") |>
        rename(Header = 1, Class = 2),
      by = c("Header", "Class")
    ) |>
    group_by(Header, Class) |>
    mutate(
      pct = n / class_n * 100,
      # Check if there are duplicate n values (ties) within this class
      has_tie = n_distinct(n) < n(),
      is_tie = n %in% n[duplicated(n)],
      pct_display = if_else(
        pct < 5, 
        "< 5%", 
        # Show .5 decimal only for tied values that round to .5
        if_else(is_tie & round(pct, 1) %% 1 == 0.5, 
                paste0(round(pct, 1), "%"), 
                paste0(floor(pct + 0.5), "%"))  # Standard round half up
      ),
      display = paste0(n, " (", pct_display, ")"),
      # Special flag to force "Other" to the end
      is_other = Subclass == "Other"
    ) |>
    ungroup() |>
    arrange(Header, Class, is_other, desc(n))
  
  # Build the hierarchical table
  table_rows <- list()
  
  for (header in header_counts$Header) {
    # Add header row (ALL CAPS, no spaces)
    table_rows[[length(table_rows) + 1]] <- tibble(
      Category = toupper(header),
      `n (%)` = header_counts$display[header_counts$Header == header]
    )
    
    # Get classes for this header, ordered by count
    classes <- class_counts |> 
      filter(Header == header) |> 
      arrange(desc(n))
    
    if (nrow(classes) > 0) {
      for (i in 1:nrow(classes)) {
      class_name <- classes$Class[i]
      
      # Add class row (4 spaces)
      table_rows[[length(table_rows) + 1]] <- tibble(
        Category = paste0("    ", class_name),
        `n (%)` = classes$display[i]
      )
      
      # Get subclasses for this class, ordered by count (with "Other" always last)
      subclasses <- subclass_counts |> 
        filter(Header == header, Class == class_name) |> 
        arrange(is_other, desc(n))
      
      if (nrow(subclasses) > 0) {
        for (j in 1:nrow(subclasses)) {
          # Add subclass row (8 spaces) with context-specific symbols
          subclass_name <- subclasses$Subclass[j]
          
          # Add symbols based on parent class context
          if (class_name == "Insecticides and Pesticides") {
            if (subclass_name == "Organophosphate") subclass_name <- "Organophosphate*"
            if (subclass_name == "Pyrethroid") subclass_name <- "Pyrethroid†"
            if (subclass_name == "Carbamate") subclass_name <- "Carbamate‡"
          }
          if (class_name == "Herbicides" && subclass_name == "Triazine") {
            subclass_name <- "Triazine¶"
          }
          
          table_rows[[length(table_rows) + 1]] <- tibble(
            Category = paste0("        ", subclass_name),
            `n (%)` = subclasses$display[j]
          )
        }
      }
    }
    }
  }
  
  # Combine all rows and apply formatting
  result <- bind_rows(table_rows) |>
    mutate(
      # Convert decimal points to middle dots (Lancet style)
      `n (%)` = str_replace_all(`n (%)`, "\\.", "·"),
      # Add symbol to Wood Preservatives (class level)
      Category = if_else(
        str_detect(Category, "^\\s*Wood Preservatives$"),
        str_replace(Category, "Wood Preservatives", "Wood Preservatives§"),
        Category
      )
    )
  
  # Export to Excel if path provided
  if (!is.null(export_path)) {
    library(openxlsx)
    
    # Create workbook
    wb <- createWorkbook()
    addWorksheet(wb, "Table")
    
    # Add custom header first (row 1)
    writeData(wb, sheet = 1, x = "GROUP: Class: Type", startRow = 1, startCol = 1, colNames = FALSE)
    writeData(wb, sheet = 1, x = "n (%)", startRow = 1, startCol = 2, colNames = FALSE)
    
    # Write data starting at row 2 without column names
    writeData(wb, sheet = 1, x = result, startRow = 2, colNames = FALSE)
    
    # Create styles
    header_style <- createStyle(fontSize = 8, fontName = "Times New Roman", textDecoration = "bold", halign = "center", valign = "bottom")
    header_style_left <- createStyle(fontSize = 8, fontName = "Times New Roman", textDecoration = "bold", halign = "left", valign = "bottom")
    table_header_style <- createStyle(fontSize = 8, fontName = "Times New Roman", textDecoration = "bold", halign = "left", valign = "bottom")
    table_header_style_center <- createStyle(fontSize = 8, fontName = "Times New Roman", textDecoration = "bold", halign = "center", valign = "bottom")
    table_class_style <- createStyle(fontSize = 8, fontName = "Times New Roman", textDecoration = "bold", halign = "left", valign = "bottom")
    table_class_style_center <- createStyle(fontSize = 8, fontName = "Times New Roman", textDecoration = "bold", halign = "center", valign = "bottom")
    subclass_style <- createStyle(fontSize = 8, fontName = "Times New Roman", textDecoration = "italic", halign = "left", valign = "bottom")
    subclass_style_center <- createStyle(fontSize = 8, fontName = "Times New Roman", textDecoration = "italic", halign = "center", valign = "bottom")
    
    # Apply header style
    addStyle(wb, sheet = 1, style = header_style_left, rows = 1, cols = 1)
    addStyle(wb, sheet = 1, style = header_style, rows = 1, cols = 2)
    
    # Identify row types and apply styles
    for (i in 1:nrow(result)) {
      category <- result$Category[i]
      row_num <- i + 1  # +1 for header row
      
      if (!grepl("^\\s", category)) {
        # No leading spaces = Table Header (ALL CAPS)
        addStyle(wb, sheet = 1, style = table_header_style, rows = row_num, cols = 1)
        addStyle(wb, sheet = 1, style = table_header_style_center, rows = row_num, cols = 2)
      } else if (grepl("^    [^ ]", category)) {
        # 4 spaces = Table Class
        addStyle(wb, sheet = 1, style = table_class_style, rows = row_num, cols = 1)
        addStyle(wb, sheet = 1, style = table_class_style_center, rows = row_num, cols = 2)
      } else if (grepl("^        ", category)) {
        # 8 spaces = Subclass
        addStyle(wb, sheet = 1, style = subclass_style, rows = row_num, cols = 1)
        addStyle(wb, sheet = 1, style = subclass_style_center, rows = row_num, cols = 2)
      }
    }
    
    # Set column widths
    setColWidths(wb, sheet = 1, cols = 1, widths = 60)
    setColWidths(wb, sheet = 1, cols = 2, widths = 15)
    
    # Save
    saveWorkbook(wb, export_path, overwrite = TRUE)
  }
  
  # Return the tibble
  result
}
