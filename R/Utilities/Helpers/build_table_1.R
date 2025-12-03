#' Build Table 1 (Demographics)
#'
#' Takes the demographic table and exports it to a formatted Excel file
#' with proper styling, order, and formatting.
#'
#' @param data The demographic tibble
#' @param export_path Path to export the Excel file
#'
#' @return The formatted tibble
#'
build_table_1 <- function(data, export_path) {
  library(openxlsx)
  
  # Get the first column name (variable column)
  var_col <- names(data)[1]
  
  # Reorder and rename variables
  data_formatted <- data |>
    mutate(
      !!var_col := case_when(
        .data[[var_col]] == "Age" ~ "Mean Age (Years)",
        .data[[var_col]] == "Sex (Female)" ~ "Sex (Female)",
        .data[[var_col]] == "2006-2009" ~ "    2006-2009",
        .data[[var_col]] == "2010-2013" ~ "    2010-2013",
        .data[[var_col]] == "2014-2017" ~ "    2014-2017",
        .data[[var_col]] == "2018-2021" ~ "    2018-2021",
        TRUE ~ .data[[var_col]]
      )
    )
  
  # Define the desired order
  ordered_vars <- c(
    "Mean Age (Years)",
    "Sex (Female)",
    "Sample Collection Timing",
    "    2006-2009",
    "    2010-2013",
    "    2014-2017",
    "    2018-2021"
  )
  
  # Reorder rows
  data_formatted <- data_formatted |>
    mutate(!!var_col := factor(.data[[var_col]], levels = ordered_vars)) |>
    arrange(.data[[var_col]]) |>
    mutate(!!var_col := as.character(.data[[var_col]]))
  
  # Insert the "Sample Collection Timing" header row
  sample_timing_row <- tibble(
    !!var_col := "Sample Collection Timing",
    Follicular = NA_character_,
    FVPTC = NA_character_,
    Papillary = NA_character_,
    Total = NA_character_
  )
  
  # Combine data with the new header row
  data_formatted <- bind_rows(
    data_formatted |> filter(.data[[var_col]] %in% c("Mean Age (Years)", "Sex (Female)")),
    sample_timing_row,
    data_formatted |> filter(str_detect(.data[[var_col]], "^    "))
  )
  
  # Replace NAs with "-" only for data rows (not Sample Collection Timing)
  for (col in names(data_formatted)[-1]) {  # Skip first column
    data_formatted[[col]] <- ifelse(
      is.na(data_formatted[[col]]) & data_formatted[[var_col]] != "Sample Collection Timing",
      "-",
      data_formatted[[col]]
    )
  }
  
  # Convert decimal points to middle dots (Lancet style)
  data_formatted <- data_formatted |>
    mutate(across(
      everything(),
      ~ str_replace_all(as.character(.x), "\\.", "·")
    ))
  
  # Store the variable names BEFORE renaming columns (for styling logic)
  var_names_for_styling <- data_formatted[[var_col]]
  
  # Update column names with asterisks
  colnames(data_formatted) <- c("Variable: Subcategory", "Follicular*", "FV-PTC*", "Papillary*", "Total*")
  
  # Create workbook
  wb <- createWorkbook()
  addWorksheet(wb, "Table 1")
  
  # Write data starting at row 2 (to leave room for header)
  writeData(wb, sheet = 1, x = data_formatted, startRow = 1, colNames = TRUE)
  
  # Create styles
  header_style <- createStyle(
    fontSize = 10.5,
    fontName = "Times New Roman",
    textDecoration = "bold",
    halign = "center",
    valign = "center"
  )
  
  bold_style <- createStyle(
    fontSize = 10.5,
    fontName = "Times New Roman",
    textDecoration = "bold",
    halign = "left",
    valign = "center"
  )
  
  italic_style <- createStyle(
    fontSize = 10.5,
    fontName = "Times New Roman",
    textDecoration = "italic",
    halign = "left",
    valign = "center"
  )
  
  plain_style_left <- createStyle(
    fontSize = 10.5,
    fontName = "Times New Roman",
    halign = "left",
    valign = "center"
  )
  
  plain_style <- createStyle(
    fontSize = 10.5,
    fontName = "Times New Roman",
    halign = "center",
    valign = "center"
  )
  
  # Apply header style (row 1, all columns)
  addStyle(wb, sheet = 1, style = header_style, rows = 1, cols = 1:5, gridExpand = TRUE, stack = FALSE)
  
  # Apply styles to data rows
  for (i in 1:nrow(data_formatted)) {
    row_num <- i + 1  # +1 for header row
    var_name <- var_names_for_styling[i]  # Use the stored names BEFORE column renaming
    
    # Column 1 (Variable) styling - check for valid var_name
    if (!is.na(var_name) && length(var_name) > 0 && var_name != "") {
      if (var_name %in% c("Mean Age (Years)", "Sex (Female)", "Sample Collection Timing")) {
        addStyle(wb, sheet = 1, style = bold_style, rows = row_num, cols = 1, stack = FALSE)
      } else if (str_detect(var_name, "^    ")) {
        addStyle(wb, sheet = 1, style = italic_style, rows = row_num, cols = 1, stack = FALSE)
      } else {
        # Plain style for any other rows in column 1
        addStyle(wb, sheet = 1, style = plain_style_left, rows = row_num, cols = 1, stack = FALSE)
      }
    } else {
      # If var_name is NA, still apply plain style
      addStyle(wb, sheet = 1, style = plain_style_left, rows = row_num, cols = 1, stack = FALSE)
    }
    
    # Columns 2-5 (data values) - plain style
    addStyle(wb, sheet = 1, style = plain_style, rows = row_num, cols = 2:5, gridExpand = TRUE, stack = FALSE)
  }
  
  # Set column widths
  setColWidths(wb, sheet = 1, cols = 1, widths = 30)
  setColWidths(wb, sheet = 1, cols = 2:5, widths = 15)
  
  # Save
  saveWorkbook(wb, export_path, overwrite = TRUE)
  
  # Return the formatted tibble
  data_formatted
}
