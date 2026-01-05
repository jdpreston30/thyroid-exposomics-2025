#' Build Supplementary Table 2 (ST2)
#'
#' Formats the metadata table with gt according to journal specifications:
#' - Times New Roman font (applied at PDF compilation)
#' - 8pt content, 10pt bold headers
#' - 0.5pt borders
#' - Hierarchical structure with LaTeX formatting embedded
#' - GROUP headers: bold
#' - Class headers: bold (2-space indent)
#' - Subclass headers: italic + underline (4-space indent)
#' - Chemical names: plain text (6-space indent)
#'
#' @param data The prepared ST2 tibble with LaTeX formatting codes
#'
#' @return A gt table object
#'
build_ST2 <- function(data) {
  
  library(gt)
  library(dplyr)
  
  # Replace NA with empty strings before creating gt table
  data <- data |>
    mutate(across(everything(), ~ ifelse(is.na(.), "", as.character(.))))
  
  # Create gt table
  gt_table <- data |>
    gt() |>
    # Column labels (will be replaced by fix_ST2_latex)
    cols_label(
      Display_Name = "GROUP Class Subclass Name",
      CAS = "CAS",
      `Potential EDC` = "Potential EDC",
      `IARC Group` = "IARC Group",
      `Superclass: Class` = "Superclass: Class"
    ) |>
    # Column alignment: left, center, center, center, left
    cols_align(
      align = "left",
      columns = c(Display_Name, `Superclass: Class`)
    ) |>
    cols_align(
      align = "center",
      columns = c(CAS, `Potential EDC`, `IARC Group`)
    ) |>
    # Table styling for LaTeX output
    tab_options(
      table.font.size = px(8),
      table.border.top.style = "solid",
      table.border.top.width = px(0.5),
      table.border.bottom.style = "solid",
      table.border.bottom.width = px(0.5),
      column_labels.border.top.width = px(0.5),
      column_labels.border.bottom.width = px(0.5),
      column_labels.font.weight = "bold",
      column_labels.font.size = px(10)
    )
  
  return(gt_table)
}
