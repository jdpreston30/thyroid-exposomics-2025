#' Build Supplementary Table 3 (ST3)
#'
#' Formats the IARC chemicals literature comparison table with gt:
#' - Times New Roman font (applied at PDF compilation)
#' - 8pt content, 10pt bold headers
#' - 0.5pt borders
#' - Hierarchical structure: bold Usage Class headers, plain chemicals (2-space indent)
#' - Single page table (no longtable)
#'
#' @param data The prepared ST3 tibble with hierarchical structure
#'
#' @return A gt table object
#'
build_ST3 <- function(data) {
  
  library(gt)
  library(dplyr)
  
  # Replace NA with empty strings before creating gt table
  data <- data |>
    mutate(across(everything(), ~ ifelse(is.na(.), "", as.character(.))))
  
  # Create gt table
  gt_table <- data |>
    gt() |>
    # Column alignment: left for Name, center for CAS/IARC, right for numeric columns
    cols_align(
      align = "left",
      columns = Name
    ) |>
    cols_align(
      align = "center",
      columns = c(CAS, `IARC Group`)
    ) |>
    cols_align(
      align = "right",
      columns = c(
        `Mean Non-Cancer Thyroid Concentration (PPB)`,
        `Mean Tumor Concentration (PPB)`,
        `Range (PPB)*`,
        `Adipose Tissue (PPB)†`,
        `Urine (PPB)‡`,
        `Serum/Plasma (PPB)‡`
      )
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
