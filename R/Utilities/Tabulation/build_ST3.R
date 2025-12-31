#' Build Supplementary Table 3 (ST3)
#'
#' Formats the chemical concentrations and literature comparison table with gt 
#' according to Lancet specifications:
#' - Times New Roman font (applied at PDF compilation)
#' - 8pt content, 10pt bold headers
#' - 0.5pt borders
#' - Middle dots for decimal points
#' - LaTeX citations embedded in data cells
#'
#' @param data The prepared ST3 tibble with literature citations
#'
#' @return A gt table object
#'
build_ST3 <- function(data) {
  
  library(gt)
  library(dplyr)
  library(stringr)
  
  # Create basic gt table
  gt_table <- data |>
    gt() |>
    # Column alignment
    cols_align(
      align = "center",
      columns = c(CAS, `IARC Group`)
    ) |>
    cols_align(
      align = "right",
      columns = c(`Mean Tumor Concentration (PPB)`,
                  `Range (PPB)*`, `Adipose Tissue (PPB)†`, 
                  `Urine (PPB)‡`, `Serum/Plasma (PPB)‡`)
    ) |>
    cols_align(
      align = "left",
      columns = Name
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
