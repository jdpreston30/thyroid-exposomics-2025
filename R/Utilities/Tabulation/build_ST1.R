#' Build Supplementary Table 1 (ST1)
#'
#' Formats the chemical library table with gt according to Lancet specifications:
#' - Times New Roman font (applied at PDF compilation)
#' - 7pt content, bold headers, italic mz columns
#' - 0.5pt borders
#' - 4 decimal places for numeric columns (padded with zeros)
#' - Middle dots for decimal points
#'
#' @param data The prepared ST1 tibble
#'
#' @return A gt table object
#'
build_ST1 <- function(data) {
  
  library(gt)
  library(dplyr)
  library(stringr)
  
  # Get mz column names
  mz_cols <- grep("^mz\\d+$", names(data), value = TRUE)
  
  # Format all numeric values to strings with 4 decimal places and middle dots
  data_formatted <- data
  
  # Format RT column
  data_formatted$`Target RT (min.)` <- sapply(data$`Target RT (min.)`, function(x) {
    if (is.na(x)) return("")
    str_replace(sprintf("%.4f", x), "\\.", "·")
  })
  
  # Format mz columns
  for (col in mz_cols) {
    data_formatted[[col]] <- sapply(data[[col]], function(x) {
      if (is.na(x)) return("")
      str_replace(sprintf("%.4f", x), "\\.", "·")
    })
  }
  
  # Format Monoisotopic Mass
  data_formatted$`Monoisotopic Mass` <- sapply(data$`Monoisotopic Mass`, function(x) {
    if (is.na(x) || x == "") return("")
    if (x == "-") return("-")
    num_val <- suppressWarnings(as.numeric(x))
    if (is.na(num_val)) return("")
    str_replace(sprintf("%.4f", num_val), "\\.", "·")
  })
  
  # Create basic gt table
  gt_table <- gt(data_formatted)
  
  return(gt_table)
}
