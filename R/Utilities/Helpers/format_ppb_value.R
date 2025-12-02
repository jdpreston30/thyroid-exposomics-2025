#' Format PPB Concentration Values for Display
#'
#' Formats chemical concentration values (in parts per billion) using context-appropriate
#' notation. Small values (< 1 ppb) and large values (≥ 1000 ppb) are displayed in
#' scientific notation, while intermediate values are rounded to integers.
#'
#' @param x Numeric value representing concentration in parts per billion (ppb).
#'   Can be NA.
#'
#' @return A character string containing the formatted concentration value:
#'   \describe{
#'     \item{< 1 ppb}{Scientific notation with 1 decimal place (e.g., "5.0e-01")}
#'     \item{1-999 ppb}{Rounded integer (e.g., "123")}
#'     \item{≥ 1000 ppb}{Scientific notation with 1 decimal place (e.g., "1.5e+03")}
#'     \item{NA}{Returns NA_character_}
#'   }
#'
#' @details
#' This function implements a three-tier formatting strategy to balance precision
#' and readability across a wide range of concentration values:
#' \itemize{
#'   \item \strong{Low concentrations (< 1 ppb):} Scientific notation prevents loss
#'     of significant figures for trace-level detections
#'   \item \strong{Medium concentrations (1-999 ppb):} Integer rounding provides
#'     clean, readable values appropriate for typical detection ranges
#'   \item \strong{High concentrations (≥ 1000 ppb):} Scientific notation avoids
#'     unwieldy large numbers while maintaining precision
#' }
#'
#' @examples
#' \dontrun{
#' # Low concentration - scientific notation
#' format_ppb_value(0.456) # Returns "4.6e-01"
#'
#' # Medium concentration - rounded integer
#' format_ppb_value(123.7) # Returns "124"
#' format_ppb_value(50.2) # Returns "50"
#'
#' # High concentration - scientific notation
#' format_ppb_value(1500) # Returns "1.5e+03"
#' format_ppb_value(50000) # Returns "5.0e+04"
#'
#' # Missing value handling
#' format_ppb_value(NA) # Returns NA_character_
#' }
#'
#' @export
format_ppb_value <- function(x) {
  if (is.na(x)) {
    return(NA_character_)
  }
  if (x < 1) {
    formatC(x, format = "e", digits = 1)
  } else if (x >= 1000) {
    formatC(x, format = "e", digits = 1)
  } else {
    round(x)
  }
}