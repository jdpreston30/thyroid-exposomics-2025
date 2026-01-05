#' Move Superscript Annotation from Mean to End of Mean ± SD Format
#'
#' Reformats statistical summary values by extracting superscript annotations from
#' mean values and appending them after the standard deviation in "mean ± SD" format.
#' This is commonly used for displaying post-hoc test groupings (compact letter displays)
#' alongside summary statistics in tables.
#'
#' @param main Character string containing the mean value, potentially with superscript
#'   annotation (e.g., "25.3ᵃ", "102.7ᵇᶜ").
#' @param sd_raw Numeric value representing the standard deviation to format.
#'
#' @return A character string in the format "mean ± SD" with any superscript annotation
#'   moved to the end (e.g., "25.3 ± 4.12ᵃ", "102.7 ± 15.00ᵇᶜ"). Standard deviation
#'   is formatted to 2 decimal places with trailing zeros preserved.
#'
#' @details
#' This function performs three operations:
#' \enumerate{
#'   \item \strong{Extract superscripts:} Identifies Unicode superscript characters
#'     including letters (ᵃ-ᶠ), numbers (⁰-⁹), and symbols (†, ‡) using regex pattern
#'     \code{[\\u1d43-\\u1d4d\\u02b0-\\u02b8\\u1d62-\\u1d6b\\u2070-\\u209f\\u2020-\\u2021]+}
#'   \item \strong{Format SD:} Rounds standard deviation to 2 decimal places with
#'     trailing zeros (e.g., 4.1 becomes "4.10")
#'   \item \strong{Recombine:} Constructs "mean ± SD" string with superscript appended
#'     at the end
#' }
#'
#' The function handles cases where no superscript is present by simply formatting
#' as "mean ± SD" without additional annotation.
#'
#' @examples
#' \dontrun{
#' # Mean with single superscript letter
#' move_superscript("25.3ᵃ", 4.123) # Returns "25.3 ± 4.12ᵃ"
#'
#' # Mean with multiple superscript letters
#' move_superscript("102.7ᵇᶜ", 15.0) # Returns "102.7 ± 15.00ᵇᶜ"
#'
#' # Mean without superscript
#' move_superscript("50.2", 8.456) # Returns "50.2 ± 8.46"
#'
#' # With superscript symbols (dagger, double-dagger)
#' move_superscript("75.1†", 12.3) # Returns "75.1 ± 12.30†"
#' }
#'
#' @seealso \code{\link{superscript}} for creating superscript annotations,
#'   \code{\link{convert_superscript}} for parsing superscript text
#'
#' @export
move_superscript <- function(main, sd_raw) {
  superscript_pattern <- "[\\u1d43-\\u1d4d\\u02b0-\\u02b8\\u1d62-\\u1d6b\\u2070-\\u209f\\u2020-\\u2021]+"
  # Extract superscript from the main value
  superscript <- str_extract(main, superscript_pattern)
  # Remove superscript from main
  main_clean <- str_remove(main, superscript_pattern)
  # Format SD to 2 decimal places
  sd_clean <- format(round(as.numeric(sd_raw), 2), nsmall = 2)
  # Recombine
  paste0(main_clean, " ± ", sd_clean, ifelse(is.na(superscript), "", superscript))
}