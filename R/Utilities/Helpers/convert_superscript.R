#' Convert Superscript Unicode Characters to Regular Letters
#'
#' Transforms Unicode superscript letters (ᵃ-ᶠ) back to their regular lowercase
#' equivalents for text processing and parsing of statistical notation.
#'
#' @param x A character string containing superscript Unicode characters to convert.
#'   Characters that are not superscripts pass through unchanged.
#'
#' @return A character string with superscript letters a-f converted to regular
#'   lowercase equivalents. Other characters remain unchanged.
#'
#' @details
#' This function reverses the transformation from \code{\link{superscript}}, converting:
#' \itemize{
#'   \item ᵃ (\u1d43) → a
#'   \item ᵇ (\u1d47) → b
#'   \item ᶜ (\u1d9c) → c
#'   \item ᵈ (\u1d48) → d
#'   \item ᵉ (\u1d49) → e
#'   \item ᶠ (\u1da0) → f
#' }
#'
#' This is useful for extracting and parsing compact letter display (CLD) labels
#' from post-hoc test results where superscripts indicate statistical groupings.
#'
#' @examples
#' \dontrun{
#' convert_superscript("ᵃ") # Returns "a"
#' convert_superscript("ᵃᵇᶜ") # Returns "abc"
#' convert_superscript("2.45ᵃ") # Returns "2.45a"
#' }
#'
#' @seealso \code{\link{superscript}} for the inverse operation
#'
#' @export
convert_superscript <- function(x) {
  str_replace_all(x, c(
    "ᵃ" = "a",
    "ᵇ" = "b",
    "ᶜ" = "c",
    "ᵈ" = "d",
    "ᵉ" = "e",
    "ᶠ" = "f"
  ))
}