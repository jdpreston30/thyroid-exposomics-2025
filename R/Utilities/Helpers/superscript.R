#' Convert Letters to Superscript Unicode Characters
#'
#' Transforms lowercase letters (a-f) to their Unicode superscript equivalents
#' for use in statistical notation and compact letter displays from post-hoc tests.
#'
#' @param letter A character string containing letters to convert to superscript.
#'   Letters outside the a-f range are returned unchanged.
#'
#' @return A character string with lowercase letters a-f converted to their
#'   Unicode superscript equivalents. Other characters pass through unchanged.
#'
#' @details
#' This function is primarily used for compact letter displays (CLD) in post-hoc
#' statistical comparisons (e.g., Tukey HSD). The function converts:
#' \itemize{
#'   \item a → ᵃ (\u1d43)
#'   \item b → ᵇ (\u1d47)
#'   \item c → ᶜ (\u1d9c)
#'   \item d → ᵈ (\u1d48)
#'   \item e → ᵉ (\u1d49)
#'   \item f → ᶠ (\u1da0)
#' }
#'
#' @examples
#' \dontrun{
#' superscript("a") # Returns "ᵃ"
#' superscript("abc") # Returns "ᵃᵇᶜ"
#' superscript("a1b2") # Returns "ᵃ1ᵇ2" (numbers unchanged)
#' }
#'
#' @export
superscript <- function(letter) {
  sapply(strsplit(letter, NULL)[[1]], function(char) {
    switch(char,
      "a" = "\u1d43",
      "b" = "\u1d47",
      "c" = "\u1d9c",
      "d" = "\u1d48",
      "e" = "\u1d49",
      "f" = "\u1da0",
      char # fallback if undefined
    )
  }) |> paste0(collapse = "")
}