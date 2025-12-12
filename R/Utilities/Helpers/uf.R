#' Reload All Utility Functions
#'
#' Quick helper to re-source all utility functions from R/Utilities/
#' Useful during development when making changes to utility functions
#'
#' @export
#' @examples
#' uf()  # or UF()
uf <- function() {
  utils_path <- "R/Utilities/"
  if (dir.exists(utils_path)) {
    purrr::walk(
      list.files(utils_path, pattern = "\\.[rR]$", full.names = TRUE, recursive = TRUE),
      source
    )
    cat("✅ Utility functions updated\n")
  } else {
    cat("⚠️  Utilities directory not found: ", utils_path, "\n")
  }
  invisible(NULL)
}

#' @rdname uf
#' @export
UF <- uf
