#' Run Diagnostics
#'
#' Quick helper to source the scratch diagnostic script (\code{.diagnostic/diagnostic.R}).
#' For one-off, exploratory analyses that should not live in the pipeline. Runs in
#' the global environment, so it can read objects already loaded (e.g. clinical_data)
#' and any objects it creates persist in the workspace.
#'
#' @export
#' @examples
#' d()  # or D()
d <- function() {
  path <- ".diagnostic/diagnostic.R"
  if (file.exists(path)) {
    source(path)  # local = FALSE -> evaluates in .GlobalEnv
  } else {
    cat("⚠️  Diagnostic script not found: ", path, "\n", sep = "")
  }
  invisible(NULL)
}

#' @rdname d
#' @export
D <- d
