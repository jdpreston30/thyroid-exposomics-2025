#' Screen Chemical Features for Association with a Covariate
#'
#' Tests every chemical feature in a wide feature table for association with a
#' single covariate, using the same per-feature test logic as the main tumor-type
#' analysis in \code{04_variant_stats.R} so that covariate-driven signal is
#' directly comparable to the reported type-driven signal (reviewer question:
#' are the tumor-type differences real, or an artifact of covariate imbalance?).
#'
#' Test selection mirrors the tumor-type analysis exactly: one uniform test per
#' mode, applied to every feature with no per-feature switching.
#' \itemize{
#'   \item \strong{Quantitative:} one-way ANOVA, \code{aov(feature ~ covariate)} --
#'     the identical call used for \code{variant}. Because \code{aov()} wraps
#'     \code{lm()}, a continuous covariate yields the equivalent linear-model
#'     F test on the slope, so no branching by covariate type is needed.
#'   \item \strong{Qualitative, categorical covariate:} Fisher's exact test on
#'     \code{table(covariate, detected)} -- the identical call used for
#'     \code{variant}.
#'   \item \strong{Qualitative, continuous covariate:} Fisher's exact test has no
#'     continuous analogue, so a logistic-regression likelihood-ratio test is used
#'     as the minimal necessary generalization (same framing as Fisher: does
#'     detection depend on the covariate). Matches the qualitative-mode adjusted
#'     models used in the confounding analysis.
#' }
#'
#' Features that cannot be tested return \code{NA_real_} rather than a p-value, so
#' they drop out of any observed-vs-chance denominator instead of counting as
#' non-significant: zero-variance features (e.g. a detection column present in
#' every sample) are skipped up front, and any other test failure is caught by
#' \code{tryCatch} as in the tumor-type analysis.
#'
#' @param features Wide feature table from \code{00c_FTs.R} carrying the sample
#'   ID: \code{tumors_quant_wt_id} or \code{tumors_qual_id}. Every column other
#'   than the ID column and \code{variant} is treated as a chemical feature.
#' @param covariates Data frame of covariates including the ID column, e.g.
#'   \code{clinical_data}. Joined to \code{features} by \code{by}; the join must
#'   match every row of \code{features} or the function errors.
#' @param covariate Name of the covariate column in \code{covariates} to test.
#' @param by Join specification passed to \code{dplyr::left_join}. Defaults to
#'   \code{c("patient_ID" = "Patient_ID")} because the feature tables use
#'   \code{patient_ID} while \code{clinical_data} uses \code{Patient_ID}.
#' @param mode One of "auto" (default), "quant", or "qual". "auto" infers
#'   qualitative mode when every feature value is 0/1 (the same binary heuristic
#'   TernTables uses); inference is done once for the whole table, not per column.
#'
#' @return A tibble with one row per feature: \code{name_sub_lib_id}, \code{mode},
#'   \code{covariate}, \code{test}, \code{p_value}, \code{effect_size},
#'   \code{effect_metric}. Summarize across features at the call site, e.g.
#'   observed significant vs. the ~5% expected by chance.
#'
#'   \code{effect_size} gives the magnitude of the covariate-feature association,
#'   which contextualizes the count of significant features (many weak
#'   associations is a different claim than a few strong ones). The metric is
#'   dictated by the test and named per row in \code{effect_metric}: eta-squared
#'   (ANOVA, fraction of feature variance explained by the covariate), Cramer's V
#'   (Fisher), or McFadden's pseudo-R-squared (logistic LRT). These are
#'   \strong{not} on a common scale, so compare them only within a mode, never
#'   across quantitative and qualitative features.
#'
#' @examples
#' \dontrun{
#' screen_covariate(tumors_quant_wt_id, clinical_data, "year")
#' screen_covariate(tumors_qual_id, clinical_data, "Sex")
#' # observed vs chance across all covariates and both modes
#' tidyr::expand_grid(
#'   tbl = list(tumors_quant_wt_id, tumors_qual_id),
#'   cov = c("Age", "year", "Sex")
#' ) |>
#'   purrr::pmap(\(tbl, cov) screen_covariate(tbl, clinical_data, cov)) |>
#'   dplyr::bind_rows() |>
#'   dplyr::summarise(
#'     n_testable = sum(!is.na(p_value)),
#'     n_sig = sum(p_value < 0.05, na.rm = TRUE),
#'     pct_sig = round(100 * n_sig / n_testable, 1),
#'     expected = round(0.05 * n_testable, 1),
#'     .by = c(mode, covariate, test)
#'   )
#' }
#'
#' @export
screen_covariate <- function(features, covariates, covariate,
                             by = c("patient_ID" = "Patient_ID"), mode = "auto") {
  id_col <- names(by)[1]
  stopifnot(id_col %in% names(features), covariate %in% names(covariates))
  feats <- setdiff(names(features), c(id_col, "variant"))
  stopifnot(length(feats) > 0)
  joined <- features |>
    select(all_of(c(id_col, feats))) |>
    left_join(covariates |> select(all_of(c(unname(by)[1], covariate))), by = by)
  cov_vec <- joined[[covariate]]
  stopifnot(!anyNA(cov_vec)) # every feature-table row must match a covariate row
  fmat <- joined[feats]
  if (identical(mode, "auto")) {
    vals <- unlist(fmat, use.names = FALSE)
    vals <- vals[!is.na(vals)]
    mode <- if (all(vals %in% c(0, 1))) "qual" else "quant"
  }
  mode <- match.arg(mode, c("quant", "qual"))
  cov_is_cat <- is.factor(cov_vec) || is.character(cov_vec) || is.logical(cov_vec)
  test_name <- if (mode == "quant") "ANOVA" else if (cov_is_cat) "Fisher exact" else "Logistic LRT"
  effect_metric <- if (mode == "quant") "eta-squared" else if (cov_is_cat) "Cramer's V" else "McFadden R2"
  res <- vapply(feats, function(f) {
    y <- fmat[[f]]
    if (length(unique(y[!is.na(y)])) < 2) return(c(NA_real_, NA_real_)) # zero-variance feature is untestable
    tryCatch(
      if (mode == "quant") {
        tid <- broom::tidy(aov(y ~ cov_vec))
        c(tid$p.value[1], tid$sumsq[1] / sum(tid$sumsq))
      } else if (cov_is_cat) {
        tab <- table(cov_vec, factor(y))
        chi <- as.numeric(suppressWarnings(chisq.test(tab)$statistic))
        c(fisher.test(tab)$p.value, sqrt(chi / (sum(tab) * (min(dim(tab)) - 1))))
      } else {
        m0 <- suppressWarnings(glm(y ~ 1, family = binomial))
        m1 <- suppressWarnings(glm(y ~ cov_vec, family = binomial))
        c(anova(m0, m1, test = "LRT")[["Pr(>Chi)"]][2], 1 - as.numeric(logLik(m1) / logLik(m0)))
      },
      error = function(e) c(NA_real_, NA_real_)
    )
  }, numeric(2))
  finite_or_na <- function(x) ifelse(is.finite(x), x, NA_real_)
  tibble(
    name_sub_lib_id = feats,
    mode = mode,
    covariate = covariate,
    test = test_name,
    p_value = unname(finite_or_na(res[1, ])),
    effect_size = unname(finite_or_na(res[2, ])),
    effect_metric = effect_metric
  )
}
