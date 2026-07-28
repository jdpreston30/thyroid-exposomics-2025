#' Test Whether Tumor-Type Effects Survive Covariate Adjustment
#'
#' Refits the tumor-type comparison for a set of chemicals under one or more
#' covariate-adjustment specifications, so the unadjusted effect can be compared
#' directly with the adjusted effect (reviewer question: do the reported
#' differences remain after accounting for age, sex, and collection timing?).
#'
#' The unadjusted model is the same test used in \code{04_variant_stats.R}, so
#' the \code{unadjusted} p-values reproduce the canonical analysis exactly.
#' Adjusted models extend it in the standard way for each mode:
#' \itemize{
#'   \item \strong{Quantitative:} \code{lm(feature ~ type + covariates)} with the
#'     type term assessed by \code{drop1()} F test. With no covariates this is
#'     identical to the canonical one-way ANOVA.
#'   \item \strong{Qualitative:} logistic regression likelihood-ratio test
#'     comparing \code{detected ~ covariates} with \code{detected ~ type +
#'     covariates}. Fisher's exact test cannot accommodate covariates, so the
#'     logistic LRT is the minimal generalization; with no covariates it tests
#'     the same hypothesis as the canonical Fisher test (asymptotically, so its
#'     unadjusted p-values are close to but not identical to Fisher's).
#' }
#'
#' Chemicals that cannot be fit under a given specification return
#' \code{NA_real_} rather than a p-value.
#'
#' @param features Wide feature table carrying the sample ID, i.e.
#'   \code{tumors_quant_wt_id} or \code{tumors_qual_id}.
#' @param covariates Data frame of covariates including the ID column, e.g.
#'   \code{clinical_data}. Must contain \code{type_var} and every covariate named
#'   in \code{adjust}.
#' @param chemicals Character vector of feature columns to test, e.g. the
#'   validated set from \code{MT_final}. Silently intersected with the columns
#'   actually present in \code{features}.
#' @param adjust Named list of adjustment specifications; each element is a
#'   character vector of covariate names (use \code{character(0)} for the
#'   unadjusted model). Names become the \code{adjustment} column.
#' @param type_var Name of the tumor-type column in \code{covariates}.
#' @param by Join specification passed to \code{dplyr::left_join}.
#' @param mode One of "auto" (default), "quant", or "qual"; "auto" infers
#'   qualitative mode when every feature value is 0/1.
#'
#' @return A long tibble with one row per chemical per adjustment:
#'   \code{name_sub_lib_id}, \code{mode}, \code{adjustment}, \code{p_value},
#'   \code{effect_size}, \code{effect_metric}. \code{effect_size} is the type
#'   term's partial eta-squared (quantitative) or McFadden's pseudo-R-squared
#'   (qualitative), letting you show whether an effect that stays significant
#'   also stays the same magnitude. Pivot wider at the call site to compare
#'   adjustments side by side.
#'
#' @examples
#' \dontrun{
#' adjust_type_effect(
#'   tumors_quant_wt_id, clinical_data,
#'   MT_final$name_sub_lib_id[MT_final$mode == "quantitative"]
#' ) |>
#'   tidyr::pivot_wider(id_cols = c(name_sub_lib_id, mode),
#'                      names_from = adjustment, values_from = p_value)
#' }
#'
#' @export
adjust_type_effect <- function(features, covariates, chemicals,
                               adjust = list(
                                 unadjusted = character(0),
                                 year = "year",
                                 full = c("year", "Age", "Sex")
                               ),
                               type_var = "Variant",
                               by = c("patient_ID" = "Patient_ID"),
                               mode = "auto") {
  id_col <- names(by)[1]
  chemicals <- intersect(chemicals, names(features))
  stopifnot(length(chemicals) > 0, id_col %in% names(features))
  cov_needed <- unique(c(type_var, unlist(adjust, use.names = FALSE)))
  stopifnot(all(cov_needed %in% names(covariates)))
  joined <- features |>
    select(all_of(c(id_col, chemicals))) |>
    left_join(covariates |> select(all_of(c(unname(by)[1], cov_needed))), by = by)
  stopifnot(!anyNA(joined[[type_var]])) # every feature-table row must match a covariate row
  fmat <- joined[chemicals]
  if (identical(mode, "auto")) {
    vals <- unlist(fmat, use.names = FALSE)
    vals <- vals[!is.na(vals)]
    mode <- if (all(vals %in% c(0, 1))) "qual" else "quant"
  }
  mode <- match.arg(mode, c("quant", "qual"))
  effect_metric <- if (mode == "quant") "partial eta-squared" else "McFadden R2"
  finite_or_na <- function(x) if (is.finite(x)) x else NA_real_
  fit_one <- function(y, adj) {
    if (length(unique(y[!is.na(y)])) < 2) return(c(NA_real_, NA_real_)) # zero-variance feature is untestable
    d <- joined
    d$.y <- y
    rhs_null <- if (length(adj) == 0) "1" else paste(sprintf("`%s`", adj), collapse = " + ")
    rhs_full <- paste(c(sprintf("`%s`", type_var), if (length(adj)) sprintf("`%s`", adj)), collapse = " + ")
    tryCatch(
      if (mode == "quant") {
        m <- lm(as.formula(paste(".y ~", rhs_full)), data = d)
        dr <- drop1(m, scope = as.formula(paste("~ `", type_var, "`", sep = "")), test = "F")
        ss <- dr[["Sum of Sq"]][2]
        rss <- sum(residuals(m)^2)
        c(finite_or_na(dr[["Pr(>F)"]][2]), finite_or_na(ss / (ss + rss)))
      } else {
        m0 <- suppressWarnings(glm(as.formula(paste(".y ~", rhs_null)), family = binomial, data = d))
        m1 <- suppressWarnings(glm(as.formula(paste(".y ~", rhs_full)), family = binomial, data = d))
        c(
          finite_or_na(anova(m0, m1, test = "LRT")[["Pr(>Chi)"]][2]),
          finite_or_na(1 - as.numeric(logLik(m1) / logLik(m0)))
        )
      },
      error = function(e) c(NA_real_, NA_real_)
    )
  }
  names(adjust) |>
    map(\(a) {
      res <- vapply(chemicals, function(cn) fit_one(fmat[[cn]], adjust[[a]]), numeric(2))
      tibble(
        name_sub_lib_id = chemicals,
        mode = mode,
        adjustment = a,
        p_value = unname(res[1, ]),
        effect_size = unname(res[2, ]),
        effect_metric = effect_metric
      )
    }) |>
    bind_rows()
}
