#' Decompose a Tumor-Type Effect into Pairwise Contrasts Under Covariate Adjustment
#'
#' Refits each pairwise tumor-type comparison for a set of chemicals under one or
#' more covariate-adjustment specifications. Where an omnibus type effect
#' attenuates after adjustment, this isolates \emph{which} contrast lost
#' significance -- which distinguishes confounding from lost signal when only some
#' type pairs are imbalanced on the covariate. In this cohort only PTC vs IEFVPTC
#' differs on collection year, so attenuation confined to that pair is consistent
#' with confounding control, whereas attenuation of a temporally balanced pair
#' (e.g. FTC vs PTC) is not.
#'
#' Models mirror \code{\link{adjust_type_effect}}, restricted to two-group
#' subsets: \code{lm()} with the type term assessed by \code{drop1()} F test for
#' quantitative features, and a logistic-regression likelihood-ratio test for
#' qualitative detection features.
#'
#' @param features Wide feature table carrying the sample ID, i.e.
#'   \code{tumors_quant_wt_id} or \code{tumors_qual_id}.
#' @param covariates Data frame of covariates including the ID column, e.g.
#'   \code{clinical_data}. Must contain \code{type_var} and every covariate named
#'   in \code{adjust}.
#' @param chemicals Character vector of feature columns to test. Silently
#'   intersected with the columns present in \code{features}.
#' @param adjust Named list of adjustment specifications; each element a character
#'   vector of covariate names (\code{character(0)} for unadjusted). Names become
#'   the \code{adjustment} column.
#' @param type_var Name of the tumor-type column in \code{covariates}.
#' @param by Join specification passed to \code{dplyr::left_join}.
#' @param mode One of "auto" (default), "quant", or "qual"; "auto" infers
#'   qualitative mode when every feature value is 0/1.
#'
#' @return A long tibble with one row per chemical per contrast per adjustment:
#'   \code{name_sub_lib_id}, \code{mode}, \code{pair}, \code{adjustment},
#'   \code{p_value}. Pivot wider on \code{adjustment} to compare side by side.
#'
#' @section Caution:
#' Each contrast is fit on a two-group subset (n = 40 here), so these are lower
#' powered than the three-group models and p-values will be correspondingly
#' larger. Read them as relative patterns across contrasts, not as significance
#' tests in their own right -- particularly for sparsely detected qualitative
#' features, where events per parameter can fall low enough to make the
#' asymptotic likelihood-ratio test unreliable.
#'
#' @examples
#' \dontrun{
#' ancova_summary |>
#'   filter(survives_year, !survives_binned) |>
#'   pull(name_sub_lib_id) |>
#'   (\(chems) pairwise_type_contrasts(tumors_quant_wt_id, clinical_data, chems))()
#' }
#'
#' @export
pairwise_type_contrasts <- function(features, covariates, chemicals,
                                    adjust = list(
                                      unadjusted = character(0),
                                      year = "year",
                                      timing_binned = "Sample_Collection_Timing"
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
  if (identical(mode, "auto")) {
    vals <- unlist(joined[chemicals], use.names = FALSE)
    vals <- vals[!is.na(vals)]
    mode <- if (all(vals %in% c(0, 1))) "qual" else "quant"
  }
  mode <- match.arg(mode, c("quant", "qual"))
  type_levels <- levels(droplevels(factor(joined[[type_var]])))
  stopifnot(length(type_levels) > 2)
  fit_pair <- function(chem, pr, adj) {
    s <- joined[joined[[type_var]] %in% pr, ]
    s[[type_var]] <- droplevels(factor(s[[type_var]]))
    s$.y <- s[[chem]]
    if (length(unique(s$.y[!is.na(s$.y)])) < 2) return(NA_real_)
    rhs_null <- if (length(adj) == 0) "1" else paste(sprintf("`%s`", adj), collapse = " + ")
    rhs_full <- paste(c(sprintf("`%s`", type_var), if (length(adj)) sprintf("`%s`", adj)), collapse = " + ")
    tryCatch(
      if (mode == "quant") {
        m <- lm(as.formula(paste(".y ~", rhs_full)), data = s)
        drop1(m, scope = as.formula(sprintf("~ `%s`", type_var)), test = "F")[["Pr(>F)"]][2]
      } else {
        anova(
          suppressWarnings(glm(as.formula(paste(".y ~", rhs_null)), family = binomial, data = s)),
          suppressWarnings(glm(as.formula(paste(".y ~", rhs_full)), family = binomial, data = s)),
          test = "LRT"
        )[["Pr(>Chi)"]][2]
      },
      error = function(e) NA_real_
    )
  }
  pairs <- combn(type_levels, 2, simplify = FALSE)
  expand_grid(chem = chemicals, adjustment = names(adjust)) |>
    pmap(\(chem, adjustment) {
      adj_spec <- adjust[[adjustment]] # resolve before tibble() recycles `adjustment` to length(pairs)
      tibble(
        name_sub_lib_id = chem,
        mode = mode,
        pair = vapply(pairs, paste, character(1), collapse = " vs "),
        adjustment = adjustment,
        p_value = vapply(pairs, \(pr) fit_pair(chem, pr, adj_spec), numeric(1))
      )
    }) |>
    bind_rows()
}
