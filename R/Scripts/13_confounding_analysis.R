#* 13: Confounding Analysis
#+ 13.1: Covariate Balance and Cross-associations
#- 13.1.1: Pare clinical data to the candidate confounders
clinical_data_confounding <- clinical_data |>
  select(Variant, Sex, Age, year, Sample_Collection_Timing)
#- 13.1.2: Test each covariate against tumor type with the SAME engine as Table 1
#-         (ternG), so test routing (continuous: normality -> Welch ANOVA / Kruskal-
#-         Wallis; categorical: expected counts -> chi-squared / Fisher) and P formatting
#-         are identical to Table 1 and cannot drift. year (continuous) and the binned
#-         Sample_Collection_Timing are both included so the 13.1.3 reconciliation shows.
.bal_tern <- ternG(
  clinical_data_confounding |> select(Variant, Age, Sex, year, Sample_Collection_Timing),
  group_var = "Variant", show_test = TRUE, show_p = TRUE, methods_doc = FALSE,
  open_doc = FALSE, citation = FALSE, print_normality = FALSE,
  indent_info_column = TRUE, smart_rename = FALSE
)
#- 13.1.2a: Harvest each variable's test + formatted P. ternG records the test on the
#-          variable row for continuous vars but on the first sub-row for categorical
#-          vars, so carry the variable name (.indent == 2) down to the first row that
#-          actually carries a test, then map back to display labels (stopifnot guards
#-          against ternG's name-cleaning breaking the match).
.bal_harvest <- .bal_tern |>
  mutate(.var = if_else(.indent == 2, Variable, NA_character_)) |>
  tidyr::fill(.var, .direction = "down") |>
  filter(!is.na(test) & trimws(test) != "" & trimws(test) != "-") |>
  distinct(.var, .keep_all = TRUE)
.nk <- function(x) gsub("[^a-z0-9]", "", tolower(x))
.bal_rows <- .bal_harvest[match(.nk(c("Age", "Sex", "year", "Sample_Collection_Timing")), .nk(.bal_harvest$.var)), ]
stopifnot(!anyNA(.bal_rows$test))
balance_by_type <- tibble(
  covariate = c("Age", "Sex", "Collection year (continuous)", "Collection timing (binned)"),
  test = .bal_rows$test,
  P = .bal_rows$P
)
#- 13.1.3: Year central tendency by type (reconciles a binned split vs a shifted center)
year_by_type <- clinical_data_confounding |>
  group_by(Variant) |>
  summarise(
    median = median(year), mean = round(mean(year), 1),
    min = min(year), max = max(year), n = n(), .groups = "drop"
  )
#- 13.1.4: Cross-associations — are the covariates entangled with each other? No ternG
#-         analogue exists (it only compares variables against the group), so these are
#-         tested directly (exact = FALSE avoids the ties warning) and formatted with the
#-         same TernTables P formatter (val_p_format) so both halves of Table S4 match.
covariate_cross <- tibble(
  pair = c("Age vs Sex", "Age vs year", "year vs Sex"),
  test = c("Wilcoxon rank-sum", "Spearman", "Wilcoxon rank-sum"),
  P = vapply(
    c(
      wilcox.test(Age ~ Sex, clinical_data_confounding, exact = FALSE)$p.value,
      suppressWarnings(cor.test(clinical_data_confounding$Age, clinical_data_confounding$year, method = "spearman")$p.value),
      wilcox.test(year ~ Sex, clinical_data_confounding, exact = FALSE)$p.value
    ),
    val_p_format, character(1)
  )
)
