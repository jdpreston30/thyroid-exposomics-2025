#* 13: Confounding Analysis
#+ 13.1: Covariate Balance and Cross-associations
#- 13.1.1: Pare clinical data to the candidate confounders
#! Taken from 01_clinical_data.R
clinical_data_confounding <- clinical_data |>
  select(Variant, Sex, Age, year, Sample_Collection_Timing)
#- 13.1.2: Covariate vs. tumor type
balance_by_type <- ternG(
  clinical_data_confounding,
  group_var = "Variant", show_test = TRUE, show_p = TRUE, methods_doc = FALSE,
  open_doc = FALSE, citation = FALSE, plain_tibble = TRUE
) |>
  filter(variable %in% c("Age", "Sex", "year", "Sample_Collection_Timing")) |>
  arrange(match(variable, c("Age", "Sex", "year", "Sample_Collection_Timing"))) |>
  mutate(variable = case_match(variable,
    "year" ~ "Collection year (continuous)",
    "Sample_Collection_Timing" ~ "Collection timing (binned)",
    .default = variable
  )) |>
  select(covariate = variable, test, P = p_fmt)
#- 13.1.3: Year central tendency by type
year_by_type <- clinical_data_confounding |>
  group_by(Variant) |>
  summarise(
    median = median(year), mean = round(mean(year), 1),
    min = min(year), max = max(year), n = n(), .groups = "drop"
  )
#+ 13.2: Cross-associations between covariates
#- 13.2.1: Determine each test based on normality testing
{
  # Normality testing on sex cross-association rows
  sex_rows <- ternG(
    clinical_data_confounding |> select(Sex, Age, year),
    group_var = "Sex", show_test = TRUE, show_p = TRUE, methods_doc = FALSE,
    open_doc = FALSE, citation = FALSE, plain_tibble = TRUE
  )
  # Normality testing on age/year
  ay_test <- classify_normality(clinical_data_confounding |> select(Age, year), group_var = NULL)$is_normal |> all() |> ifelse("Pearson", "Spearman")
}
#- 13.2.2: Determine cross-associations with appropriate tests
covariate_cross <- tibble(
  pair = c("Age vs Sex", "Age vs year", "year vs Sex"),
  test = c(sex_rows$test[1], ay_test, sex_rows$test[2]),
  P = c(
    sex_rows$p_fmt[1],
    val_p_format(suppressWarnings(cor.test(clinical_data_confounding$Age, clinical_data_confounding$year, method = tolower(ay_test))$p.value)),
    sex_rows$p_fmt[2]
  )
)
#+ 13.3: Associations of covariates with chemicals
#- 13.3.1: Screen every feature against each covariate in both modes
#! year tested continuously (1 df linear trend) and binned (3 df) to separate monotone drift from batch structure
covariate_chemical_screens <- c("Age", "year", "Sex", "Sample_Collection_Timing") |>
  map(\(cov) bind_rows(
    screen_covariate(tumors_quant_wt_id, clinical_data, cov),
    screen_covariate(tumors_qual_id, clinical_data, cov)
    #! Both feature tables prepared in 00c_FTs.R
  )) |>
  bind_rows()
#- 13.3.2: Observed significant vs the ~5% expected by chance
covariate_chemical_summary <- covariate_chemical_screens |>
  group_by(mode, covariate, test) |>
  summarise(
    n_testable = sum(!is.na(p_value)),
    n_sig = sum(p_value < 0.05, na.rm = TRUE),
    pct_sig = round(100 * n_sig / n_testable, 1),
    expected = round(0.05 * n_testable, 1),
    effect_metric = first(effect_metric),
    median_effect = round(median(effect_size, na.rm = TRUE), 3),
    median_effect_sig = round(median(effect_size[p_value < 0.05], na.rm = TRUE), 3),
    .groups = "drop"
  ) |>
  mutate(vs_chance = case_when(
    n_sig <= expected * 1.5 ~ "~chance",
    n_sig >= expected * 2 ~ "above chance",
    .default = "modest excess"
  )) |>
  arrange(mode, covariate)
#- 13.3.3: Benchmark tumor type itself on the same effect-size scale as the covariates
variant_benchmark <- bind_rows(
  screen_covariate(tumors_quant_wt_id, clinical_data, "Variant"),
  screen_covariate(tumors_qual_id, clinical_data, "Variant")
) |>
  group_by(mode, covariate, test) |>
  summarise(
    n_testable = sum(!is.na(p_value)),
    n_sig = sum(p_value < 0.05, na.rm = TRUE),
    pct_sig = round(100 * n_sig / n_testable, 1),
    effect_metric = first(effect_metric),
    median_effect = round(median(effect_size, na.rm = TRUE), 3),
    median_effect_sig = round(median(effect_size[p_value < 0.05], na.rm = TRUE), 3),
    .groups = "drop"
  )
#- 13.3.4: Validate the screen reproduces the canonical counts from 04_variant_stats.R
#! Same data, test, and grouping as the canonical analysis, so counts must match exactly
variant_benchmark_check <- variant_benchmark |>
  select(mode, screen_n_sig = n_sig) |>
  left_join(
    tibble(
      mode = c("quant", "qual"),
      canonical_n_sig = c(length(anova_results_sig), n_distinct(fisher_results_i$name_sub_lib_id))
    ),
    by = "mode"
  ) |>
  mutate(match = screen_n_sig == canonical_n_sig)
#+ 13.4.1: Targeted ANCOVA
#- 13.4.1: Specify ANCOVA adjustment variables
#! Continuous year is the primary adjustment; binned timing is a sensitivity check on the Table 1 parameterization
ancova_adjustments <- list(
  unadjusted = character(0),
  year = "year",
  full = c("year", "Age", "Sex"),
  timing_binned = "Sample_Collection_Timing"
)
#- 13.4.2: Refit each validated chemical under each adjustment
ancova_res <- bind_rows(
  adjust_type_effect(tumors_quant_wt_id, clinical_data,
    MT_final$name_sub_lib_id[MT_final$mode == "quantitative"], adjust = ancova_adjustments),
  adjust_type_effect(tumors_qual_id, clinical_data,
    MT_final$name_sub_lib_id[MT_final$mode == "qualitative"], adjust = ancova_adjustments)
)
#- 13.4.3: Detection counts for qualitative chemicals
#! Logistic fits are constrained by min(detected, not detected); few events per parameter makes adjusted p-values unreliable
qual_detection <- tumors_qual_id |>
  select(any_of(MT_final$name_sub_lib_id[MT_final$mode == "qualitative"])) |>
  summarise(across(everything(), sum)) |>
  pivot_longer(everything(), names_to = "name_sub_lib_id", values_to = "n_detected")
#- 13.4.4: Compare unadjusted vs adjusted; flag which lose significance
ancova_summary <- ancova_res |>
  pivot_wider(id_cols = c(name_sub_lib_id, mode), names_from = adjustment,
              values_from = c(p_value, effect_size)) |>
  left_join(MT_final |> select(name_sub_lib_id, short_name), by = "name_sub_lib_id") |>
  left_join(qual_detection, by = "name_sub_lib_id") |>
  mutate(
    survives_year = p_value_year < 0.05,
    survives_binned = p_value_timing_binned < 0.05,
    flips = p_value_unadjusted < 0.05 & !survives_year,
    low_detection = n_detected < 10
  )


#+ 13.5: Why the binned sensitivity attenuates
#- 13.5.1: Association of tumor type with binned collection timing
timing_cramer <- table(clinical_data$Variant, clinical_data$Sample_Collection_Timing) |>
  (\(tab) sqrt(as.numeric(suppressWarnings(chisq.test(tab)$statistic)) / (sum(tab) * (min(dim(tab)) - 1))))()
#- 13.5.2: Variance inflation of the tumor-type term under each adjustment
#! VIF depends only on the predictors, so any chemical as outcome gives the same result
{
  quant_joined <- tumors_quant_wt_id |> left_join(clinical_data, by = c("patient_ID" = "Patient_ID"))
  qual_joined <- tumors_qual_id |> left_join(clinical_data, by = c("patient_ID" = "Patient_ID"))
  vif_outcome <- sprintf("`%s`", MT_final$name_sub_lib_id[MT_final$mode == "quantitative"][1])
  collinearity_check <- c("year", "Sample_Collection_Timing") |>
    map(\(adj) {
      v <- car::vif(lm(as.formula(sprintf("%s ~ Variant + %s", vif_outcome, adj)), data = quant_joined))
      tibble(adjustment = adj, variant_gvif = v["Variant", 1], variant_gvif_scaled = v["Variant", 3])
    }) |>
    bind_rows()
}
#- 13.5.3: Pairwise contrasts for chemicals that attenuate only under binned adjustment
#! Only PTC vs IEFVPTC is temporally imbalanced, so attenuation confined to that contrast is consistent with confounding rather than lost signal
attenuating <- ancova_summary |> filter(survives_year, !survives_binned)
attenuating_pairwise <- bind_rows(
  if (any(attenuating$mode == "quant")) {
    pairwise_type_contrasts(tumors_quant_wt_id, clinical_data, attenuating$name_sub_lib_id[attenuating$mode == "quant"])
  },
  if (any(attenuating$mode == "qual")) {
    pairwise_type_contrasts(tumors_qual_id, clinical_data, attenuating$name_sub_lib_id[attenuating$mode == "qual"])
  }
) |>
  left_join(MT_final |> select(name_sub_lib_id, short_name), by = "name_sub_lib_id") |>
  pivot_wider(id_cols = c(short_name, mode, pair), names_from = adjustment, values_from = p_value)



#+ 13.6: Common support
#! Year adjustment is interpolation only where types share collection years; no overlap would make it extrapolation
common_support <- combn(levels(clinical_data$Variant), 2, simplify = FALSE) |>
  map(\(pr) {
    a <- clinical_data$year[clinical_data$Variant == pr[1]]
    b <- clinical_data$year[clinical_data$Variant == pr[2]]
    lo <- max(min(a), min(b))
    hi <- min(max(a), max(b))
    tibble(
      type_1 = pr[1], type_2 = pr[2],
      overlap = if (lo <= hi) paste0(lo, "-", hi) else "none",
      n_type_1 = sum(a >= lo & a <= hi), n_type_2 = sum(b >= lo & b <= hi),
      total_1 = length(a), total_2 = length(b)
    )
  }) |>
  bind_rows()
