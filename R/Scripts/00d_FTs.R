#* 0d: Feature Table Import and Preprocessing
#+ 0d.1: Import data, clean, preprocess
#- 0d.1.1: Import raw feature table data
tumor_raw <- read_excel(config$paths$primary_data, sheet = "lib.subject.summary")
#- 0d.1.2: Import tumor sequence/variant data
tumor_seq <- read_excel(config$paths$primary_data, sheet = "tumors_sequence") |>
  select(ID, variant) |>
  unique()
#- 0d.1.3: Import feature metadata
feature_metadata <- read_excel(config$paths$chemical_metadata, sheet = "feature_metadata")
#- 0d.1.4: Import and clean tissue weights
weights <- read_excel(config$paths$primary_data, sheet = "tissue_weights") |>
  filter(samples == "Tumor") |>
  select(ID, weight_mg) |>
  rename(patient_ID = ID)
#- 0d.1.5: Import absolute quant data
conc_raw <- read_excel(config$paths$primary_data, sheet = "lib.subject.qsummary", col_type = "text") |>
  select(name_sub_lib_id, F1:F20) |>
  mutate(across(-name_sub_lib_id, as.numeric)) |>
  mutate(across(-name_sub_lib_id, ~ .x * (0.47 / 0.5)))
# Apply correction factor: original pipeline used 0.5 ng/mL, actual is 0.47 ng/mL
#- 0d.1.6: Import cadaver control tissue weights; clean
cadaver_tissue_wts <- read_excel(config$paths$primary_data, sheet = "tissue_weights") |>
  filter(samples == "Control") |>
  select(ID, weight_mg) |>
  rename(control_ID = ID)
#- 0d.1.7: Import cadaver absolute quant, clean
cadaver_qraw_i <- read_excel(config$paths$primary_data, sheet = "lib.subject.qsummary.cadaver") |>
  select(name_sub_lib_id, T001:T009) |>
  mutate(across(-name_sub_lib_id, as.numeric)) |>
  mutate(across(-name_sub_lib_id, ~ .x * (0.47 / 0.5))) |>
  pivot_longer(cols = -name_sub_lib_id, names_to = "control_ID", values_to = "Ce") |>
  mutate(C = as.numeric(gsub(",", "", Ce)))
# Apply correction factor: original pipeline used 0.5 ng/mL, actual is 0.47 ng/mL
#- 0d.1.8: Import in fragment quality info, clean
fragment_quality_info <- read_excel(config$paths$primary_data, sheet = "lib.subject.qsummary") |>
  arrange(cas) |>
  select(name_sub_lib_id, iMean) |>
  rename(iMean_tumors = iMean)
#- 0d.1.9: Import clean absolute quant for control IARCs
IARC_controls_i <- read_excel(config$paths$primary_data, sheet = "lib.subject.qsummary.cadaver")
#- 0d.1.10: Import clean absolute quant for tumor IARCs
IARC_tumors_i <- read_excel(config$paths$primary_data, sheet = "lib.subject.qsummary")
#- 0d.1.11: Import and clean library
ST1_import <- read_excel(config$paths$primary_data, sheet = "library") |>
  filter(Disposition != "Endogenous") |>
  mutate(subid_col = paste0("mz", subid)) |>
  select(id, name, short_display_name, trt, monoisotopic, cas, formula, Disposition, subid_col, tmz) |>
  distinct() |>
  pivot_wider(
    names_from = subid_col,
    values_from = tmz
  ) |>
  arrange(cas)
#- 0d.1.12: Import File List
file_list <- read_excel(config$paths$primary_data, sheet = "file_list") |>
  select(file, ID, replicate, study, type) |>
  group_by(ID, study, type) |>
  summarize(files = paste(file, collapse = ", "), .groups = "drop")
#- 0d.1.13: Import GC2 feature list
GC2_features <- read_csv(config$paths$gc2_features)
#- 0d.1.14: Import expanded library
expanded_lib <- read_csv(config$paths$gc2_expanded)
#- 0d.1.15: Validation
validation_check <- read_xlsx(config$paths$variant_validation, sheet = "validation")
#+ 0d.2: Structure data
#- 0d.2.1: Pull the tumor columns
tumor_column <- tumor_raw |>
  select(name_sub_lib_id, "F1":"F20")
#- 0d.2.2: Pivot longer/wider to get into analysis format
tumor <- tumor_column |>
  pivot_longer(-name_sub_lib_id, names_to = "ID", values_to = "Value") |>
  pivot_wider(names_from = name_sub_lib_id, values_from = Value) |>
  left_join(tumor_seq, by = "ID") |>
  select(ID, variant, everything(), -ID)
#+ 0d.3: Split into quantitative and qualitative features based on missingness
#- 0d.3.1: Calculate proportion of missing values per feature
{
  na_threshold <- 0.3
  PMV <- colMeans(is.na(tumor))
  cols_quant <- names(PMV[PMV <= na_threshold])
  cols_qual <- names(PMV[PMV > na_threshold])
  cols_quant <- unique(c("variant", cols_quant))
  cols_qual <- setdiff(cols_qual, "variant")
}
#- 0d.3.2: Create quantitative table with 1/2 minimum imputation and log2 transform
tumors_quant <- tumor |>
  select(all_of(cols_quant)) |>
  mutate(across(-variant, ~ ifelse(is.na(.), 0.5 * min(., na.rm = TRUE), .))) |>
  mutate(across(-variant, ~ log2(.)))
#- 0d.3.3: Create qualitative table with binary presence/absence
tumors_qual <- tumor |>
  select(all_of(c("variant", cols_qual))) |>
  mutate(across(-variant, ~ ifelse(is.na(.), 0, 1))) |>
  select(where(~ any(. != 0) | is.character(.)))
#+ 0d.4: Clean up targeted feature tables
#- 0d.4.1: Normalize by weights
tumors_quant_wt_i <- tumor_column |>
  pivot_longer(-name_sub_lib_id, names_to = "ID", values_to = "Value") |>
  pivot_wider(names_from = name_sub_lib_id, values_from = Value) |>
  mutate(patient_ID = ID) |>
  left_join(tumor_seq, by = "ID") |>
  select(patient_ID, all_of(cols_quant)) |>
  mutate(across(-c(variant, patient_ID), ~ as.numeric(.))) |>
  mutate(across(-c(variant, patient_ID), ~ ifelse(is.na(.), 0.5 * min(., na.rm = TRUE), .))) |>
  left_join(weights, by = "patient_ID") |> # join tissue weights in
  mutate(across(where(is.numeric) & !matches("weight_mg"), ~ . / weight_mg)) |>
  select(-c(weight_mg, patient_ID)) |>
  mutate(across(-variant, ~ log2(.)))
#- 0d.4.2: Pull endogenous features
endog_cas <- read_excel(config$paths$chemical_metadata, sheet = "endogenous_excluded_features") |>
  pull(cas)
#- 0d.4.3: Get endogenous feature key
cas_key_endog <- tumor_raw |>
  select(name_sub_lib_id, cas) |>
  filter(cas %in% endog_cas) |>
  pull(name_sub_lib_id)
#- 0d.4.4: Remove endogenous features from quantitative weighted table
tumors_quant_wt <- tumors_quant_wt_i |>
  select(-any_of(cas_key_endog))
