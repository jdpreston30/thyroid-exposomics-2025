#* 0d: Feature Table Import and Preprocessing
#+ 0d.1: Import raw data and sequence data
tumor_raw <- read_excel(config$paths$primary_data, sheet = "lib.subject.summary")
tumor_seq <- read_excel(config$paths$primary_data, sheet = "tumors_sequence") %>%
  select(ID, variant) %>%
  unique()
#+ 0d.2: Structure data
#- 0d.2.1: Pull the tumor columns
tumor_column <- tumor_raw %>%
  select(name_sub_lib_id, "F1":"F20")
#- 0d.2.2: Pivot longer/wider to get into analysis format
tumor <- tumor_column %>%
  pivot_longer(-name_sub_lib_id, names_to = "ID", values_to = "Value") %>%
  pivot_wider(names_from = name_sub_lib_id, values_from = Value) %>%
  left_join(tumor_seq, by = "ID") %>%
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
tumors_quant <- tumor %>%
  select(all_of(cols_quant)) %>%
  mutate(across(-variant, ~ ifelse(is.na(.), 0.5 * min(., na.rm = TRUE), .))) %>%
  mutate(across(-variant, ~ log2(.)))
#- 0d.3.3: Create qualitative table with binary presence/absence
tumors_qual <- tumor %>%
  select(all_of(c("variant", cols_qual))) %>%
  mutate(across(-variant, ~ ifelse(is.na(.), 0, 1))) %>%
  select(where(~ any(. != 0) | is.character(.)))