#* 6: Tumor/Cadaver Comparison
#+ 6.1: Estimate PPM/PPB for cadaver thyroid
#- 6.1.1: Normalize by tissue weight, compute PPM/PPB
cadaver_qraw <- cadaver_qraw_i |>
left_join(cadaver_tissue_wts, by = "control_ID") |>
  mutate(
    PPM = (Ce * 10^2) / weight_mg,
    PPB = (Ce * 10^5) / weight_mg
  ) |>
  select(name_sub_lib_id, control_ID, PPM, PPB)
#- 6.1.2: Create PPM feature table
ppm_raw_ctrl <- cadaver_qraw |>
  select(name_sub_lib_id, control_ID, PPM) |>
  pivot_wider(names_from = control_ID, values_from = PPM)
#- 6.1.3: Create PPB feature table
ppb_raw_ctrl <- cadaver_qraw |>
  select(name_sub_lib_id, control_ID, PPB) |>
  pivot_wider(names_from = control_ID, values_from = PPB)
#- 6.1.4: Use to estimate PPM and PPB per chemical in DETECTED samples
ppm_ppb_summary_detected_ctrl <- ppm_raw_ctrl |>
  rowwise() |>
  mutate(
    pct_det_ctrl = sum(!is.na(c_across(-c(name_sub_lib_id)))) / (ncol(ppm_raw_ctrl) - 1) * 100,
    PPM_mean = mean(c_across(-name_sub_lib_id), na.rm = TRUE),
    PPM_sd = sd(c_across(-name_sub_lib_id), na.rm = TRUE),
    PPM_ctrl = paste0(format(PPM_mean, digits = 3), " ± ", format(PPM_sd, digits = 3))
  ) |>
  select(name_sub_lib_id, PPM_ctrl, pct_det_ctrl) |>
  left_join(
    ppb_raw_ctrl |>
      rowwise() |>
      mutate(
        PPB_mean = mean(c_across(-name_sub_lib_id), na.rm = TRUE),
        PPB_sd = sd(c_across(-name_sub_lib_id), na.rm = TRUE),
        PPB_ctrl = paste0(format(PPB_mean, digits = 3), " ± ", format(PPB_sd, digits = 3))
      ) |>
      select(name_sub_lib_id, PPB_ctrl),
    by = "name_sub_lib_id"
  ) |>
  left_join(cas_key_2, by = "name_sub_lib_id") |>
  select(cas, everything(), -c(Name))
#+ 6.2: Estimate PPM/PPB for tumor samples (detected)
#- 6.2.1: Pull the differing features by variant
{
  diff_by_var <- MTi |>
    pull(name_sub_lib_id)
  diff_by_var_quant <- MTi |>
    filter(mode == "quantitative") |>
    pull(name_sub_lib_id)
  diff_by_var_qual <- MTi |>
    filter(mode == "qualitative") |>
    pull(name_sub_lib_id)
}
#- 6.2.2: Pivot long and normalize by tissue weight
ppm_ppb_long <- conc_raw |>
  pivot_longer(cols = -name_sub_lib_id, names_to = "patient_ID", values_to = "Ce") |>
  mutate(C = as.numeric(gsub(",", "", Ce))) |>
  left_join(weights, by = "patient_ID") |>
  mutate(
    PPM = (Ce * 10^2) / weight_mg,
    PPB = (Ce * 10^5) / weight_mg
  ) |>
  select(name_sub_lib_id, patient_ID, PPM, PPB)
#- 6.2.3: Create PPM feature table
ppm_raw <- ppm_ppb_long |>
  select(name_sub_lib_id, patient_ID, PPM) |>
  pivot_wider(names_from = patient_ID, values_from = PPM)
#- 6.2.4: Create PPB feature table
ppb_raw <- ppm_ppb_long |>
  select(name_sub_lib_id, patient_ID, PPB) |>
  pivot_wider(names_from = patient_ID, values_from = PPB)
#- 6.2.5: Use to estimate PPM and PPB per chemical in DETECTED samples
ppm_ppb_summary_detected <- ppm_raw |>
  rowwise() |>
  mutate(
    pct_det = sum(!is.na(c_across(-c(name_sub_lib_id)))) / (ncol(ppm_raw) - 1) * 100,
    PPM_mean = mean(c_across(-name_sub_lib_id), na.rm = TRUE),
    PPM_sd = sd(c_across(-name_sub_lib_id), na.rm = TRUE),
    PPM = paste0(format(PPM_mean, digits = 3), " ± ", format(PPM_sd, digits = 3))
  ) |>
  select(name_sub_lib_id, PPM, pct_det) |>
  left_join(
    ppb_raw |>
      rowwise() |>
      mutate(
        PPB_mean = mean(c_across(-name_sub_lib_id), na.rm = TRUE),
        PPB_sd = sd(c_across(-name_sub_lib_id), na.rm = TRUE),
        PPB = paste0(format(PPB_mean, digits = 3), " ± ", format(PPB_sd, digits = 3))
      ) |>
      select(name_sub_lib_id, PPB),
    by = "name_sub_lib_id"
  ) |>
  left_join(cas_key_2, by = "name_sub_lib_id") |>
  select(cas, everything(), -c(Name))
#+ 6.3: Merge the control and tumor PPM/PPB with the master table
#- 6.3.1: Pare down control to differing by variant
joiner_ctrl <- ppm_ppb_summary_detected_ctrl |>
  filter(name_sub_lib_id %in% diff_by_var) |>
  arrange(name_sub_lib_id)
#- 6.3.2: Pare down tumor to differing by variant
joiner_tumor <- ppm_ppb_summary_detected |>
  filter(name_sub_lib_id %in% diff_by_var) |>
  select(-cas)
#- 6.3.3: Inner join these, use later below
joiner <- joiner_ctrl |>
  inner_join(joiner_tumor, by = "name_sub_lib_id") |>
  select(-name_sub_lib_id)
#+ 6.4: All Inner join (retain features in TUMORS that match controls)
#- 6.4.1: Define the tumor variant columns for taking means
{
ptc_cols <- names(ppm_raw)[grepl("^P\\d+$", names(ppm_raw))]
ftc_cols <- names(ppm_raw)[grepl("^F\\d+$", names(ppm_raw))]
fvptc_cols <- names(ppm_raw)[grepl("^FVPTC\\d+$", names(ppm_raw))]
}
#- 6.4.2: Fully join and gather metadata for detection and means
ppm_full_table <- ppm_raw_ctrl |>
  inner_join(ppm_raw, by = "name_sub_lib_id") |>
  left_join(cas_key_2, by = "name_sub_lib_id") |>
  filter(!cas %in% endog_cas) |>
  arrange(Name) |>
  rowwise() |>
  mutate(
    half_min = (min(c_across(where(is.numeric)), na.rm = TRUE)) / 2,
    max_value = max(c_across(where(is.numeric)), na.rm = TRUE)
  ) |>
  mutate(
    pct_det_ctrl = sum(!is.na(c_across(T001:T009))) / length(c_across(T001:T009)),
    pct_det_tumor = sum(!is.na(c_across(F1:F20))) / length(c_across(F1:F20))
  ) |> # Add info for ideal fragment based on pct_det_tumor, then ctrl, then iMean_tumors below
  left_join(fragment_quality_info, by = "name_sub_lib_id") |>
  arrange(desc(pct_det_tumor), desc(pct_det_ctrl), desc(iMean_tumors)) |>
  group_by(cas) |>
  slice_head(n = 1) |>
  ungroup() |>
  select(-iMean_tumors) |>
  mutate(across(
    T001:F20,
    ~ if_else(is.na(.x), half_min, .x)
  )) |>
  rowwise() |>
  mutate(
    mean_ctrl = mean(c_across(T001:T009), na.rm = TRUE),
    mean_tumor = mean(c_across(F1:F20), na.rm = TRUE),
    mean_PTC = mean(c_across(all_of(ptc_cols)), na.rm = TRUE),
    mean_FTC = mean(c_across(all_of(ftc_cols)), na.rm = TRUE),
    mean_FVPTC = mean(c_across(all_of(fvptc_cols)), na.rm = TRUE)
  ) |>
  ungroup() |>
  left_join(short_name |> select(cas, annot_ident, IARC_Group), by = "cas") |>
  select(Name, cas, annot_ident, IARC_Group, pct_det_ctrl, pct_det_tumor, half_min, max_value, name_sub_lib_id, mean_ctrl:mean_FVPTC, T001:F20) |>
  arrange(cas)
#- 6.4.3: Create a pivoted version for subsequent analysis
full_joiner_i <- ppm_full_table |>
  select(name_sub_lib_id, T001:F20) |>
  pivot_longer(
    cols = -name_sub_lib_id,
    names_to = "sample_ID",
    values_to = "value"
  ) |>
  mutate(
    variant = case_when(
      str_starts(sample_ID, "T") ~ "Ctrl",
      str_starts(sample_ID, "P") ~ "PTC",
      str_starts(sample_ID, "FVPTC") ~ "FV_PTC",
      str_starts(sample_ID, "F") ~ "FTC"
    )
  ) |>
  pivot_wider(names_from = name_sub_lib_id, values_from = value) |>
  mutate(
    tumor_vs_ctrl = ifelse(variant == "Ctrl", "Control", "Tumor")
  ) |>
    
#- 6.4.4: Create version of ppm_full_table except no filtering
ppm_ppb_inclusive <- ppm_raw_ctrl |>
  inner_join(ppm_raw, by = "name_sub_lib_id") |>
  left_join(cas_key_2, by = "name_sub_lib_id") |>
  filter(!cas %in% endog_cas) |>
  arrange(Name) |>
  rowwise() |>
  mutate(
    half_min = (min(c_across(where(is.numeric)), na.rm = TRUE)) / 2,
    max_value = max(c_across(where(is.numeric)), na.rm = TRUE)
  ) |>
  mutate(
    pct_det_ctrl = sum(!is.na(c_across(T001:T009))) / length(c_across(T001:T009)),
    pct_det_tumor = sum(!is.na(c_across(F1:F20))) / length(c_across(F1:F20))
  ) |>
  mutate(across(
    T001:F20,
    ~ if_else(is.na(.x), half_min, .x)
  )) |>
  rowwise() |>
  mutate(
    mean_ctrl_PPB = mean(c_across(T001:T009), na.rm = TRUE),
    mean_tumor_PPB = mean(c_across(F1:F20), na.rm = TRUE),
    mean_PTC_PPB = mean(c_across(all_of(ptc_cols)), na.rm = TRUE),
    mean_FTC_PPB = mean(c_across(all_of(ftc_cols)), na.rm = TRUE),
    mean_FVPTC_PPB = mean(c_across(all_of(fvptc_cols)), na.rm = TRUE)
  ) |>
  ungroup() |>
  arrange(cas) |>
  select(name_sub_lib_id, pct_det_ctrl, pct_det_tumor, mean_ctrl_PPB:mean_FVPTC_PPB) |>
  mutate(across(mean_ctrl_PPB:mean_FVPTC_PPB, ~ .x * 1000)) # convert to ppb
#- 6.4.6: Join temp table with master table
MT_final_ii <- MTi |>
  left_join(ppm_ppb_inclusive, by = "name_sub_lib_id") |>
  select(short_name, cas, name_sub_lib_id,subid, mode, annot_ident, p_value, highest, FTC:PTC, FTC_let:PTC_let, Carcinogenicity:GHS_var_diff_only, Potential_EDC, usage_class:Table_Class, pct_det_ctrl:mean_FVPTC_PPB)
#- 6.4.7: Pull metadata for fragments to IDs
ID_to_frag_cadaver <- IARC_controls_ii |>
  select(id, name_sub_lib_id)
ID_to_frag_tumor <- tumor_raw |>
  select(id, name_sub_lib_id)
#- 6.4.7: Create initial MT_final_i pending validation
MT_final_i <- MT_final_ii |>
  left_join(ID_to_frag_tumor, by = "name_sub_lib_id") |>
  relocate(id, .after = name_sub_lib_id) |>
  mutate(id_subid = paste0(id, "_", subid))
#! Did QC to validate top 5 quant and IARC group 1 detected in both cadaver and tumor
#+ 6.5: Pull IARC Group 1 Features
#- 6.5.1: Get the CAS numbers for IARC 1
iarc_1 <- feature_metadata |>
  filter(IARC_Group == "1") |>
  pull(cas)
#- 6.5.2: Pull the IARC1s in controls
IARC_controls_i <- IARC_controls_ii |>
  filter(cas %in% iarc_1) |>
  rename(pct_NA_ctrl = pct_NA, iMean_ctrl = iMean)
#- 6.5.3: Pull the IARC1s in tumors
IARC_tumors_i <- IARC_tumors_ii |>
  filter(cas %in% iarc_1) |>
  arrange(name)