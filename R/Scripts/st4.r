full_joiner_i

ST3_tibble <- ppm_full_table |>
  select()

ppm_ppb_inclusive





#+ Determine est PPM/PPB for detected samples and controls, join (ST3/4)
#- Import the quantified data
#_Pull the differing features by variant
diff_by_var <- MTi %>%
  pull(name_sub_lib_id)
diff_by_var_quant <- MTi %>%
  filter(mode == "quantitative") %>%
  pull(name_sub_lib_id)
diff_by_var_qual <- MTi %>%
  filter(mode == "qualitative") %>%
  pull(name_sub_lib_id)
#_Import raw concentration data
conc_raw <- read_excel("Data and Metadata Files/primary_data.xlsx", sheet = "lib.subject.qsummary",col_type = "text") %>%
  select(name_sub_lib_id,F1:F20) %>%
  mutate(across(-name_sub_lib_id, as.numeric))
#- Create a mean PPM/PPB in detected samples value for each chemical, add to master table
#_Pivot long and normalize by weight
ppm_ppb_long <- conc_raw %>%
  pivot_longer(cols = -name_sub_lib_id, names_to = "patient_ID", values_to = "Ce") %>%
  mutate(C = as.numeric(gsub(",", "", Ce))) %>%
  left_join(weights, by = "patient_ID") %>%
  mutate(
    PPM = (Ce * 10^2) / weight_mg,
    PPB = (Ce * 10^5) / weight_mg
  ) %>%
  select(name_sub_lib_id, patient_ID, PPM, PPB)
#_Convert to PPM and PPB feature tables
ppm_raw <- ppm_ppb_long %>%
  select(name_sub_lib_id, patient_ID, PPM) %>%
  pivot_wider(names_from = patient_ID, values_from = PPM)
ppb_raw <- ppm_ppb_long %>%
  select(name_sub_lib_id, patient_ID, PPB) %>%
  pivot_wider(names_from = patient_ID, values_from = PPB)
#_Now, use these to make estimated PPM and PPB per chemical in DETECTED samples
ppm_ppb_summary_detected <- ppm_raw %>%
  rowwise() %>%
  mutate(
    pct_det = sum(!is.na(c_across(-c(name_sub_lib_id)))) / (ncol(.) - 1) * 100,
    PPM_mean = mean(c_across(-name_sub_lib_id), na.rm = TRUE),
    PPM_sd = sd(c_across(-name_sub_lib_id), na.rm = TRUE),
    PPM = paste0(format(PPM_mean, digits = 3), " ± ", format(PPM_sd, digits = 3))
  ) %>%
  select(name_sub_lib_id, PPM,pct_det) %>%
  left_join(
    ppb_raw %>%
      rowwise() %>%
      mutate(
        PPB_mean = mean(c_across(-name_sub_lib_id), na.rm = TRUE),
        PPB_sd = sd(c_across(-name_sub_lib_id), na.rm = TRUE),
        PPB = paste0(format(PPB_mean, digits = 3), " ± ", format(PPB_sd, digits = 3))
      ) %>%
      select(name_sub_lib_id, PPB),
    by = "name_sub_lib_id"
  ) %>%
  left_join(cas_key_2, by = "name_sub_lib_id") %>%
  select(cas, everything(),-c(Name))
#- Determine PPM and PPB for controls
#_Import control tissue weights
cadaver_tissue_wts <- read_excel("Data and Metadata Files/primary_data.xlsx", sheet = "tissue_weights") %>%
  filter(samples == "Control") %>%
  select(ID, weight_mg) %>%
  rename(control_ID = ID)
#_Import control conc data, normalize by tissue weight, compute PPM/PPB
cadaver_qraw <- read_excel("Data and Metadata Files/primary_data.xlsx", sheet = "lib.subject.qsummary.cadaver") %>%
  select(name_sub_lib_id, T001:T009) %>%
  mutate(across(-name_sub_lib_id, as.numeric)) %>%
  pivot_longer(cols = -name_sub_lib_id, names_to = "control_ID", values_to = "Ce") %>%
  mutate(C = as.numeric(gsub(",", "", Ce))) %>%
  left_join(cadaver_tissue_wts, by = "control_ID") %>%
  mutate(
    PPM = (Ce * 10^2) / weight_mg,
    PPB = (Ce * 10^5) / weight_mg
  ) %>%
  select(name_sub_lib_id, control_ID, PPM, PPB)
#_Convert to PPM and PPB feature tables
ppm_raw_ctrl <- cadaver_qraw %>%
  select(name_sub_lib_id, control_ID, PPM) %>%
  pivot_wider(names_from = control_ID, values_from = PPM)
ppb_raw_ctrl <- cadaver_qraw %>%
  select(name_sub_lib_id, control_ID, PPB) %>%
  pivot_wider(names_from = control_ID, values_from = PPB)
#_Now, use these to make estimated PPM and PPB per chemical in DETECTED samples
ppm_ppb_summary_detected_ctrl <- ppm_raw_ctrl %>%
  rowwise() %>%
  mutate(
    pct_det_ctrl = sum(!is.na(c_across(-c(name_sub_lib_id)))) / (ncol(.) - 1) * 100,
    PPM_mean = mean(c_across(-name_sub_lib_id), na.rm = TRUE),
    PPM_sd = sd(c_across(-name_sub_lib_id), na.rm = TRUE),
    PPM_ctrl = paste0(format(PPM_mean, digits = 3), " ± ", format(PPM_sd, digits = 3))
  ) %>%
  select(name_sub_lib_id, PPM_ctrl, pct_det_ctrl) %>%
  left_join(
    ppb_raw_ctrl %>%
      rowwise() %>%
      mutate(
        PPB_mean = mean(c_across(-name_sub_lib_id), na.rm = TRUE),
        PPB_sd = sd(c_across(-name_sub_lib_id), na.rm = TRUE),
        PPB_ctrl = paste0(format(PPB_mean, digits = 3), " ± ", format(PPB_sd, digits = 3))
      ) %>%
      select(name_sub_lib_id, PPB_ctrl),
    by = "name_sub_lib_id"
  ) %>%
  left_join(cas_key_2, by = "name_sub_lib_id") %>%
  select(cas, everything(), -c(Name))
#- Merge the control and tumor PPM/PPB with the master table
#_Create pared down filtered control and tumor for only variant diff
joiner_ctrl <- ppm_ppb_summary_detected_ctrl %>%
  filter(name_sub_lib_id %in% diff_by_var) %>%
  arrange(name_sub_lib_id)
joiner_tumor <- ppm_ppb_summary_detected %>%
  filter(name_sub_lib_id %in% diff_by_var) %>%
  select(-cas)
#_Inner join these, use later below
joiner <- joiner_ctrl %>%
  inner_join(joiner_tumor, by = "name_sub_lib_id") %>%
  select(-name_sub_lib_id)
#+ All Inner join (only retain features in TUMORS that match controls) (ST3/4)
#- Define the tumor variant columns for taking means
ptc_cols <- names(ppm_raw)[grepl("^P\\d+$", names(ppm_raw))]
ftc_cols <- names(ppm_raw)[grepl("^F\\d+$", names(ppm_raw))]
fvptc_cols <- names(ppm_raw)[grepl("^FVPTC\\d+$", names(ppm_raw))]
#- Bring in appropriate metadata to make prioritization decisions
fragment_quality_info <- read_excel("Data and Metadata Files/primary_data.xlsx", sheet = "lib.subject.qsummary") %>%
arrange(cas) %>%
select(name_sub_lib_id, iMean) %>%
rename(iMean_tumors = iMean)
#- Fully join and gather metadata for detection and means
SF4_ppm_fulltable <- ppm_raw_ctrl %>%
inner_join(ppm_raw, by = "name_sub_lib_id") %>%
left_join(cas_key_2, by = "name_sub_lib_id") %>%
filter(!cas %in% endog_cas) %>%
  arrange(Name) %>%
  rowwise() %>%
  mutate(
    half_min = (min(c_across(where(is.numeric)), na.rm = TRUE)) / 2,
    max_value = max(c_across(where(is.numeric)), na.rm = TRUE)
  ) %>%
  mutate(
    pct_det_ctrl = sum(!is.na(c_across(T001:T009))) / length(c_across(T001:T009)),
    pct_det_tumor = sum(!is.na(c_across(F1:F20))) / length(c_across(F1:F20))
  ) %>% #Filter to ideal fragment based on pct_det_tumor, then ctrl, then iMean_tumors below
left_join(fragment_quality_info, by = "name_sub_lib_id") %>%
arrange(desc(pct_det_tumor), desc(pct_det_ctrl), desc(iMean_tumors)) %>%
group_by(cas) %>%
slice_head(n = 1) %>%
ungroup() %>%
select(-iMean_tumors) %>%
mutate(across(
  T001:F20,
  ~ if_else(is.na(.x), half_min, .x)
)) %>%
rowwise() %>%
  mutate(
    mean_ctrl = mean(c_across(T001:T009), na.rm = TRUE),
    mean_tumor = mean(c_across(F1:F20), na.rm = TRUE),
    mean_PTC = mean(c_across(all_of(ptc_cols)), na.rm = TRUE),
    mean_FTC = mean(c_across(all_of(ftc_cols)), na.rm = TRUE),
    mean_FVPTC = mean(c_across(all_of(fvptc_cols)), na.rm = TRUE)
  ) %>%
ungroup() %>%
left_join(short_name %>% select(cas, annot_ident,IARC_Group), by = "cas") %>%
select(Name, cas, annot_ident, IARC_Group, pct_det_ctrl, pct_det_tumor, half_min, max_value, name_sub_lib_id, mean_ctrl:mean_FVPTC,T001:F20) %>%
arrange(cas)
#- Format and export to excel
format_ppb_value <- function(x) {
if (is.na(x)) {
  return(NA_character_)
}
if (x < 1) {
  formatC(x, format = "e", digits = 1)
} else if (x >= 1000) {
  formatC(x, format = "e", digits = 1)
} else {
  round(x)
}
}
PPB_version <- SF4_ppm_fulltable %>%
select(Name:mean_tumor,mean_FTC,mean_FVPTC,mean_PTC,-c(T001:F20,name_sub_lib_id)) %>%
mutate(across(half_min:mean_FVPTC, ~ .x * 1000)) %>%
  mutate(across(
    which(names(.) == "half_min"):ncol(.),
    ~ sapply(.x, format_ppb_value)))
write_xlsx(PPB_version, "SF4_ppb_fulltable.xlsx")
#- Create a pivoted version for subsequent analysis
full_joiner <- SF4_ppm_fulltable %>%
select(name_sub_lib_id,T001:F20) %>%
pivot_longer(
  cols = -name_sub_lib_id,
  names_to = "sample_ID",
  values_to = "value"
) %>%
mutate(
  variant = case_when(
    str_starts(sample_ID, "T") ~ "Ctrl",
    str_starts(sample_ID, "P") ~ "PTC",
    str_starts(sample_ID, "FVPTC") ~ "FV_PTC",
    str_starts(sample_ID, "F") ~ "FTC"
  )
) %>%
pivot_wider(names_from = name_sub_lib_id, values_from = value) %>%
mutate(
  tumor_vs_ctrl = ifelse(variant == "Ctrl", "Control", "Tumor"))
#+ Create final master variant table (ST3/4)
#- Create version of above SF4 ppm table except no filtering
SF4_ppb_inclusive <- ppm_raw_ctrl %>%
inner_join(ppm_raw, by = "name_sub_lib_id") %>%
left_join(cas_key_2, by = "name_sub_lib_id") %>%
filter(!cas %in% endog_cas) %>%
arrange(Name) %>%
rowwise() %>%
mutate(
  half_min = (min(c_across(where(is.numeric)), na.rm = TRUE)) / 2,
  max_value = max(c_across(where(is.numeric)), na.rm = TRUE)
) %>%
mutate(
  pct_det_ctrl = sum(!is.na(c_across(T001:T009))) / length(c_across(T001:T009)),
  pct_det_tumor = sum(!is.na(c_across(F1:F20))) / length(c_across(F1:F20))
) %>% 
mutate(across(
  T001:F20,
  ~ if_else(is.na(.x), half_min, .x)
)) %>%
rowwise() %>%
mutate(
  mean_ctrl_PPB = mean(c_across(T001:T009), na.rm = TRUE),
  mean_tumor_PPB = mean(c_across(F1:F20), na.rm = TRUE),
  mean_PTC_PPB = mean(c_across(all_of(ptc_cols)), na.rm = TRUE),
  mean_FTC_PPB = mean(c_across(all_of(ftc_cols)), na.rm = TRUE),
  mean_FVPTC_PPB = mean(c_across(all_of(fvptc_cols)), na.rm = TRUE)
) %>%
ungroup() %>%
arrange(cas) %>% 
select(name_sub_lib_id,pct_det_ctrl,pct_det_tumor,mean_ctrl_PPB:mean_FVPTC_PPB) %>%
mutate(across(mean_ctrl_PPB:mean_FVPTC_PPB, ~ .x * 1000)) # convert to ppb
#- Join with master table, pare down to final version
MT_final <- MTi %>%
left_join(SF4_ppb_inclusive, by = "name_sub_lib_id") %>%
select(short_name, cas, mode, annot_ident, p_value, highest, FTC:PTC, FTC_let:PTC_let, Carcinogenicity:GHS_var_diff_only, Potential_EDC, usage_class:Table_Class, pct_det_ctrl:mean_FVPTC_PPB)
#- Pull superscript letters off and place on +/- columns
#_Write function to do this
move_superscript <- function(main, sd_raw) {
  superscript_pattern <- "[\\u1d43-\\u1d4d\\u02b0-\\u02b8\\u1d62-\\u1d6b\\u2070-\\u209f\\u2020-\\u2021]+"
  # Extract superscript from the main value
  superscript <- str_extract(main, superscript_pattern)
  # Remove superscript from main
  main_clean <- str_remove(main, superscript_pattern)
  # Format SD to 2 decimal places
  sd_clean <- format(round(as.numeric(sd_raw), 2), nsmall = 2)
  # Recombine
  paste0(main_clean, " ± ", sd_clean, ifelse(is.na(superscript), "", superscript))
}
#_Apply function to the relevant columns
MT_export <- MT_final %>%
  mutate(
    FTC_let = case_when(
      mode == "quantitative" & !is.na(FTC_let) & str_detect(FTC, "±") ~
        move_superscript(FTC_let, str_trim(str_extract(FTC, "(?<=±).*"))),
      mode == "qualitative" & is.na(FTC_let) ~ FTC,
      TRUE ~ FTC_let
    ),
    FV_PTC_let = case_when(
      mode == "quantitative" & !is.na(FV_PTC_let) & str_detect(FV_PTC, "±") ~
        move_superscript(FV_PTC_let, str_trim(str_extract(FV_PTC, "(?<=±).*"))),
      mode == "qualitative" & is.na(FV_PTC_let) ~ FV_PTC,
      TRUE ~ FV_PTC_let
    ),
    PTC_let = case_when(
      mode == "quantitative" & !is.na(PTC_let) & str_detect(PTC, "±") ~
        move_superscript(PTC_let, str_trim(str_extract(PTC, "(?<=±).*"))),
      mode == "qualitative" & is.na(PTC_let) ~ PTC,
      TRUE ~ PTC_let
    )
  ) %>%
  select(-c(FTC:PTC)) %>%
  rename(FTC = FTC_let, FV_PTC = FV_PTC_let, PTC = PTC_let)
#- Export
write_xlsx(MT_export, "MT_final_SF4.xlsx")