#* 10: Post Validation Cleaning
#+ 10.1: Subset to validated compounds
#- 10.1.1: Validated compounds
validated_compounds <- validation_check_files |>
  select(id, quality, cas, source, short_name) |>
  filter(quality %in% c(1, 2))
#- 10.1.2: Validated Variant Differences list
validated_variant <- validated_compounds |>
  filter(source %in% c("VD", "IARC and VD")) |>
  pull(id)
#- 10.1.2: Validated iarc
validated_iarc <- validated_compounds |>
  filter(source %in% c("IARC", "IARC and VD")) |>
  pull(id)
#- 10.1.2: Validated variant joiner tibble
validated_list <- validated_compounds |>
  select(id, quality)
#+ 10.2: Modify MT_final per validation results
#- 10.2.1: Subset MT_final to validated variant differences
MT_final <- MT_final_i |>
  filter(id %in% validated_variant) |>
  inner_join(validated_list, by = "id") |>
  relocate(quality, .before = annot_ident)
#- 10.2.2: Pull name_sub_lib_id list
MT_final_namesub_list <- MT_final |>
  pull(name_sub_lib_id)
#- 10.2.2: Pull cas list
MT_final_cas_list <- MT_final |>
  pull(cas)
#+ 10.3: Modify cadaver/tumor comparisons by validation results
#- 10.3.1: Subset controls to validated IARCs
IARC_controls <- IARC_controls_i |>
  filter(id %in% validated_iarc) |>
  left_join(validated_compounds |>
    select(id, short_name), by = "id") |>
  select(pct_NA_ctrl, id, short_name, name_sub_lib_id, iMean_ctrl, subid)
#- 10.3.2: Subset tumors to validated IARCs
IARC_tumors <- IARC_tumors_i |>
  filter(id %in% validated_iarc) |>
  left_join(validated_compounds |>
    select(id, short_name), by = "id") |>
  select(pct_NA_tumor = pct_NA, id, short_name, name_sub_lib_id, iMean_tumor = iMean, subid)
#- 10.3.4: Do a full join of the detection percentages of both
IARC_combined <- IARC_controls |>
  full_join(IARC_tumors, by = c("name_sub_lib_id", "id", "subid", "short_name"), suffix = c("_ctrl", "_tumor")) |>
  mutate(
    pct_detected_tumor = 1 - replace_na(pct_NA_tumor, 1),
    pct_detected_ctrl = 1 - replace_na(pct_NA_ctrl, 1)
  ) |>
  filter(pct_detected_tumor >= 0.7 & pct_detected_ctrl >= 0.7) |>
  arrange(short_name, desc(pct_detected_tumor)) |>
  select(name_sub_lib_id, short_name, id, subid, pct_detected_tumor, pct_detected_ctrl, iMean_ctrl, iMean_tumor) |>
  arrange(desc(pct_detected_tumor), desc(pct_detected_ctrl))
#- 10.3.5: Create list of those remaining IARCs
IARC_namesub_pull <- IARC_combined |>
  pull(name_sub_lib_id)
#- 10.3.2: Use these to complete full joiner (ppm ppb estimates)
full_joiner <- full_joiner_i |>
  select(sample_ID, variant, tumor_vs_ctrl, any_of(IARC_namesub_pull))
#+ 10.4: Modify the posthoc and quant sig to only validated compounds
#- 10.4.1: Subset quant sig to validated compounds
tumors_quant_sig <- tumors_quant_sig_i |>
  select(variant, any_of(MT_final_namesub_list))
#- 10.4.2: Subset posthoc table to validated compounds
posthoc_table_pvalues <- posthoc_table_pvalues_i |>
  filter(cas %in% MT_final_cas_list) |>
  arrange(p_value)
#- 10.4.4: Subset fisher results to validated
fisher_results <- fisher_results_i |>
  filter(name_sub_lib_id %in% MT_final_namesub_list)
