#* 10: Post Validation Cleaning
#+ 10.1: Subset to validated compounds
#- 10.1.1: Validated compounds
validated_compounds <- validation_check_files |>
  select(id, quality) |>
  filter(quality %in% c(1, 2))
#- 10.1.2: Validated list
validated_list <- validated_compounds |>
  pull(id)
#+ 10.2: Modify MT_final per validation results
MT_final <- MT_final_i |>
  filter(id %in% validated_list) |>
  inner_join(validated_compounds, by = "id") |>
  relocate(quality, .before = annot_ident)



#+ 10.3: Modify cadaver/tumor comparisons by validation results
cadaver_iarc_keep <- cadaver_iarc |>
  filter(cas == "95-53-4") |>
  pull(name_sub_lib_id)
#- 10.3.2: full_joiner_i
full_joiner <- full_joiner_i |>
  select(sample_ID, variant, tumor_vs_ctrl, any_of(cadaver_iarc_keep))
  