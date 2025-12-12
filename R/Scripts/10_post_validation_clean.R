#+ 8.6: Modify MT_final per validation results
MT_final <- MT_final_i |>
  select()
#+ 8.7: Modify cadaver/tumor comparisons by validation results
cadaver_iarc_keep <- cadaver_iarc |>
  filter(cas == "95-53-4") |>
  pull(name_sub_lib_id)
#- 8.6.2: full_joiner_i
full_joiner <- full_joiner_i |>
  select(sample_ID, variant, tumor_vs_ctrl, any_of(cadaver_iarc_keep))
  