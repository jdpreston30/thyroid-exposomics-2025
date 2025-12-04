
#- 6.5.5: Pull actual features out, transpose for graphing
IARC_comp <- full_joiner |>
  select(sample_ID, any_of(IARC_tumors_ctrl_filtered)) |>
  mutate(across(
    .cols = -c(sample_ID),
    .fns = ~ {
      col_min <- min(.x, na.rm = TRUE)
      replace(.x, is.na(.x), 0.5 * col_min)
    }
  )) |>
  pivot_longer(-sample_ID, names_to = "name_sub_lib_id", values_to = "value") |>
  pivot_wider(names_from = sample_ID, values_from = value) |>
  left_join(cas_key_2, by = "name_sub_lib_id") |>
  left_join(short_name |> select(cas, annot_ident), by = "cas") |>
  mutate(Name = if_else(
    annot_ident == "Annotation",
    paste0(Name, "\u2020"),
    Name
  )) |>
  select(Name, T001:F20)
#+ 6.6: IARC Stats (ttest on log transformed)
# ! Data not shown but listed in manuscript results
IARC_ttests <- full_joiner |>
  select(sample_ID, any_of(IARC_tumors_ctrl_filtered)) |>
  mutate(across(
    .cols = -c(sample_ID),
    .fns = ~ {
      col_min <- min(.x, na.rm = TRUE)
      replace(.x, is.na(.x), 0.5 * col_min)
    }
  )) |>
  mutate(across(where(is.numeric), ~ log2(. + 1e-6))) |>
  mutate(group = if_else(str_starts(sample_ID, "T00"), "Control", "Tumor")) |>
  pivot_longer(
    cols = where(is.numeric) & !any_of("group"),
    names_to = "chemical", values_to = "value"
  ) |>
  group_by(chemical) |>
  summarise(
    p_value = t.test(value ~ group)$p.value,
    mean_control = mean(value[group == "Control"]),
    mean_tumor = mean(value[group == "Tumor"]),
    .groups = "drop"
  ) |>
  arrange(p_value)
#+ 6.7: Make MT Final with QC
#---------
# make mt final with removed if needed