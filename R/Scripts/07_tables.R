#* 7: Supplementary Tables
#+ 7.0: Make final master table, pare down to final version

#+ 7.1: ST1- Chemical library (pivoted subid)
ST1 <- read_excel(config$paths$primary_data, sheet = "library") |>
  filter(Disposition != "Endogenous") |>
  mutate(subid_col = paste0("mz", subid)) |>
  select(id, name, short_display_name, trt, monoisotopic, cas, formula, Disposition, subid_col, tmz) |>
  distinct() |>
  pivot_wider(
    names_from = subid_col,
    values_from = tmz
  ) |>
  arrange(cas)
      #! In excel, then pared down and formatted, but reimporting here to double check the molecular formulas
#+ 7.4: ST4- Full ppm/ppb table for tumor-cadaver inner join
ST4 <- ppm_full_table |>
  select(Name:mean_tumor, mean_FTC, mean_FVPTC, mean_PTC, -c(T001:F20, name_sub_lib_id)) |>
  mutate(across(half_min:mean_FVPTC, ~ .x * 1000)) |>
  mutate(across(
    half_min:mean_PTC,
    ~ sapply(.x, format_ppb_value)
  ))
