#* 11: Supplementary Tables
#+ 11.4: ST1- Chemical library (pivoted subid)
ST1_import
      #! In excel, then pared down and formatted, but reimporting here to double check the molecular formulas
#+ 11.7: ST4- Full ppm/ppb table for tumor-cadaver inner join
ST4 <- ppm_full_table |>
  select(Name:mean_tumor, mean_FTC, mean_FVPTC, mean_PTC, -c(T001:F20, name_sub_lib_id)) |>
  mutate(across(half_min:mean_FVPTC, ~ .x * 1000)) |>
  mutate(across(
    half_min:mean_PTC,
    ~ sapply(.x, format_ppb_value)
  ))
