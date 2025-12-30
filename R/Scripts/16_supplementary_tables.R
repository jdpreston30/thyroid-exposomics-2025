#* 11: Supplementary Tables
#+ 11.4: ST1- Chemical library (pivoted subid)
ST1_import
      #! In excel, then pared down and formatted, but reimporting here to double check the molecular formulas

      #! •	SUPPLEMENTARY TABLE 1. The full library of xenobiotic chemicals employed for chemical identification. The library of 710 confirmed xenobiotic chemicals employed for chemical identification. All chemicals were present in pooled reference plasma at a concentration of 0·5 ng/mL. There are 710 total unique chemicals (i.e., unique CAS numbers), but for some chemicals, there are multiple fragments from different standards used for identification, thus resulting in 892 total rows in the table. The individual mz columns indicate typical fragments observed for the given chemical.
#+ 11.7: ST4- Full ppm/ppb table for tumor-cadaver inner join
ST4 <- ppm_full_table |>
  select(Name:mean_tumor, mean_FTC, mean_FVPTC, mean_PTC, -c(T001:F20, name_sub_lib_id)) |>
  mutate(across(half_min:mean_FVPTC, ~ .x * 1000)) |>
  mutate(across(
    half_min:mean_PTC,
    ~ sapply(.x, format_ppb_value)
  ))
