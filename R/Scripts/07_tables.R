#* 7: Supplementary Tables
#+ 7.1: Build Table 1 (with function); export
table_1 <- build_table_1(
  data = demographic_table,
  export_path = "Outputs/Tables/T1.xlsx"
)
#+ 7.2: Build Table 2 (with function); export
table_2 <- build_table_2(
  data = feature_metadata,
  header_col = "Table_Header",
  class_col = "Table_Class",
  subclass_col = "Table_Subclass",
  export_path = "Outputs/Tables/T2.xlsx"
)
#+ 7.3: Build Table 3 (with function); export
#- 7.3.1: Create table columns and format
table_3 <- MT_final |>
  mutate(
    Table_Class = case_when(
      Table_Class == "Insecticides and Pesticides" ~ "Insecticide/Pesticide",
      str_detect(Table_Class, "Dye intermediates") ~ "Dye intermediate",
      str_detect(Table_Class, "Chemical Synthesis Intermediates") ~ "Chemical Synthesis Intermediate",
      str_detect(Table_Class, "Carcinogenic Research Chemicals") ~ "Carcinogenic Research Chemical",
      str_detect(Table_Class, "Combustion Byproducts") ~ "Combustion Byproduct",
      str_detect(Table_Class, "Side-Reaction Byproducts") ~ "Side-Reaction Byproduct",
      str_detect(Table_Class, "Fungicides") ~ "Fungicide",
      str_detect(Table_Class, "Herbicides") ~ "Herbicide",
      str_detect(Table_Class, "Disinfectant Breakdown Products") ~ "Disinfectant Breakdown Product",
      str_detect(Table_Class, "Flavoring or Fragrance Agents") ~ "Flavoring or Fragrance Agent",
      str_detect(Table_Class, "Plasticizers and Plastic Additives") ~ "Plasticizer",
      str_detect(Table_Class, "Preservatives") ~ "Preservative",
      str_detect(Table_Class, "Flame Retardants") ~ "Flame Retardant",
      str_detect(Table_Class, "Plant Growth Regulators") ~ "Plant Growth Regulator",
      str_detect(Table_Class, "Humectants") ~ "Humectant",
      TRUE ~ Table_Class # Leave PFAS, Organic UV Filters, etc. unchanged
    )
  ) |>
  left_join(feature_metadata |> select(cas, Table_Subclass, Table_Qualifier), by = "cas") |>
  mutate(
    Table_Class = ifelse(
      !is.na(Table_Subclass) & Table_Subclass != "",
      paste0(Table_Class, " (", Table_Subclass, ")"),
      Table_Class
    )
  ) |>
  mutate(
    short_name = paste0(
      short_name,
      if_else(annot_ident == "Annotation", "*", ""),
      if_else(
        Carcinogenicity %in% c("Likely Carcinogen", "Possible Carcinogen", "Known Carcinogen"),
        "\u2020", ""
      ),
      ifelse(Potential_EDC == "Y", "\u2021", "") # <-- FIXED HERE
    ),
    FTC_let = coalesce(FTC_let, FTC),
    FV_PTC_let = coalesce(FV_PTC_let, FV_PTC),
    PTC_let = coalesce(PTC_let, PTC),
    p_value = round(p_value, 3)
  ) |>
  mutate(short_name = str_replace(short_name, "NA$", "")) |>
  select(short_name, cas, annot_ident, Table_Class, FTC_let:PTC_let, p_value)
#- 7.3.2: Build Table 3 with function

#+ 7.4: ST1- Chemical library (pivoted subid)
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
#+ 7.7: ST4- Full ppm/ppb table for tumor-cadaver inner join
ST4 <- ppm_full_table |>
  select(Name:mean_tumor, mean_FTC, mean_FVPTC, mean_PTC, -c(T001:F20, name_sub_lib_id)) |>
  mutate(across(half_min:mean_FVPTC, ~ .x * 1000)) |>
  mutate(across(
    half_min:mean_PTC,
    ~ sapply(.x, format_ppb_value)
  ))
