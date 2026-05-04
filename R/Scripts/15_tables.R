#* 15: Tables
#+ 15.1: Build Table 1 (with function); export
table_1 <- build_table_1(
  data = demographic_table,
  export_path = "Outputs/Tables/T1.xlsx"
)
#+ 15.2: Build Table 2 (with function); export
table_2 <- build_table_2(
  data = feature_metadata,
  header_col = "Table_Header",
  class_col = "Table_Class",
  subclass_col = "Table_Subclass",
  export_path = "Outputs/Tables/T2.xlsx"
)
#+ 15.3: Build Table 3 (with function); export
#- 15.3.1: Create table columns and format
table_3_tibble <- MT_final |>
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
    has_asterisk = str_detect(short_name, "\\*"),
    short_name = str_replace_all(short_name, "\\*", ""),
    short_name = paste0(
      short_name,
      if_else(
        Carcinogenicity %in% c("Likely Carcinogen", "Possible Carcinogen", "Known Carcinogen"),
        "\u2020", ""
      ),
      ifelse(Potential_EDC == "Y", "\u2021", ""),
      ifelse(has_asterisk, "\u00b6", "")
    ),
    FTC_let = coalesce(FTC_let, FTC),
    FV_PTC_let = coalesce(FV_PTC_let, FV_PTC),
    PTC_let = coalesce(PTC_let, PTC),
    p_value = sprintf("%.3f", p_value)
  ) |>
  mutate(short_name = str_replace(short_name, "NA([\u2020\u2021\u00b6]*)$", "\\1")) |>
  # Correct marker ordering (ensure † ‡ ¶ order)
  mutate(
    short_name = case_when(
      str_detect(short_name, "\u2021\u2020") ~ str_replace(short_name, "\u2021\u2020", "\u2020\u2021"),
      str_detect(short_name, "\u00b6\u2020") ~ str_replace(short_name, "\u00b6\u2020", "\u2020\u00b6"),
      str_detect(short_name, "\u00b6\u2021\u2020") ~ str_replace(short_name, "\u00b6\u2021\u2020", "\u2020\u2021\u00b6"),
      str_detect(short_name, "\u00b6\u2020\u2021") ~ str_replace(short_name, "\u00b6\u2020\u2021", "\u2020\u2021\u00b6"),
      str_detect(short_name, "\u2021\u00b6") ~ str_replace(short_name, "\u2021\u00b6", "\u2021\u00b6"),
      str_detect(short_name, "\u2020\u00b6") ~ str_replace(short_name, "\u2020\u00b6", "\u2020\u00b6"),
      TRUE ~ short_name
    )
  ) |>
  arrange(p_value) |>
  select(`Chemical Name` = short_name, `Usage Class (Type)` = Table_Class, FTC = FTC_let, `IEFVPTC` = FV_PTC_let, PTC = PTC_let, `p-value` = p_value)
#- 15.3.2: Build Table 3 with function
table_3 <- build_table_3(
  data = table_3_tibble,
  export_path = "Outputs/Tables/T3.xlsx"
)
#+ 15.4: Build Table 4 (with function); export
table_4 <- build_table_4(
  ppm_full_table = ppm_full_table,
  ST1_tibble = ST1_tibble,
  literature_ST3 = literature_ST3,
  validation_check_files_unfiltered = validation_check_files_unfiltered,
  export_path = "Outputs/Tables/T4.xlsx"
)

								