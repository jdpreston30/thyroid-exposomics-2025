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
        "ᵃ", ""
      ),
      ifelse(Potential_EDC == "Y", "ᵇ", ""),
      ifelse(has_asterisk, "ᶜ", "")
    ),
    FTC_let = coalesce(FTC_let, FTC),
    FV_PTC_let = coalesce(FV_PTC_let, FV_PTC),
    PTC_let = coalesce(PTC_let, PTC),
    p_value = sprintf("%.3f", p_value)
  ) |>
  mutate(short_name = str_replace(short_name, "NA([ᵃᵇᶜ]*)$", "\\1")) |>
  # Correct marker ordering (ensure abc order)
  mutate(
    short_name = case_when(
      str_detect(short_name, "ᵇᵃ") ~ str_replace(short_name, "ᵇᵃ", "ᵃᵇ"),
      str_detect(short_name, "ᶜᵃ") ~ str_replace(short_name, "ᶜᵃ", "ᵃᶜ"),
      str_detect(short_name, "ᶜᵇᵃ") ~ str_replace(short_name, "ᶜᵇᵃ", "ᵃᵇᶜ"),
      str_detect(short_name, "ᶜᵃᵇ") ~ str_replace(short_name, "ᶜᵃᵇ", "ᵃᵇᶜ"),
      str_detect(short_name, "ᵇᶜ") ~ str_replace(short_name, "ᵇᶜ", "ᵇᶜ"),
      str_detect(short_name, "ᵃᶜ") ~ str_replace(short_name, "ᵃᶜ", "ᵃᶜ"),
      TRUE ~ short_name
    )
  ) |>
  arrange(p_value) |>
  select(`Chemical Name` = short_name, `Usage Class (Type)` = Table_Class, FTC = FTC_let, `FV-PTC` = FV_PTC_let, PTC = PTC_let, `p-value` = p_value) |>
  mutate(across(c(FTC, `FV-PTC`, PTC), ~ str_replace_all(.x, c("ᵃ" = "ᵈ", "ᵇ" = "ᵉ", "ᶜ" = "ᶠ"))))
#- 15.3.2: Build Table 3 with function
table_3 <- build_table_3(
  data = table_3_tibble,
  export_path = "Outputs/Tables/T3.xlsx"
)
