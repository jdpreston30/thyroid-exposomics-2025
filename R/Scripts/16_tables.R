#* 16: Tables
#+ 16.1: Pre-build ST1_tibble (needed by build_table_4)
#- 16.1.1: Get list of Y/N expanded fragments
expanded_chemicals <- expanded_validation |>
  select(id, expanded) |>
  filter(expanded == "Y") |>
  pull(id)
#- 16.1.2: Prepare ST1_tibble
ST1_tibble <- ST1_import |>
  mutate(
    base_num = as.numeric(str_extract(id, "\\d+(?=_|$)")),
    suffix_num = as.numeric(str_extract(id, "(?<=_)\\d+"))
  ) |>
  mutate(suffix_num = replace_na(suffix_num, 0)) |>
  group_by(cas) |>
  arrange(base_num, suffix_num, .by_group = TRUE) |>
  mutate(
    is_first_cas = row_number() == 1,
    cas_group_size = n()
  ) |>
  ungroup() |>
  # Sort alphabetically by name for global order
  arrange(name, base_num, suffix_num) |>
  # Recalculate grouping after alphabetical sort
  group_by(cas) |>
  mutate(
    is_first_cas = row_number() == 1,
    cas_group_size = n()
  ) |>
  ungroup() |>
  # Replace values with dashes for non-first occurrences of same CAS
  mutate(
    name = if_else(!is_first_cas, "-", name),
    cas = if_else(!is_first_cas, "-", cas),
    monoisotopic = if_else(!is_first_cas, "-", as.character(monoisotopic))
  ) |>
  # Replace asterisks with dagger for endogenous chemicals
  mutate(name = gsub("\\*", "\u2020", name)) |>
  # Add superscript b for expanded chemicals
  mutate(
    name = if_else(
      id %in% expanded_chemicals & name != "-",
      paste0(name, "\u2021"),
      name
    )
  ) |>
  # Final cleanup - add Index column
  mutate(Index = row_number(), .before = 1) |>
  select(Index, id, name, cas, monoisotopic, trt, starts_with("mz")) |>
  rename(
    `Library ID` = id,
    Name = name,
    CAS = cas,
    `Monoisotopic Mass` = monoisotopic,
    `Target RT (min)` = trt
  )
#+ 16.2: Build Table 1 (with TernTables); Export
table1 <- ternG(
  clinical_data,
  exclude_vars = c("Patient_ID", "year"),
  group_var = "Variant",
  output_docx = "Outputs/Tables/T1.docx",
  methods_doc = FALSE,
  print_normality = FALSE,
  show_test = FALSE,
  show_p = FALSE,
  open_doc = FALSE,
  citation = FALSE,
  font_family = "Times New Roman",
  table_font_size = 10.5,
  abbreviation_footnote = "Abbreviations: FTC = follicular thyroid carcinoma; IEFVPTC = invasive encapsulated follicular variant of papillary thyroid carcinoma; M = metastasis; N = lymph node; PTC = papillary thyroid carcinoma; T = tumor",
  category_start = c(
    "Demographics" = "Sex",
    "Staging" = "T Category"
  ))
#+ 16.3: Build Table 2 (with function); Export
table_2 <- build_table_2(
  data = feature_metadata,
  header_col = "Table_Header",
  class_col = "Table_Class",
  subclass_col = "Table_Subclass",
  export_path = "Outputs/Tables/T2.xlsx"
)
#+ 16.4: Build Table 3 (with function); Export
#- 16.4.1: Create table columns and format
table_3_tibble <- MT_final |>
  mutate(
    Table_Class = case_when(
      Table_Class == "Insecticides and Pesticides" ~ "Insecticide/Pesticide",
      str_detect(Table_Class, "Dye Intermediates") ~ "Dye Intermediate",
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
#! Keep the numeric P for sorting: sprintf() coerces to character, and arranging the string puts 3-decimal ties (e.g. the four 0.039 rows) in arbitrary order rather than true ascending order.
    p_sort = p_value,
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
  arrange(p_sort) |>
  select(`Chemical Name` = short_name, `Usage Class (Type)` = Table_Class, FTC = FTC_let, `IEFVPTC` = FV_PTC_let, PTC = PTC_let, `P` = p_value)
#- 16.4.2: Build table 3 with function
table_3 <- build_table_3(
  data = table_3_tibble,
  export_path = "Outputs/Tables/T3.xlsx"
)
#+ 16.5: Build Table 4 (with function); Export
table_4 <- build_table_4(
  ppm_full_table = ppm_full_table,
  ST1_tibble = ST1_tibble,
  literature_ST3 = literature_ST3,
  literature_long = literature_long,
  validation_check_files_unfiltered = validation_check_files_unfiltered,
  export_path = "Outputs/Tables/T4.xlsx"
)

								
