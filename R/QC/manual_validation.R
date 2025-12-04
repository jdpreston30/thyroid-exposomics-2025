#* Manual Spectral Validation QC
#+ Read in sequences
#- Tumors
tumor_files <- read_excel(config$paths$primary_data, sheet = "tumors_sequence") |>
  select(file, ID, replicate) |>
  group_by(ID) |>
  summarize(files = paste(file, collapse = ", "), .groups = "drop")
#- Cadaver Thyroid
cadaver_files <- read_csv(config$paths$cadaver_sequence) |>
  select(file = Filename, ID = `Sample ID`, sample_name = `Sample name`) |>
  filter(str_detect(ID, "^[Tt]-\\d+")) |>
  mutate(
    replicate = as.numeric(str_extract(sample_name, "(?<=_)\\d+$")),
    ID = str_replace(ID, "-", "")  # Remove hyphen to match column names (T001, T002, etc.)
  ) |>
  select(-sample_name) |>
  group_by(ID) |>
  summarize(files = paste(file, collapse = ", "), .groups = "drop")
#+ Pull features that need validation
#- Subset to relevant features (top 5 quant)
top_5_quant_cas <- MT_final_i |>
  filter(mode == "quantitative") |>
  arrange(p_value) |>
  slice_head(n = 5) |>
  pull(cas)
#- Subset to relevant IARCs (group 1)
iarc_1_validate <- tumor_raw |>
  filter(name_sub_lib_id %in% IARC_tumors_ctrl_filtered) |>
  pull(cas)
#+ Create tumor validation table
tumor_top5_iarc <- tumor_raw |>
  filter(cas %in% top_5_quant_cas | cas %in% iarc_1_validate) |>
  mutate(
    validate_for = case_when(
      cas %in% top_5_quant_cas & cas %in% iarc_1_validate ~ "quant and iarc",
      cas %in% top_5_quant_cas ~ "quant",
      cas %in% iarc_1_validate ~ "iarc",
      TRUE ~ NA_character_
    )
  ) |>
  select(cas, short_display_name, lib, mzMed, rtMin, rtMed, rtMax, tmz, trt, id, subid, id_subid, name_sub_lib_id, validate_for) |>
  bind_cols(
    tumor_raw |>
      filter(cas %in% top_5_quant_cas | cas %in% iarc_1_validate) |>
      select(matches("^F[0-9]+$")) |>
      apply(1, function(row) {
        if (all(is.na(row))) {
          tibble(check1 = NA_character_, check2 = NA_character_, check3 = NA_character_,
                 check4 = NA_character_, check5 = NA_character_)
        } else {
          top_indices <- order(row, decreasing = TRUE, na.last = TRUE)
          top_5 <- head(names(row)[top_indices][!is.na(row[top_indices])], 5)
          top_5 <- c(top_5, rep(NA_character_, 5 - length(top_5)))
          tibble(check1 = top_5[1], check2 = top_5[2], check3 = top_5[3],
                 check4 = top_5[4], check5 = top_5[5])
        }
      }) |>
      bind_rows()
  ) |>
  left_join(tumor_files |> select(ID, file1 = files), by = c("check1" = "ID")) |>
  left_join(tumor_files |> select(ID, file2 = files), by = c("check2" = "ID")) |>
  left_join(tumor_files |> select(ID, file3 = files), by = c("check3" = "ID")) |>
  left_join(tumor_files |> select(ID, file4 = files), by = c("check4" = "ID")) |>
  left_join(tumor_files |> select(ID, file5 = files), by = c("check5" = "ID")) |>
  relocate(file1, .after = check1) |>
  relocate(file2, .after = check2) |>
  relocate(file3, .after = check3) |>
  relocate(file4, .after = check4) |>
  relocate(file5, .after = check5) |>
  arrange(cas, id_subid) |>
  arrange(desc(validate_for))
#+ Create controls validation table
cadaver_top5_iarc <- IARC_controls_i |>
  filter(cas %in% iarc_1_validate) |>
  left_join(feature_metadata |> select(cas, Short_display_name), by = "cas") |>
  mutate(across(matches("^T[0-9]+$"), as.numeric)) |>
  mutate(validate_for = "iarc") |>
  select(cas, short_display_name = Short_display_name, lib, mzMed, rtMin, rtMed, rtMax, tmz, trt, id, subid, id_subid, name_sub_lib_id, validate_for) |>
  bind_cols(
    IARC_controls_i |>
      filter(cas %in% iarc_1_validate) |>
      mutate(across(matches("^T[0-9]+$"), as.numeric)) |>
      select(matches("^T[0-9]+$")) |>
      apply(1, function(row) {
        if (all(is.na(row))) {
          tibble(check1 = NA_character_, check2 = NA_character_, check3 = NA_character_,
                 check4 = NA_character_, check5 = NA_character_)
        } else {
          top_indices <- order(row, decreasing = TRUE, na.last = TRUE)
          top_5 <- head(names(row)[top_indices][!is.na(row[top_indices])], 5)
          top_5 <- c(top_5, rep(NA_character_, 5 - length(top_5)))
          tibble(check1 = top_5[1], check2 = top_5[2], check3 = top_5[3],
                 check4 = top_5[4], check5 = top_5[5])
        }
      }) |>
      bind_rows()
  ) |>
  left_join(cadaver_files |> select(ID, file1 = files), by = c("check1" = "ID")) |>
  left_join(cadaver_files |> select(ID, file2 = files), by = c("check2" = "ID")) |>
  left_join(cadaver_files |> select(ID, file3 = files), by = c("check3" = "ID")) |>
  left_join(cadaver_files |> select(ID, file4 = files), by = c("check4" = "ID")) |>
  left_join(cadaver_files |> select(ID, file5 = files), by = c("check5" = "ID")) |>
  relocate(file1, .after = check1) |>
  relocate(file2, .after = check2) |>
  relocate(file3, .after = check3) |>
  relocate(file4, .after = check4) |>
  relocate(file5, .after = check5) |>
  arrange(cas, id_subid) |>
  arrange(desc(validate_for))
#+ Create mz reference table
#- Subset to metadata for tumors
tumor_features <- tumor_top5_iarc |>
  select(cas, short_display_name, lib, id, subid, tmz, trt)
#- Subset to metadata for cadavers
cadaver_features <- cadaver_top5_iarc |>
  select(cas, short_display_name, lib, id, subid, tmz, trt)
#- Bind rows and get distinct
mz_reference_table <- bind_rows(tumor_features, cadaver_features) |>
  distinct() |>
  arrange(short_display_name, id, subid) |>
  # Pivot to create mz0, mz1, mz2, etc columns
  pivot_wider(
    id_cols = c(cas, short_display_name, lib, id, trt),
    names_from = subid,
    values_from = tmz,
    names_prefix = "mz",
    values_fill = NA_real_,
    names_sort = TRUE
  ) |>
  arrange(short_display_name, id)
#+ Export to Excel with formatting
#- Create workbook
wb <- createWorkbook()
#- Add tumor sheet
addWorksheet(wb, "Tumor_Validation")
writeData(wb, "Tumor_Validation", tumor_top5_iarc, startRow = 1, startCol = 1)
#- Add cadaver sheet
addWorksheet(wb, "Cadaver_Validation")
writeData(wb, "Cadaver_Validation", cadaver_top5_iarc, startRow = 1, startCol = 1)
#- Add mz reference sheet
addWorksheet(wb, "MZ_Reference")
writeData(wb, "MZ_Reference", mz_reference_table, startRow = 1, startCol = 1)
#- Style for header row (bold)
header_style <- createStyle(textDecoration = "bold", border = "TopBottomLeftRight", borderStyle = "thin")
#- Style for data cells (borders only)
data_style <- createStyle(border = "TopBottomLeftRight", borderStyle = "thin")
#- Apply styles to tumor sheet
addStyle(wb, "Tumor_Validation", header_style, rows = 1, cols = 1:ncol(tumor_top5_iarc), gridExpand = TRUE)
addStyle(wb, "Tumor_Validation", data_style, rows = 2:(nrow(tumor_top5_iarc) + 1), cols = 1:ncol(tumor_top5_iarc), gridExpand = TRUE)
#- Apply styles to cadaver sheet
addStyle(wb, "Cadaver_Validation", header_style, rows = 1, cols = 1:ncol(cadaver_top5_iarc), gridExpand = TRUE)
addStyle(wb, "Cadaver_Validation", data_style, rows = 2:(nrow(cadaver_top5_iarc) + 1), cols = 1:ncol(cadaver_top5_iarc), gridExpand = TRUE)
#- Apply styles to mz reference sheet
addStyle(wb, "MZ_Reference", header_style, rows = 1, cols = 1:ncol(mz_reference_table), gridExpand = TRUE)
addStyle(wb, "MZ_Reference", data_style, rows = 2:(nrow(mz_reference_table) + 1), cols = 1:ncol(mz_reference_table), gridExpand = TRUE)
#- Save workbook
saveWorkbook(wb, "metadata_files/manual_validation.xlsx", overwrite = TRUE)

