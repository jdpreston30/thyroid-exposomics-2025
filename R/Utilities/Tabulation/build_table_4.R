#' Build Table 4 (Observed Concentrations vs. Literature Values)
#'
#' Constructs and exports the supplementary table comparing tumor/non-cancer
#' thyroid concentrations to literature values for select IARC-classified chemicals.
#' Encapsulates all data wrangling (steps 16.3.0-16.3.9) and Excel export.
#'
#' @param ppm_full_table ppm data table from the pipeline
#' @param ST1_tibble ST1 tibble (used to join chemical names)
#' @param literature_ST3 Literature comparison data (from primary_data.xlsx, sheet "literature_comp_pared")
#' @param validation_check_files_unfiltered Unfiltered validation check files
#' @param export_path Path to export the Excel file
#'
#' @return The final ST3 tibble (hierarchical, formatted)
#'
build_table_4 <- function(ppm_full_table, ST1_tibble, literature_ST3,
                           validation_check_files_unfiltered, export_path) {
  library(openxlsx)
#! PCB-138 (35065-28-2) and PCB-153 (35065-27-1) are legitimately quantified and stay in Table 4. An older standard list flags them as mixture-sourced and therefore uncalibratable; that list is out of date. The preparation actually used was a newer 10 ug/mL per-congener standard, which propagates correctly through the BP1 dilution to the same reference concentration as every other EC (nominal 0.5 ng/mL, corrected to 0.47 by the 0.94 factor applied at import in 00c_FTs.R). Their ppb estimates rest on the same calibration as everything else in this table -- do not exclude them.
#_ Vector of selected CAS numbers
selected_ST3_cas <- c(
  "83-32-9", "208-96-8", "120-12-7", "56-55-3", "205-99-2", "207-08-9", "50-32-8", "218-01-9",
  "53-70-3", "206-44-0", "91-20-3", "85-01-8", "129-00-0", "91-59-8", "92-67-1", "92-87-5",
  "95-53-4", "32598-14-4", "35065-28-2", "35065-27-1", "52663-74-8", "35065-29-3", "52663-68-0"
)
#_ Create initial tibble
ST3_tibble_i <- ppm_full_table |>
  filter(!is.na(IARC_Group)) |>
  select(CAS = cas, everything(), -Name) |>
  mutate(
    across(c(mean_ctrl, mean_tumor, half_min, max_value),
           ~ . * 10^3,
           .names = "{.col}")
  ) |>
  mutate(
    across(c(mean_ctrl, mean_tumor),
           ~ if_else(. < 1, "< 1", format(round(.), big.mark = ",")),
           .names = "{.col}"),
    half_min = if_else(half_min < 1, "< 1", format(round(half_min), big.mark = ",")),
    max_value = if_else(max_value < 1, "< 1", format(round(max_value), big.mark = ","))
  ) |>
  left_join(ST1_tibble |> select(Name, CAS), by = "CAS") |>
  mutate(range = paste0(trimws(half_min), " - ", trimws(max_value))) |>
  group_by(CAS) |>
  mutate(
    best_fragment = if_else(pct_det_tumor == max(pct_det_tumor), "★", ""),
    tumor_ctrl_det = paste0(round(pct_det_tumor * 100), "/", round(pct_det_ctrl * 100))
  ) |>
  ungroup()
#_ Inspection tibble: filter to selected CAS and detection threshold
ST3_tibble_inspect <- ST3_tibble_i |>
  select(name_sub_lib_id, best_fragment, tumor_ctrl_det, range, IARC_Group, CAS, everything()) |>
  filter(CAS %in% selected_ST3_cas) |>
  filter(!(pct_det_tumor < 0.5 & pct_det_ctrl < 0.5)) |>
  group_by(CAS) |>
  mutate(Tie = if_else(sum(best_fragment == "★") > 1, "Y", "")) |>
  ungroup() |>
  arrange(name_sub_lib_id)
#_ Tiebreaker fragments (manual inspection)
tiebreak_ST3_fragments <- c(
  "Pyrene_3_BP3.GC2_CP3078",
  "Fluoranthene_3_BP3.GC2_CP3057",
  "Chrysene_0_BP3.GC2_CP3044"
)
#_ Apply tiebreaker
ST3_tibble_2 <- ST3_tibble_inspect |>
  filter(best_fragment == "★") |>
  filter((Tie == "Y" & name_sub_lib_id %in% tiebreak_ST3_fragments) | Tie == "")
#_ Validation failure check (informational)
ST3_fail_check <- validation_check_files_unfiltered |>
  select(id, cas, short_name, quality, state, modification) |>
  arrange(quality) |>
  arrange(short_name) |>
  filter(cas %in% ST3_tibble_2$CAS)
#_ Failed chemicals vector
failed_ST3 <- c(
  "50-32-8", # Benzo[a]pyrene
  "207-08-9", # Benzo(k)fluoranthene
  "92-87-5"  # Benzidine
)
#_ Build flat tibble with literature values
# citation letters (reading order): ᵃ=mlyczynska2023, ᵇ=wang2010, ᶜ=riffelmann1995, ᵈ=cdc2024, ᵉ=maier2022
ref_map <- c(
  mlyczynska2023 = "\u1d43",  # ᵃ
  maier2022      = "\u1d49",  # ᵉ
  riffelmann1995 = "\u1d9c",  # ᶜ
  cdc2024        = "\u1d48",  # ᵈ
  wang2010       = "\u1d47"   # ᵇ
)
resolve_ref <- function(ref) {
  dplyr::coalesce(ref_map[ref], ref)
}
ST3_tibble_flat <- ST3_tibble_2 |>
  filter(!CAS %in% failed_ST3) |>
  left_join(
    literature_ST3 |> select(CAS, `Usage Class`, AT_manuscript, AT_ref,
                              plasma_manuscript, plasma_ref, urine_manuscript, urine_ref),
    by = "CAS"
  ) |>
  mutate(
    Name = gsub("[\u2020\u2021]", "", Name),
    Name = gsub("\\*", "", Name),
    `Adipose Tissue (ppb)` = case_when(
      !is.na(AT_manuscript) & !is.na(AT_ref) ~
        paste0(AT_manuscript, resolve_ref(AT_ref)),
      !is.na(AT_manuscript) ~ AT_manuscript,
      TRUE ~ NA_character_
    ),
    `Urine (ppb)` = case_when(
      !is.na(urine_manuscript) & !is.na(urine_ref) ~
        paste0(urine_manuscript, resolve_ref(urine_ref)),
      !is.na(urine_manuscript) ~ urine_manuscript,
      TRUE ~ NA_character_
    ),
    `Serum/Plasma (ppb)` = case_when(
      !is.na(plasma_manuscript) & !is.na(plasma_ref) ~
        paste0(plasma_manuscript, resolve_ref(plasma_ref)),
      !is.na(plasma_manuscript) ~ plasma_manuscript,
      TRUE ~ NA_character_
    ),
    across(
      c(`Adipose Tissue (ppb)`, `Serum/Plasma (ppb)`),
      ~ str_replace_all(.x, "\u00a7", "\u2016")
    ),
    # Strip smoker flag from urine cells — ‖ moved to column header
    `Urine (ppb)` = str_replace_all(`Urine (ppb)`, "\u00a7", "")
  ) |>
  select(
    `Usage Class`,
    Name,
    CAS,
    `IARC Group` = IARC_Group,
    `Mean Non-Cancer Thyroid Conc. (ppb)` = mean_ctrl,
    `Mean Tumor Conc. (ppb)` = mean_tumor,
    `Range (ppb)†` = range,
    `Adipose Tissue (ppb)‡` = `Adipose Tissue (ppb)`,
    `Urine (ppb)¶‖` = `Urine (ppb)`,
    `Serum/Plasma (ppb)¶` = `Serum/Plasma (ppb)`
  ) |>
  arrange(`Usage Class`, Name)
#_ Build hierarchical structure with Usage Class headers
ST3_list <- list()
all_usage_classes <- sort(unique(ST3_tibble_flat$`Usage Class`))
for (usage_idx in seq_along(all_usage_classes)) {
  usage_class <- all_usage_classes[usage_idx]
  ST3_list[[length(ST3_list) + 1]] <- tibble(
    Name = usage_class,
    CAS = NA_character_,
    `IARC Group` = NA_character_,
    `Mean Non-Cancer Thyroid Conc. (ppb)` = NA_character_,
    `Mean Tumor Conc. (ppb)` = NA_character_,
    `Range (ppb)†` = NA_character_,
    `Adipose Tissue (ppb)‡` = NA_character_,
    `Urine (ppb)¶‖` = NA_character_,
    `Serum/Plasma (ppb)¶` = NA_character_
  )
  usage_data <- ST3_tibble_flat |>
    filter(`Usage Class` == usage_class) |>
    mutate(Name = paste0("  ", Name)) |>
    select(-`Usage Class`)
  ST3_list[[length(ST3_list) + 1]] <- usage_data
  if (usage_idx < length(all_usage_classes)) {
    ST3_list[[length(ST3_list) + 1]] <- tibble(
      Name = "", CAS = "", `IARC Group` = "",
      `Mean Non-Cancer Thyroid Conc. (ppb)` = "", `Mean Tumor Conc. (ppb)` = "",
      `Range (ppb)†` = "", `Adipose Tissue (ppb)‡` = "",
      `Urine (ppb)¶‖` = "", `Serum/Plasma (ppb)¶` = ""
    )
  }
}
#_ Combine and fill missing literature values with en-dash
ST3_tibble <- bind_rows(ST3_list) |>
  mutate(
    across(
      c(`Adipose Tissue (ppb)‡`, `Urine (ppb)¶‖`, `Serum/Plasma (ppb)¶`),
      ~ case_when(
        Name != "" & !is.na(CAS) & (is.na(.) | . == "") ~ "\u2013",
        TRUE ~ .
      )
    )
  )
#_ Export to Excel
n_cols <- ncol(ST3_tibble)
wb <- createWorkbook()
addWorksheet(wb, "Table 4")
writeData(wb, sheet = 1, x = ST3_tibble, startRow = 1, colNames = TRUE)
header_style_left <- createStyle(
  fontSize = 7.5, fontName = "Times New Roman", textDecoration = "bold",
  halign = "left", valign = "center", fgFill = "#D9D9D9"
)
header_style_center <- createStyle(
  fontSize = 7.5, fontName = "Times New Roman", textDecoration = "bold",
  halign = "center", valign = "center", fgFill = "#D9D9D9"
)
data_style_left <- createStyle(
  fontSize = 7, fontName = "Times New Roman", halign = "left", valign = "bottom"
)
data_style_center <- createStyle(
  fontSize = 7, fontName = "Times New Roman", halign = "center", valign = "bottom"
)
group_style <- createStyle(
  fontSize = 7, fontName = "Times New Roman", textDecoration = "bold",
  halign = "left", valign = "bottom"
)
addStyle(wb, 1, header_style_left,   rows = 1, cols = 1,          gridExpand = TRUE, stack = FALSE)
addStyle(wb, 1, header_style_center, rows = 1, cols = 2:n_cols,   gridExpand = TRUE, stack = FALSE)
data_rows <- 2:(nrow(ST3_tibble) + 1)
addStyle(wb, 1, data_style_left,   rows = data_rows, cols = 1,        gridExpand = TRUE, stack = FALSE)
addStyle(wb, 1, data_style_center, rows = data_rows, cols = 2:n_cols, gridExpand = TRUE, stack = FALSE)
# Bold Usage Class header rows
group_rows <- which(!is.na(ST3_tibble$CAS) == FALSE & ST3_tibble$Name != "") + 1
addStyle(wb, 1, group_style, rows = group_rows, cols = 1, gridExpand = TRUE, stack = FALSE)
# Header borders
addStyle(wb, 1, createStyle(border = "top",    borderColour = "black", borderStyle = "thin"),
         rows = 1, cols = 1:n_cols, gridExpand = TRUE, stack = TRUE)
addStyle(wb, 1, createStyle(border = "bottom", borderColour = "black", borderStyle = "thin"),
         rows = 1, cols = 1:n_cols, gridExpand = TRUE, stack = TRUE)
# Bottom double border
bottom_row <- nrow(ST3_tibble) + 2
writeData(wb, 1, x = "", startRow = bottom_row, startCol = 1)
addStyle(wb, 1, createStyle(border = "bottom", borderColour = "black", borderStyle = "double"),
         rows = bottom_row, cols = 1:n_cols, gridExpand = TRUE, stack = TRUE)
# Footnote
footnote_row <- bottom_row + 1
mergeCells(wb, 1, cols = 1:n_cols, rows = footnote_row)
footnote_text <- paste(
  "\u2020 The lower bound of the range is the \u2018theoretical minimum,\u2019 which is one half the lowest detected value. This value was used for imputing missing values to determine means.",
  "\u2021 PAH measurements were obtained from a study that analyzed multiple adipose tissue depots across two geographically distinct cohorts\u1d43, while PCB measurements were derived from a separate study that also included two geographically distinct cohorts\u1d47. Thus, multiple means were listed for various subgroups within these studies. However, for each compound, the highest mean value reported across any cohort or cohort-tissue depot combination is presented.",
  "\u00b6 Some urine or plasma values are derived from the CDC\u2019s Biomonitoring Data Tables for Environmental Chemicals\u1d48, which is based on data from NHANES cohorts. In the cases where compounds had data listed from multiple cohorts, the maximum geometric mean across any given cohort is listed.",
  "\u2016 All reported urine values were derived from measurements listed separately for smokers and non-smokers; values from smokers are presented, as these were the highest reported.",
  "References: \u1d43 = (Ml\u0079czy\u0144ska et al., 2023); \u1d47 = (Wang et al., 2010); \u1d9c = (Riffelmann et al., 1995); \u1d48 = (CDC, 2024); \u1d49 = (Maier, 2022)",
  "Abbreviations: CDC = Centers for Disease Control and Prevention; IARC = International Agency for Research on Cancer; NHANES = National Health and Nutrition Examination Survey; PAH = polycyclic aromatic hydrocarbon; PCB = polychlorinated biphenyl; ppb = parts per billion",
  sep = "\n"
)
writeData(wb, 1, x = footnote_text, startRow = footnote_row, startCol = 1)
addStyle(wb, 1,
  createStyle(fontSize = 8, fontName = "Times New Roman", textDecoration = "italic",
              halign = "left", valign = "center", wrapText = TRUE,
              border = "bottom", borderColour = "black", borderStyle = "double"),
  rows = footnote_row, cols = 1:n_cols, gridExpand = TRUE, stack = TRUE)
setRowHeights(wb, 1, rows = footnote_row, heights = 100)
# Column widths
setColWidths(wb, 1, cols = 1,         widths = 30)  # Name
setColWidths(wb, 1, cols = 2,         widths = 12)  # CAS
setColWidths(wb, 1, cols = 3,         widths = 10)  # IARC Group
setColWidths(wb, 1, cols = 4:5,       widths = 20)  # Mean concentrations
setColWidths(wb, 1, cols = 6,         widths = 18)  # Range
setColWidths(wb, 1, cols = 7:9,       widths = 16)  # Literature values
saveWorkbook(wb, export_path, overwrite = TRUE)
cat("Table 4 exported:", export_path, "\n")
ST3_tibble
}
