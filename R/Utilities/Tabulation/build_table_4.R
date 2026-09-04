#' Build Table 4 (Observed Concentrations vs. Literature Values)
#'
#' Constructs and exports the supplementary table comparing tumor/non-cancer
#' thyroid concentrations to literature values for select IARC-classified chemicals.
#' Encapsulates all data wrangling (steps 16.3.0-16.3.9) and Excel export.
#'
#' @param ppm_full_table ppm data table from the pipeline
#' @param ST1_tibble ST1 tibble (used to join chemical names)
#' @param literature_ST3 Usage-class lookup (from primary_data.xlsx, sheet "literature_comp_pared")
#' @param literature_long Audited long-form literature values (sheet "literature_long"); supplies the
#'   three literature columns, one row per (reference, compound, matrix), already in ppb
#' @param validation_check_files_unfiltered Unfiltered validation check files
#' @param export_path Path to export the Excel file
#'
#' @return The final ST3 tibble (hierarchical, formatted)
#'
build_table_4 <- function(ppm_full_table, ST1_tibble, literature_ST3, literature_long,
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
#! trimws() on every format() call: format() pads a vector to a common width, so a 4-digit value like 1,248 silently left-pads every shorter cell in the same column. The intentional spaces in "< 1" and in the "a - b" range are inside the strings and unaffected.
  mutate(
    across(c(mean_ctrl, mean_tumor),
           ~ if_else(. < 1, "< 1", trimws(format(round(.), big.mark = ","))),
           .names = "{.col}"),
    half_min = if_else(half_min < 1, "< 1", trimws(format(round(half_min), big.mark = ","))),
    max_value = if_else(max_value < 1, "< 1", trimws(format(round(max_value), big.mark = ",")))
  ) |>
  left_join(ST1_tibble |> select(Name, CAS), by = "CAS") |>
#! Unspaced en dash: number ranges close up (< 1–81), and this matches the pasted Table 4 so a re-paste cannot silently revert it.
  mutate(range = paste0(trimws(half_min), "–", trimws(max_value))) |>
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
  "207-08-9", # Benzo[k]fluoranthene
  "92-87-5"  # Benzidine
)
#_ Build flat tibble with literature values
# citation letters (reading order): ᵃ=mlyczynska2023, ᵇ=wang2010, ᶜ=riffelmann1995, ᵈ=cdc2024, ᵉ=maier2022
#! Keyed on BOTH the reference slug and the single letter: literature_long carries the letter directly, while the older literature_comp_pared columns carry the slug. Same superscript either way.
ref_map <- c(
  mlyczynska2023 = "\u1d43", a = "\u1d43",  # ᵃ
  maier2022      = "\u1d49", e = "\u1d49",  # ᵉ
  riffelmann1995 = "\u1d9c", c = "\u1d9c",  # ᶜ
  cdc2024        = "\u1d48", d = "\u1d48",  # ᵈ
  wang2010       = "\u1d47", b = "\u1d47"   # ᵇ
)
resolve_ref <- function(ref) {
  dplyr::coalesce(ref_map[ref], ref)
}
#_ Literature columns, derived from the audited long-form sheet
#! literature_long holds one row per (reference, compound, matrix), each already the highest subgroup mean within its own source. Here we take the max ACROSS sources per (CAS, column) and carry that winner's citation letter, so the superscript always names the study the printed number came from. Hb adducts are excluded: they measure cumulative covalent binding to erythrocyte protein over the ~120-day red-cell lifespan, not a circulating concentration, so they share no axis with the ng/g and ng/mL values here. Riffelmann reports 4-aminobiphenyl ONLY as an adduct, so its serum/plasma cell is correctly empty; its urine comparator comes from CDC. Plasma and Serum both feed the single "Serum/Plasma" column; AT and Urine map 1:1.
lit_col <- c(AT = "Adipose Tissue (ppb)", Urine = "Urine (ppb)",
             Plasma = "Serum/Plasma (ppb)", Serum = "Serum/Plasma (ppb)")
literature_wide <- literature_long |>
  filter(matrix != "Hb (adduct)", !is.na(value_ppb)) |>
  mutate(column = unname(lit_col[matrix])) |>
#! Ties get BOTH citations rather than an arbitrary winner. Fluoranthene plasma is the only case: Mlyczynska 0.90 and Maier 0.9 are numerically equal, and both studies genuinely report that value, so citing one would imply the other did not. Letters are sorted so the pair is deterministic.
  group_by(CAS, column) |>
  filter(value_ppb == max(value_ppb)) |>
  summarise(value_ppb = first(value_ppb),
#! Comma-separate tied citation letters: bare "ᵃᵉ" reads as one glyph pair, "ᵃ,ᵉ" reads as two refs.
            refs = paste0(resolve_ref(sort(unique(letter))), collapse = ","),
            .groups = "drop") |>
#! Same rounding rule as the study columns: "< 1" below unity, nearest integer otherwise. trimws() because format() pads a vector to a common width, which would leave stray leading spaces.
  mutate(cell = paste0(if_else(value_ppb < 1, "< 1",
                               trimws(format(round(value_ppb), big.mark = ","))),
                       refs)) |>
  select(CAS, column, cell) |>
  pivot_wider(names_from = column, values_from = cell)
ST3_tibble_flat <- ST3_tibble_2 |>
  filter(!CAS %in% failed_ST3) |>
  left_join(literature_ST3 |> select(CAS, `Usage Class`), by = "CAS") |>
  left_join(literature_wide, by = "CAS") |>
  mutate(
    Name = gsub("[†‡]", "", Name),
    Name = gsub("\\*", "", Name)
  ) |>
  select(
    `Usage Class`,
    Name,
    CAS,
    `IARC Group` = IARC_Group,
    `Mean Non-Cancer Thyroid Conc. (ppb)` = mean_ctrl,
    `Mean Tumor Conc. (ppb)` = mean_tumor,
    `Study Range (ppb)†` = range,
    `Adipose Tissue (ppb)‡` = `Adipose Tissue (ppb)`,
    `Urine (ppb)‡` = `Urine (ppb)`,
    `Serum/Plasma (ppb)‡` = `Serum/Plasma (ppb)`
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
    `Study Range (ppb)†` = NA_character_,
    `Adipose Tissue (ppb)‡` = NA_character_,
    `Urine (ppb)‡` = NA_character_,
    `Serum/Plasma (ppb)‡` = NA_character_
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
      `Study Range (ppb)†` = "", `Adipose Tissue (ppb)‡` = "",
      `Urine (ppb)‡` = "", `Serum/Plasma (ppb)‡` = ""
    )
  }
}
#_ Combine and fill missing literature values with en-dash
ST3_tibble <- bind_rows(ST3_list) |>
  mutate(
    across(
      c(`Adipose Tissue (ppb)‡`, `Urine (ppb)‡`, `Serum/Plasma (ppb)‡`),
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
  "\u2021 In several of the studies cited, values were reported for multiple subgroups within a given biological matrix. This included non-smoking vs. smoking, occupationally exposed vs. non-exposed, geographically distinct cohorts, varying adipose tissue depots, and different NHANES survey years. In all these cases, the highest reported mean (or geometric mean for NHANES) for any subgroup is summarized in this table. When more than one study reported values for the same compound in the same matrix, the highest value is also listed in this table among those studies. All NHANES values represent total population (values for demographic subgroups not used).",
  "References: \u1d43 = (Ml\u0079czy\u0144ska et al., 2023); \u1d47 = (Wang et al., 2010); \u1d9c = (Riffelmann et al., 1995); \u1d48 = (\u201cBiomonitoring Data Tables for Environmental Chemicals | CDC,\u201d 2024); \u1d49 = (Vermillion Maier et al., 2022)",
  "Abbreviations: CDC = Centers for Disease Control and Prevention; IARC = International Agency for Research on Cancer; NHANES = National Health and Nutrition Examination Survey; PAH = polycyclic aromatic hydrocarbon; PCB = polychlorinated biphenyl; ppb = parts per billion",
  sep = "\n"
)
writeData(wb, 1, x = footnote_text, startRow = footnote_row, startCol = 1)
addStyle(wb, 1,
  createStyle(fontSize = 8, fontName = "Times New Roman",
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
