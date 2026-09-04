#* 0c: Feature Table Import and Preprocessing
#+ 0c.1: Import data, clean, Preprocess
#- 0c.1.1: Import raw feature table data
tumor_raw <- read_excel(config$paths$primary_data, sheet = "lib.subject.summary")
#! ============================================================================
#! TERMINOLOGY: "variant" (code) == "type" (manuscript / pathology consensus)
#! ----------------------------------------------------------------------------
#! Throughout this pipeline, the three differentiated thyroid cancer types
#! analyzed -- papillary (PTC), follicular (FTC), and the invasive encapsulated
#! follicular variant of papillary thyroid carcinoma (IEFVPTC) -- are named
#! "variant" in column names, objects, functions, scripts, and file paths. This
#! is an intentional, isolated departure from the manuscript's terminology.
#! Per the 2022 WHO Classification of Thyroid Tumours (Jung et al., 2022; WHO
#! Editorial Board, 2022), "variant" is reserved for genetic variants and
#! IEFVPTC is now a distinct entity rather than a subtype of PTC; the manuscript
#! therefore uses "type"/"tumor type" throughout. We deliberately retain
#! "variant" as an internal identifier solely to preserve the integrity and
#! reproducibility of a validated, working pipeline -- renaming risks
#! introducing errors into otherwise verified analysis code. In every case,
#! read "variant" in the code as "type" in the manuscript; the discrepancy is
#! purely nominal and confined to code-level naming.
#! ============================================================================
#- 0c.1.2: Import tumor sequence/variant data
tumor_seq <- read_excel(config$paths$primary_data, sheet = "tumors_sequence") |>
  select(ID, variant) |>
  unique()
#- 0c.1.3: Import feature metadata
#! Source data has "Dye intermediates" with a lowercase i while every sibling class is Title Case; normalize at import so Tables 2/3/4, ST2 and the figures all agree.
feature_metadata <- read_excel(config$paths$chemical_metadata, sheet = "feature_metadata") |>
  mutate(across(any_of(c("Table_Class", "Graph_Class")),
                \(x) str_replace(x, "^Dye intermediates$", "Dye Intermediates"))) |>
#! Graph_Class carries "Surfactants or Detergent" (singular) while Table_Class, the manuscript prose,
#! Table 2 and ST2 all carry "Detergents". Only the Figure 2A bar label was wrong. Anchored, so it
#! cannot touch the already-correct Table_Class value.
  mutate(across(any_of(c("Table_Class", "Graph_Class")),
                \(x) str_replace(x, "^Surfactants or Detergent$", "Surfactants or Detergents"))) |>
#! 07_validation_prep.R reads Short_display_name off this object rather than off the short_name tibble built in 04, so the shared nomenclature fixes have to land here too or the validation figure titles keep the source spelling.
  mutate(across(any_of(c("name", "Short_display_name")), fix_chem_nomenclature))
#- 0c.1.4: Import and clean tissue weights
weights <- read_excel(config$paths$primary_data, sheet = "tissue_weights") |>
  filter(samples == "Tumor") |>
  select(ID, weight_mg) |>
  rename(patient_ID = ID)
#- 0c.1.5: Import absolute quant data
conc_raw <- read_excel(config$paths$primary_data, sheet = "lib.subject.qsummary", col_type = "text") |>
  select(name_sub_lib_id, F1:F20) |>
  mutate(across(-name_sub_lib_id, as.numeric)) |>
  mutate(across(-name_sub_lib_id, ~ .x * (0.47 / 0.5)))
# Apply correction factor: original pipeline used 0.5 ng/mL, actual is 0.47 ng/mL
#- 0c.1.6: Import cadaver control tissue weights; clean
cadaver_tissue_wts <- read_excel(config$paths$primary_data, sheet = "tissue_weights") |>
  filter(samples == "Control") |>
  select(ID, weight_mg) |>
  rename(control_ID = ID)
#- 0c.1.7: Import cadaver absolute quant, clean
cadaver_qraw_i <- read_excel(config$paths$primary_data, sheet = "lib.subject.qsummary.cadaver") |>
  select(name_sub_lib_id, T001:T009) |>
  mutate(across(-name_sub_lib_id, as.numeric)) |>
  mutate(across(-name_sub_lib_id, ~ .x * (0.47 / 0.5))) |>
  pivot_longer(cols = -name_sub_lib_id, names_to = "control_ID", values_to = "Ce") |>
  mutate(C = as.numeric(gsub(",", "", Ce)))
# Apply correction factor: original pipeline used 0.5 ng/mL, actual is 0.47 ng/mL
#- 0c.1.8: Import in fragment quality info, clean
fragment_quality_info <- read_excel(config$paths$primary_data, sheet = "lib.subject.qsummary") |>
  arrange(cas) |>
  select(name_sub_lib_id, iMean) |>
  rename(iMean_tumors = iMean)
#- 0c.1.9: Import clean absolute quant for control IARCs
IARC_controls_ii <- read_excel(config$paths$primary_data, sheet = "lib.subject.qsummary.cadaver")
#- 0c.1.10: Import clean absolute quant for tumor IARCs
IARC_tumors_ii <- read_excel(config$paths$primary_data, sheet = "lib.subject.qsummary")
#- 0c.1.11: Import and clean library
ST1_import <- read_excel(config$paths$primary_data, sheet = "library") |>
  filter(Disposition != "Endogenous") |>
  mutate(subid_col = paste0("mz", subid)) |>
  select(id, name = st1_name, trt, monoisotopic, cas, formula, Disposition, subid_col, tmz) |>
  #! Display-casing corrections pending a source fix in the library sheet; root name is Title-Case in tables per convention
  mutate(
    name = str_replace(name, "^propoxur", "Propoxur"),
    name = str_replace(name, "benzenehexachloride", "Benzenehexachloride"),
    name = str_replace(name, "HXCDD", "HxCDD"),
    name = fix_chem_nomenclature(name)
  ) |>
  mutate(
    name = if_else(
      Disposition == "Exogenous and Endogenous" & !is.na(name),
      paste0(name, "*"),
      name
    )
  ) |>
  distinct() |>
  pivot_wider(
    names_from = subid_col,
    values_from = tmz
  ) |>
  arrange(cas)
#- 0c.1.12: ST1 abbreviations
#! Lives in-repo (not primary_data.xlsx) because it is display metadata, not measurement data; version-controlled so dictionary edits are diffable. Use == "Y" selects the subset shown in the supplement
ST_abbrevs <- read_tsv(here::here("Supplementary", "Components", "abbreviations.tsv"), show_col_types = FALSE, na = "") |>
  mutate(formatted = paste(Abbrev, "=", Name)) |>
  filter(Use == "Y") |>
  arrange(tolower(Abbrev)) |>
  select(formatted)
#- 0c.1.13: Import file list
file_list <- read_excel(config$paths$primary_data, sheet = "file_list") |>
  select(file, ID, replicate, study, type) |>
  group_by(ID, study, type) |>
  summarize(files = paste(file, collapse = ", "), .groups = "drop")
#- 0c.1.14: Import GC2 feature list
GC2_features <- read_csv(config$paths$gc2_features)
#- 0c.1.15: Import expanded library
expanded_lib <- read_csv(config$paths$gc2_expanded)
#- 0c.1.16: Validation file (Unfiltered)
validation_check_files_unfiltered <- read_xlsx(config$paths$validation, sheet = "validation")
#- 0c.1.17: Validation file (Filtered)
validation_check_files <- validation_check_files_unfiltered |>
  filter(!state %in% c("failed", "alternate")) |>
  mutate(rt_range = (rtu-rtl)/2) |>
  select(state, everything(), -c(modification, note, rtl, rtu))
#- 0c.1.18: Validation plot metadata
validation_plot_metadata_ordered <- read_xlsx(config$paths$validation, sheet = "figure_order") %>%
  arrange(order) %>%
  mutate(
    sf_sub = paste(figure, subfigure, sep = "."),
    full_path = here::here(plot),
    grob = map(full_path, readRDS)
  ) %>%
  select(order, id, short_name, figure, subfigure, sf_sub, panel, plot, full_path, grob)
#- 0c.1.19: ST3 literature review
literature_ST3 <- read_excel(config$paths$primary_data, sheet = "literature_comp_pared") |>
  select(CAS, `Usage Class`, AT_manuscript, AT_ref, urine_manuscript, urine_ref, plasma_manuscript, plasma_ref)
#- 0c.1.20: Literature review, long form -- one row per (reference, compound, matrix)
#! Audited row by row against the source PDFs; supersedes literature_comp_pared as the source of the three literature columns in Table 4. "-" marks "not reported by that source" and becomes NA here. ppb is ng/g for solids and ng/mL for fluids; ug/L equals ng/mL at density ~1, so those are already ppb, while pg/mL and ng/L are 1000x smaller.
LIT_TO_PPB <- c("ng/g" = 1, "ng/mL" = 1, "ug/L" = 1, "pg/mL" = 1e-3, "ng/L" = 1e-3)
literature_long <- read_excel(config$paths$primary_data, sheet = "literature_long") |>
  mutate(across(where(is.character), \(x) na_if(x, "-")),
         value = as.numeric(value),
         value_ppb = value * unname(LIT_TO_PPB[unit]))
#! Fails loudly on an unrecognised unit rather than silently dropping that row from the maxima.
stopifnot("literature_long: unit with no ppb conversion" =
            all(is.na(literature_long$value) | !is.na(literature_long$value_ppb)))
#- 0c.1.21: Defective library entries excluded from ALL analyses
#! CP2302 ("tris(tribromoneopentyl)", CAS 21850-44-2, C21H20Br8O2) is removed because its recorded target ions cannot belong to the annotated compound. The four targets (647.4252 / 648.4294 / 649.4321 / 664.4547) are one species: the first three are a 13C ladder (spacings 1.0042, 1.0027; the +2 ion is 2x13C at 0.3 ppm, NOT 81Br, which misses by 13.8 ppm -- outside the file's own +/- 6 ppm window), and 664.4547 is the [M+NH4]+ partner of 647.4252 as [M+H]+ (delta 17.0295 = NH3). No 79/81Br isotope envelope is present anywhere, yet 647.4252 exceeds the heaviest halogen-free subformula of C21H20Br8O2 (304.15), so the ion MUST contain Br to be a fragment of the parent. An octabrominated compound also has its base isotope peak at M+8 with 2 Da spacing, not the observed 1 Da ladder. The claim is that the recorded ions are inconsistent with the annotated compound -- NOT that the compound is absent. Manual spectral validation cannot catch this: it confirms reproducibility, not identity. Exclusion lives here (in code, at the single external-data gateway) rather than in the OneDrive spreadsheets so it is version-controlled, diffable, and reversible in one line.
EXCLUDED_LIB_IDS <- c("CP2302")
EXCLUDED_LIB_CAS <- c("21850-44-2")
#' Drop excluded library entries from any per-chemical frame
#' @details Matches on whichever key the frame carries: `id`, `cas`/`CAS` (case-insensitive),
#'   or `name_sub_lib_id` (which embeds the library id, so all subids are caught).
drop_excluded <- function(df) {
  nm <- names(df)
  if ("id" %in% nm) df <- filter(df, !id %in% EXCLUDED_LIB_IDS)
  cas_col <- nm[tolower(nm) == "cas"]
  if (length(cas_col)) df <- filter(df, !as.character(.data[[cas_col[1]]]) %in% EXCLUDED_LIB_CAS)
  if ("name_sub_lib_id" %in% nm)
    df <- filter(df, !str_detect(name_sub_lib_id, paste(EXCLUDED_LIB_IDS, collapse = "|")))
  df
}
#! Applied to every imported frame that carries per-chemical rows. The loop reports each drop so a silent no-op (e.g. a renamed key column) is visible in the run log.
cat("\n-- Excluded library entries (", paste(EXCLUDED_LIB_IDS, collapse = ", "), ") --\n", sep = "")
for (.obj in c("tumor_raw", "feature_metadata", "ST1_import", "conc_raw", "cadaver_qraw_i",
               "fragment_quality_info", "IARC_controls_ii", "IARC_tumors_ii", "GC2_features",
               "expanded_lib", "validation_check_files_unfiltered", "validation_check_files",
               "validation_plot_metadata_ordered", "literature_ST3", "literature_long")) {
  if (!exists(.obj)) { cat(sprintf("  %-34s MISSING (skipped)\n", .obj)); next }
  .before <- nrow(get(.obj)); assign(.obj, drop_excluded(get(.obj))); .after <- nrow(get(.obj))
  cat(sprintf("  %-34s %5d -> %5d  (%d dropped)\n", .obj, .before, .after, .before - .after))
}
rm(.obj, .before, .after)
#- 0c.1.22: Drop the BLANK spacer orphaned by the exclusion
#! The figure_order sheet is HAND-CURATED, not a simple alternating list: a chemical with two plots occupies top AND bottom of the same page (main plot above its isolated-fragment plot), and single-plot chemicals are padded with explicit BLANK spacer rows. Do NOT re-derive `panel` by row position -- that reorders curated pairs and would place isolated-fragment panels above their parents.
#! CP2302 sat at order 43 (panel=top) with a BLANK at order 43.5 (panel=bottom), i.e. it held a page of its own. Removing the compound orphans that spacer, so it goes too. Net effect: S2.1 loses one whole page (22 -> 21) and downstream page numbers shift by 1.
validation_plot_metadata_ordered <- validation_plot_metadata_ordered |>
  filter(!(sf_sub == "2.1" & short_name == "BLANK" & order == 43.5))
#! Every page holds exactly two rows, so an odd count in any sf_sub means a broken pairing. Fail loudly here rather than silently emitting a half-blank page into the supplement.
.pp <- validation_plot_metadata_ordered |> count(sf_sub, name = "rows") |> mutate(pages = rows / 2)
cat("-- Validation figure layout after exclusion --\n")
print(as.data.frame(.pp), row.names = FALSE)
if (any(.pp$rows %% 2 != 0))
  stop("figure_order: odd row count in sf_sub ",
       paste(.pp$sf_sub[.pp$rows %% 2 != 0], collapse = ", "),
       " -- top/bottom pairing is broken; fix before building the supplement.")
rm(.pp)
#+ 0c.2: Structure Data
#- 0c.2.1: Pull the tumor columns
tumor_column <- tumor_raw |>
  select(name_sub_lib_id, "F1":"F20")
#- 0c.2.2: Check for all-zero rows in tumor_raw
# Get column names matching the pattern (F, P, or FVPTC followed by numbers)
variant_cols <- colnames(tumor_raw)[grepl("^(F|P|FVPTC)\\d+$", colnames(tumor_raw))]
# Check each row has at least one non-zero value in variant columns
has_value <- tumor_raw %>%
  select(all_of(variant_cols)) %>%
  mutate(across(everything(), ~ replace_na(., 0))) %>%
  rowwise() %>%
  mutate(has_nonzero = any(c_across(everything()) != 0)) %>%
  ungroup() %>%
  pull(has_nonzero)
# Summary
cat("Total rows:", nrow(tumor_raw), "\n")
cat("Rows with at least one non-zero value:", sum(has_value), "\n")
cat("Rows with all zeros:", sum(!has_value), "\n")
#- 0c.2.3: Pivot longer/wider to get into analysis format
tumor <- tumor_column |>
  pivot_longer(-name_sub_lib_id, names_to = "ID", values_to = "Value") |>
  pivot_wider(names_from = name_sub_lib_id, values_from = Value) |>
  left_join(tumor_seq, by = "ID") |>
  select(ID, variant, everything(), -ID)
#+ 0c.3: Split Into Quantitative and Qualitative Features Based on Missingness
#- 0c.3.1: Calculate proportion of missing values per feature
{
  na_threshold <- 0.3
  PMV <- colMeans(is.na(tumor))
  cols_quant <- names(PMV[PMV <= na_threshold])
  cols_qual <- names(PMV[PMV > na_threshold])
  cols_quant <- unique(c("variant", cols_quant))
  cols_qual <- setdiff(cols_qual, "variant")
}
#- 0c.3.2: Create quantitative table with 1/2 minimum imputation and log2 transform
tumors_quant <- tumor |>
  select(all_of(cols_quant)) |>
  mutate(across(-variant, ~ ifelse(is.na(.), 0.5 * min(., na.rm = TRUE), .))) |>
  mutate(across(-variant, ~ log2(.)))
#- 0c.3.3: Create qualitative table with binary presence/absence
tumors_qual <- tumor |>
  select(all_of(c("variant", cols_qual))) |>
  mutate(across(-variant, ~ ifelse(is.na(.), 0, 1))) |>
  select(where(~ any(. != 0) | is.character(.)))
#- 0c.3.4: Qualitative table with patient_ID retained, for script 13's covariate join
{
  tumor_id <- tumor_column |>
    pivot_longer(-name_sub_lib_id, names_to = "ID", values_to = "Value") |>
    pivot_wider(names_from = name_sub_lib_id, values_from = Value) |>
    left_join(tumor_seq, by = "ID") |>
    rename(patient_ID = ID)
  tumors_qual_id <- tumor_id |>
    select(all_of(c("patient_ID", "variant", cols_qual))) |>
    mutate(across(-c(variant, patient_ID), ~ ifelse(is.na(.), 0, 1))) |>
    select(where(~ any(. != 0) | is.character(.)))
}
#+ 0c.4: Clean Up Targeted Feature Tables
#- 0c.4.1: Normalize by weights
tumors_quant_wt_i <- tumor_column |>
  pivot_longer(-name_sub_lib_id, names_to = "ID", values_to = "Value") |>
  pivot_wider(names_from = name_sub_lib_id, values_from = Value) |>
  mutate(patient_ID = ID) |>
  left_join(tumor_seq, by = "ID") |>
  select(patient_ID, all_of(cols_quant)) |>
  mutate(across(-c(variant, patient_ID), ~ as.numeric(.))) |>
  mutate(across(-c(variant, patient_ID), ~ ifelse(is.na(.), 0.5 * min(., na.rm = TRUE), .))) |>
  left_join(weights, by = "patient_ID") |> # join tissue weights in
  mutate(across(where(is.numeric) & !matches("weight_mg"), ~ . / weight_mg)) |>
  select(-weight_mg) |>
  mutate(across(-c(variant, patient_ID), ~ log2(.)))
#- 0c.4.2: Pull endogenous features
endog_cas <- read_excel(config$paths$chemical_metadata, sheet = "endogenous_excluded_features") |>
  pull(cas)
#- 0c.4.3: Get endogenous feature key
cas_key_endog <- tumor_raw |>
  select(name_sub_lib_id, cas) |>
  filter(cas %in% endog_cas) |>
  pull(name_sub_lib_id)
#- 0c.4.4: Remove endogenous features from quantitative weighted table
tumors_quant_wt_id <- tumors_quant_wt_i |>
  select(-any_of(cas_key_endog))
#- 0c.4.5: Strip patient ID for downstream analysis
tumors_quant_wt <- tumors_quant_wt_id |>
  select(-patient_ID)
