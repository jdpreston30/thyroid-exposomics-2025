#* 17: Supplementary Tables
#+ 17.1: ST1: Chemical Library (pivoted subid)
#- 17.1.1: Build and format ST1 gt table
gt_ST1 <- build_ST1(ST1_tibble)
#- 17.1.2: Save ST1 as LaTeX (without table wrapper) to Supplementary/Components/Tables
latex_code <- gt::as_latex(gt_ST1) |> as.character()
latex_code <- fix_latex_header_fill(latex_code)
# Remove table wrapper for direct inclusion in supplementary
latex_lines <- strsplit(latex_code, "\n")[[1]]
latex_lines <- latex_lines[-c(1, length(latex_lines))]  # Remove \begin{table} and \end{table}
latex_code <- paste(latex_lines, collapse = "\n")
writeLines(latex_code, "Supplementary/Components/Tables/ST1.tex")
#+ 17.2: ST2
#- 17.2.1: Prepare base ST2 data
ST2_base <- feature_metadata |>
  mutate(`Superclass: Class` = if_else(is.na(Superclass), Class, paste(Superclass, Class, sep = ": "))) |>
  select(GROUP = Table_Header, Class = Table_Class, Subclass = Table_Subclass, CAS = cas, `Potential EDC` = Potential_EDC, `IARC Group` = IARC_Group, `Superclass: Class`) |>
  mutate(GROUP = toupper(GROUP)) |>
  left_join(ST1_tibble |> select(Name, CAS), by = "CAS") |>
  mutate(Name = gsub("\u2021", "", Name)) |>
  mutate(Name = gsub("\\*", "\u2020", Name)) |>
  arrange(GROUP, Class, Subclass, Name)
#- 17.2.2: Build hierarchical structure with proper nesting
ST2_list <- list()
all_groups <- unique(ST2_base$GROUP)
for (group_idx in seq_along(all_groups)) {
  group <- all_groups[group_idx]
  # Add GROUP header
  ST2_list[[length(ST2_list) + 1]] <- tibble(
    Display_Name = group,
    row_type = "GROUP",
    CAS = NA_character_,
    `Potential EDC` = NA_character_,
    `IARC Group` = NA_character_,
    `Superclass: Class` = NA_character_
  )
  # Get classes for this group
  group_data <- ST2_base |> filter(GROUP == group)
  all_classes <- unique(group_data$Class)
  for (class_idx in seq_along(all_classes)) {
    class <- all_classes[class_idx]
    # Add Class header with 2 spaces
    ST2_list[[length(ST2_list) + 1]] <- tibble(
      Display_Name = paste0("  ", class),
      row_type = "Class",
      CAS = NA_character_,
      `Potential EDC` = NA_character_,
      `IARC Group` = NA_character_,
      `Superclass: Class` = NA_character_
    )
    # Get subclasses for this class (sort with "Other" last)
    class_data <- group_data |> filter(Class == class)
    subclasses <- unique(class_data$Subclass)
    # Remove NA/empty subclasses - these will be listed directly under Class
    subclasses <- subclasses[!is.na(subclasses) & subclasses != ""]
    
    if (length(subclasses) == 0) {
      # No meaningful subclasses - add chemicals directly under Class with 6 spaces
      class_chemicals <- class_data |>
        mutate(
          Display_Name = paste0("      ", Name),
          row_type = "Chemical"
        ) |>
        select(Display_Name, row_type, CAS, `Potential EDC`, `IARC Group`, `Superclass: Class`)
      ST2_list[[length(ST2_list) + 1]] <- class_chemicals
    } else {
      # Has meaningful subclasses - use full hierarchy
      # Move "Other" to the end
      other_subclasses <- subclasses[subclasses == "Other"]
      non_other_subclasses <- subclasses[subclasses != "Other"]
      subclasses <- c(sort(non_other_subclasses), other_subclasses)
      
      for (subclass in subclasses) {
        # Add Subclass header with 4 spaces
        ST2_list[[length(ST2_list) + 1]] <- tibble(
          Display_Name = paste0("    ", subclass),
          row_type = "Subclass",
          CAS = NA_character_,
          `Potential EDC` = NA_character_,
          `IARC Group` = NA_character_,
          `Superclass: Class` = NA_character_
        )
        # Add chemicals for this subclass with 6 spaces
        subclass_data <- class_data |> 
          filter(Subclass == subclass) |>
          mutate(
            Display_Name = paste0("      ", Name),
            row_type = "Chemical"
          ) |>
          select(Display_Name, row_type, CAS, `Potential EDC`, `IARC Group`, `Superclass: Class`)
        ST2_list[[length(ST2_list) + 1]] <- subclass_data
      }
    }
    
    # Add one blank row after each Class (except the last Class in a GROUP)
    if (class_idx < length(all_classes)) {
      ST2_list[[length(ST2_list) + 1]] <- tibble(
        Display_Name = "",
        row_type = "Spacer",
        CAS = "",
        `Potential EDC` = "",
        `IARC Group` = "",
        `Superclass: Class` = ""
      )
    }
    # Add extra spacing row after specific classes for page break control
    if (class %in% c("Herbicides", "Organic UV Filters")) {
      ST2_list[[length(ST2_list) + 1]] <- tibble(
        Display_Name = "",
        row_type = "Spacer",
        CAS = "",
        `Potential EDC` = "",
        `IARC Group` = "",
        `Superclass: Class` = ""
      )
    }
  }
  # Add two blank rows after each GROUP section (except the last one)
  if (group_idx < length(all_groups)) {
    for (i in 1:2) {
      ST2_list[[length(ST2_list) + 1]] <- tibble(
        Display_Name = "",
        row_type = "Spacer",
        CAS = "",
        `Potential EDC` = "",
        `IARC Group` = "",
        `Superclass: Class` = ""
      )
    }
    # Add page break marker after GROUP spacing
    ST2_list[[length(ST2_list) + 1]] <- tibble(
      Display_Name = "PAGEBREAK",
      row_type = "PageBreak",
      CAS = "",
      `Potential EDC` = "",
      `IARC Group` = "",
      `Superclass: Class` = ""
    )
  }
}
#- 17.2.3: Combine and format
ST2_tibble <- bind_rows(ST2_list) |>
  mutate(
    `Potential EDC` = case_when(
      `Potential EDC` == "Y" ~ "✓",
      row_type == "Chemical" & (is.na(`Potential EDC`) | `Potential EDC` == "") ~ "–",
      TRUE ~ `Potential EDC`
    ),
    `IARC Group` = case_when(
      row_type == "Chemical" & (is.na(`IARC Group`) | `IARC Group` == "") ~ "–",
      TRUE ~ `IARC Group`
    )
  ) |>
  select(-row_type)  # Remove helper column
#- 17.2.4: Build and format ST2 gt table
gt_ST2 <- build_ST2(ST2_tibble)
#- 17.2.5: Save ST2 as LaTeX (without table wrapper) to Supplementary/Components/Tables
latex_code <- gt::as_latex(gt_ST2) |> as.character()
latex_code <- fix_ST2_latex(latex_code)
# Remove table wrapper for direct inclusion in supplementary
latex_lines <- strsplit(latex_code, "\n")[[1]]
latex_lines <- latex_lines[-c(1, length(latex_lines))]  # Remove \begin{table} and \end{table}
latex_code <- paste(latex_lines, collapse = "\n")
writeLines(latex_code, "Supplementary/Components/Tables/ST2.tex")
#+ 17.3: Abbreviations Dictionary
#- 17.3.1: Build abbreviations list
abbrev_list <- ST_abbrevs |>
  arrange(formatted)
#- 17.3.2: Convert to LaTeX itemize list with reduced spacing
abbrev_latex <- c(
  "\\begin{itemize}",
  "\\setlength{\\itemsep}{2pt}",
  "\\setlength{\\parskip}{0pt}",
  "\\setlength{\\parsep}{0pt}",
  paste0("  \\item ", abbrev_list$formatted),
  "\\end{itemize}"
)
#- 17.3.3: Save abbreviations to Supplementary/Components/Sections
writeLines(abbrev_latex, "Supplementary/Components/Sections/abbreviations.tex")
#+ 17.4: ST3: Cohort Demographics (DTC Tumors vs Cadaver controls)
#- 17.4.1: Assemble combined demographic frame (age, sex, sample collection timing)
#! Collection year is binned (same scheme as Table 1's Sample_Collection_Timing) so it
#! is summarized categorically as n (%) per bin rather than as a spurious numeric median.
.yr_breaks <- seq(2006, 2022, length.out = 5)
.yr_labels <- c("2006-2009", "2010-2013", "2014-2017", "2018-2021")
demo_dtc <- tumor_pathology_raw |>
  filter(str_detect(Patient_ID, "^(F|P|FVPTC)\\d+$")) |>
  transmute(
    Cohort = "Thyroid Tumor",
    Sex = if_else(as.numeric(Sex) == 1, "Female", "Male"),
    Age = as.numeric(Age),
    `Sample Collection Timing` = cut(as.numeric(year), breaks = .yr_breaks, labels = .yr_labels, include.lowest = TRUE)
  )
demo_cad <- cadaver_metadata |>
  transmute(
    Cohort = "Cadaver Thyroid",
    Sex = if_else(sex == "F", "Female", "Male"),
    Age = as.numeric(age),
    `Sample Collection Timing` = cut(as.numeric(collection_year), breaks = .yr_breaks, labels = .yr_labels, include.lowest = TRUE)
  )
demo_S3_data <- bind_rows(demo_dtc, demo_cad) |>
  mutate(
    Cohort = factor(Cohort, levels = c("Thyroid Tumor", "Cadaver Thyroid")),  # column order: tumor (primary) first
    Sex = factor(Sex, levels = c("Female", "Male"))
  )
#- 17.4.2: Summarize with TernTables (returns display tibble; also writes word version)
ST3_tern <- ternG(
  demo_S3_data,
  group_var = "Cohort",
  methods_doc = FALSE,
  print_normality = FALSE,
  show_test = FALSE,
  show_p = TRUE,
  open_doc = FALSE,
  citation = FALSE,
  font_family = "Times New Roman",
  table_font_size = 10.5,
  indent_info_column = TRUE
)  
#- 17.4.3: Build ST3 LaTeX (self-contained tabular) and save to Supplementary/Components/Tables
st3_latex <- build_ST3(ST3_tern)
writeLines(st3_latex, "Supplementary/Components/Tables/ST3.tex")
#+ 17.5: ST4: Candidate Confounder Assessment (reviewer #2 — confounding)
#! Two same-structured test blocks from script 13 (covariate vs. type; covariate vs.
#! covariate) merged under underlined section headers. The collection-year distribution
#! reconciliation is emitted to ST4_caption.tex so its numbers are wired from R (not
#! hand-typed) — build_ST4() styles the tabular to match ST3 (see its roxygen).
#- 17.5.1: Merge the balance + cross-association blocks under two section headers
#! Script 13 supplies numeric P; formatting happens here so ST3, ST4 and ST5 all round identically
.fmt_p4 <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3)))
ST4_data <- tibble(
  Comparison = c(
    "Covariate vs. Tumor Type", balance_by_type$covariate,
    "Covariate vs. Covariate", covariate_cross$pair
  ),
  Test = c("", balance_by_type$test, "", covariate_cross$test),
  P = c("", .fmt_p4(balance_by_type$P), "", .fmt_p4(covariate_cross$P)),
  .section = c(TRUE, rep(FALSE, nrow(balance_by_type)), TRUE, rep(FALSE, nrow(covariate_cross)))
)
#- 17.5.2: Build ST4 LaTeX (self-contained tabular) and save
writeLines(build_ST4(ST4_data), "Supplementary/Components/Tables/ST4.tex")
#- 17.5.3: Wire the caption deterministically from year_by_type (medians/ranges from R)
.type_order <- c("PTC", "FTC", "IEFVPTC")
.yb <- year_by_type[match(.type_order, as.character(year_by_type$Variant)), ]
.med_fmt <- ifelse(.yb$median == floor(.yb$median),
  formatC(.yb$median, format = "d"), formatC(.yb$median, format = "f", digits = 1)
)
.med_phrase <- paste(sprintf("%s (%s)", .med_fmt, .yb$Variant), collapse = ", ")
.range_lo <- min(year_by_type$min)
.range_hi <- max(year_by_type$max)
st4_caption <- paste0(
  "\\textbf{Candidate confounder assessment for the tumor-type chemical comparisons.} ",
  "Each candidate confounder (age, sex, and collection year) was tested against tumor type, and ",
  "the confounders were tested against one another; the test used for each comparison is given in the ",
  "table, and the basis for its selection is described in the Methods. No covariate differed across ",
  "tumor type except collection timing in its binned form, which reflects the arbitrary display bins of ",
  "Table 1: treated as the continuous variable it is, collection-year central tendency did not differ ",
  "across types (medians ", .med_phrase, "; fully overlapping ranges, ", .range_lo, "--", .range_hi, "). ",
  "The candidate confounders were also mutually independent."
)
writeLines(st4_caption, "Supplementary/Components/Tables/ST4_caption.tex")
#+ 17.6: ST5: Covariate-chemical association screen (reviewer #2 — do covariates track the chemicals)
#- 17.6.1: Reshape the per-covariate screen counts, with tumor type as an on-scale reference row
#! Effect size is quantitative-only: the two qualitative tests yield different metrics (Cramer's V, McFadden R-squared) that are not comparable
.st5_rows <- bind_rows(
  covariate_chemical_summary |> select(mode, covariate, n_sig, pct_sig, median_effect_sig),
  variant_benchmark |> select(mode, covariate, n_sig, pct_sig, median_effect_sig)
) |>
  mutate(cell = sprintf("%d (%.1f%%)", n_sig, pct_sig)) |>
  select(covariate, mode, cell, median_effect_sig) |>
  pivot_wider(names_from = mode, values_from = c(cell, median_effect_sig)) |>
  mutate(covariate = case_match(covariate,
    "year" ~ "Collection Year (Continuous)",
    "Sample_Collection_Timing" ~ "Collection Timing (Binned)",
    "Variant" ~ "Tumor Type",
    .default = covariate
  )) |>
  arrange(match(covariate, c("Collection Year (Continuous)", "Collection Timing (Binned)", "Sex", "Age", "Tumor Type"))) |>
  transmute(
    Covariate = covariate,
    `Significant|Quantitative` = cell_quant,
    `Significant|Qualitative` = cell_qual,
    `Median $\\eta^{2}$|Quantitative` = if_else(is.na(median_effect_sig_quant), "-", sprintf("%.3f", median_effect_sig_quant)),
    .section = FALSE
  )
#! Tumor type is split into its own section so it reads as the on-scale reference contrast rather than as another candidate confounder
ST5_data <- bind_rows(
  tibble(Covariate = "Candidate Confounders", .section = TRUE),
  .st5_rows |> filter(Covariate != "Tumor Type"),
  tibble(Covariate = "Reference", .section = TRUE),
  .st5_rows |> filter(Covariate == "Tumor Type")
) |>
  mutate(across(-c(Covariate, .section), \(x) replace_na(x, "")))
#- 17.6.2: Build ST5 LaTeX (self-contained tabular) and save
writeLines(build_supp_tabular(ST5_data), "Supplementary/Components/Tables/ST5.tex")
#- 17.6.3: Wire the caption deterministically from the screen summaries
.n_quant_tested <- max(covariate_chemical_summary$n_testable[covariate_chemical_summary$mode == "quant"])
.n_qual_tested <- max(covariate_chemical_summary$n_testable[covariate_chemical_summary$mode == "qual"])
st5_caption <- paste0(
  "\\textbf{Association of candidate confounders with individual chemical features.} ",
  "Every annotated feature was tested against each candidate confounder in both analytical modes ",
  "(", .n_quant_tested, " quantitative and ", .n_qual_tested, " qualitative features), and the number ",
  "of features significant at $P < 0.05$ is reported as \\textit{n} (\\%). Approximately 5\\% of features ",
  "are expected to reach this threshold by chance alone, so the relevant comparison for each covariate is ",
  "against that 5\\% baseline rather than against zero. Collection year appears twice: as a continuous ",
  "measure and as the categorical intervals used in Table 1. Tumor type is included as a reference row so ",
  "that the covariate associations can be read against the exposure contrast of interest using the same ",
  "metric; it is not a confounder, and the comparison is not adjusted for the differing degrees of freedom ",
  "consumed by each variable. Effect size is reported for quantitative features only, as a median across ",
  "the features reaching significance. Test selection and effect-size metrics are detailed in ",
  "Supplementary Note 2."
)
writeLines(st5_caption, "Supplementary/Components/Tables/ST5_caption.tex")
#+ 17.7: ST6: Covariate-adjusted tumor-type effects (reviewer #2 — do findings survive adjustment)
#- 17.7.1: Format p-values and effect sizes
.fmt_p6 <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3)))
.fmt_es6 <- function(e) ifelse(is.na(e), "-", formatC(e, format = "f", digits = 3))
#- 17.7.2: One block per mode, ordered by unadjusted significance within each
.st6_block <- function(md, header) {
  rows <- ancova_summary |>
    filter(mode == md) |>
    arrange(p_value_unadjusted) |>
    transmute(
      #! MT_final appends * for quality == 2; converted to the dagger/double-dagger convention used in S1/S2
      Chemical = paste0(
        str_remove(short_name, "\\*$"),
        if_else(str_detect(short_name, "\\*$"), "\\textsuperscript{\\textdagger}", ""),
        if_else(!is.na(n_detected) & n_detected <= 10, "\\textsuperscript{\\textdaggerdbl}", "")
      ),
      `n Detected` = if (md == "qual") formatC(n_detected, format = "d") else "-",
      #! Model 1/2/3 rather than spelled-out adjustment sets; the sets are defined in the caption and the numbering keeps the columns narrow
      `P|Model 1` = .fmt_p6(p_value_unadjusted),
      `P|Model 2` = .fmt_p6(p_value_year),
      `P|Model 3` = .fmt_p6(p_value_full),
      `Effect|Model 1` = .fmt_es6(effect_size_unadjusted),
      `Effect|Model 2` = .fmt_es6(effect_size_year),
      `Effect|Model 3` = .fmt_es6(effect_size_full),
      .section = FALSE
    )
  bind_rows(tibble(Chemical = header, .section = TRUE), rows)
}
#! Test and effect-size metric are stated in the caption rather than the section headers, which otherwise force the first column wide
ST6_data <- bind_rows(
  .st6_block("quant", "Quantitative features"),
  .st6_block("qual", "Qualitative features")
) |>
  mutate(across(-c(Chemical, .section), \(x) replace_na(x, "")))
#- 17.7.3: Build ST6 LaTeX (self-contained tabular) and save
writeLines(build_supp_tabular(ST6_data), "Supplementary/Components/Tables/ST6.tex")
#- 17.7.4: Wire the caption deterministically from ancova_summary
.n_total <- nrow(ancova_summary)
.n_binned <- sum(ancova_summary$survives_binned, na.rm = TRUE)
.sparse_min <- min(ancova_summary$n_detected, na.rm = TRUE)
st6_caption <- paste0(
  "\\textbf{Tumor-type differences in validated chemicals before and after covariate adjustment.} ",
  "Each of the ", .n_total, " validated type-differential chemicals was refit under three nested models: ",
  "Model 1, unadjusted; Model 2, adjusted for collection year; Model 3, adjusted for collection year, age, ",
  "and sex. Effect size is partial $\\eta^{2}$ for quantitative features and McFadden's pseudo-$R^{2}$ for ",
  "qualitative features, and should be compared only within a mode. Unadjusted \\textit{P} values for ",
  "qualitative features differ slightly from the exact tests reported in Table 3, because a ",
  "likelihood-ratio test is used throughout this table so that unadjusted and adjusted models remain ",
  "directly comparable; every chemical is significant under both. The chemicals that lost significance ",
  "under adjustment were marginal beforehand and are not among the findings emphasized in the manuscript. ",
  "In a sensitivity analysis adjusting instead for the categorical collection intervals of Table 1, ",
  .n_binned, " of ", .n_total, " remained significant; because that parameterization discards ",
  "within-interval variation it is reported here for completeness rather than as a primary specification. ",
  "Full model specification is given in Supplementary Note 2. ",
  "\\textit{$\\dagger$ Level 2 identification; $\\ddagger$ detected in 10 or fewer of the 60 samples ",
  "(minimum ", .sparse_min, "), for which covariate-adjusted estimates should be interpreted with caution.}"
)
writeLines(st6_caption, "Supplementary/Components/Tables/ST6_caption.tex")
