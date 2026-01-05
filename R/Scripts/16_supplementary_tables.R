#* 16: Supplementary Tables
#+ 16.1: ST1: Chemical library (pivoted subid)
#- 16.1.0: Get list of Y/N expanded fragments
expanded_chemicals <- expanded_validation |>
  select(id, expanded) |>
  filter(expanded == "Y") |>
  pull(id)
#- 16.1.1: Prepare ST1 tibble
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
  # Add dagger superscript for expanded chemicals
  mutate(
    name = if_else(
      id %in% expanded_chemicals & name != "-",
      paste0(name, "$^\\dagger$"),
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
    `Target RT (min.)` = trt
  )
#- 16.1.2: Build and format ST1 gt table
gt_ST1 <- build_ST1(ST1_tibble)
#- 16.1.3: Save ST1 as LaTeX (without table wrapper) to Supplementary/Components/Tables
latex_code <- gt::as_latex(gt_ST1) |> as.character()
latex_code <- fix_latex_header_fill(latex_code)
# Remove table wrapper for direct inclusion in supplementary
latex_lines <- strsplit(latex_code, "\n")[[1]]
latex_lines <- latex_lines[-c(1, length(latex_lines))]  # Remove \begin{table} and \end{table}
latex_code <- paste(latex_lines, collapse = "\n")
writeLines(latex_code, "Supplementary/Components/Tables/ST1.tex")
#+ 16.2: ST2
#- 16.2.1: Prepare base ST2 data 
ST2_base <- feature_metadata |>
  mutate(`Superclass: Class` = paste(Superclass, Class, sep = ": ")) |>
  select(GROUP = Table_Header, Class = Table_Class, Subclass = Table_Subclass, CAS = cas, `Potential EDC` = Potential_EDC, `IARC Group` = IARC_Group, `Superclass: Class`) |>
  mutate(GROUP = toupper(GROUP)) |>
  left_join(ST1_tibble |> select(Name, CAS), by = "CAS") |>
  mutate(Name = gsub("\\$\\^\\\\dagger\\$", "", Name)) |>
  arrange(GROUP, Class, Subclass, Name)
#- 16.2.2: Build hierarchical structure with proper nesting
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
#- 16.2.3: Combine and format
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
#- 16.2.4: Build and format ST2 gt table
gt_ST2 <- build_ST2(ST2_tibble)
#- 16.2.5: Save ST2 as LaTeX (without table wrapper) to Supplementary/Components/Tables
latex_code <- gt::as_latex(gt_ST2) |> as.character()
latex_code <- fix_ST2_latex(latex_code)
# Remove table wrapper for direct inclusion in supplementary
latex_lines <- strsplit(latex_code, "\n")[[1]]
latex_lines <- latex_lines[-c(1, length(latex_lines))]  # Remove \begin{table} and \end{table}
latex_code <- paste(latex_lines, collapse = "\n")
writeLines(latex_code, "Supplementary/Components/Tables/ST2.tex")
#+ 16.3: ST3
#- 16.3.0: Vector of selected CAS
selected_ST3_cas <- c(
  "83-32-9", "208-96-8", "120-12-7", "56-55-3", "205-99-2", "207-08-9", "50-32-8", "218-01-9",
  "53-70-3", "206-44-0", "91-20-3", "85-01-8", "129-00-0", "91-59-8", "92-67-1", "92-87-5",
  "95-53-4", "32598-14-4", "35065-28-2", "35065-27-1", "52663-74-8", "35065-29-3", "52663-68-0"
)
#- 16.3.1: Create ST3 tibble (initial)
ST3_tibble_i <- ppm_full_table |>
  filter(!is.na(IARC_Group)) |>
  select(CAS = cas, everything(), -Name) |>
  mutate(
    across(c(mean_ctrl, mean_tumor, half_min, max_value), 
           ~ . * 10^3, 
           .names = "{.col}")
  ) |>
  mutate(
    # Format with commas, then handle < 1 cases
    across(c(mean_ctrl, mean_tumor),
           ~ if_else(. < 1, "< 1", format(round(.), big.mark = ",")),
           .names = "{.col}"),
    # For range components, format separately
    half_min = if_else(half_min < 1, "< 1", format(round(half_min), big.mark = ",")),
    max_value = if_else(max_value < 1, "< 1", format(round(max_value), big.mark = ","))
  ) |>
  left_join(ST1_tibble |> select(Name, CAS), by = "CAS") |>
  mutate(range = paste0(half_min, "-", max_value)) |>
  group_by(CAS) |>
  mutate(
    best_fragment = if_else(pct_det_tumor == max(pct_det_tumor), "★", ""),
    tumor_ctrl_det = paste0(round(pct_det_tumor * 100), "/", round(pct_det_ctrl * 100))
  ) |>
  ungroup()
#- 16.3.2: Create inspection tibble; filter drop bad ones
ST3_tibble_inspect <- ST3_tibble_i |>
  select(name_sub_lib_id, best_fragment, tumor_ctrl_det, range, IARC_Group, CAS, everything()) |>
  filter(CAS %in% selected_ST3_cas) |>
  filter(!(pct_det_tumor < 0.5 & pct_det_ctrl < 0.5)) |>
  group_by(CAS) |>
  mutate(
    Tie = if_else(sum(best_fragment == "★") > 1, "Y", "")
  ) |>
  ungroup() |>
  arrange(name_sub_lib_id)
#- 16.3.3: Set tiebreaker vector based on manual inspection
tiebreak_ST3_fragments <- c(
  "Pyrene_3_BP3.GC2_CP3078",
  "Fluoranthene_3_BP3.GC2_CP3057",
  "Chrysene_0_BP3.GC2_CP3044"
)
#- 16.3.4: Near final ST3 tibble using tiebreaker
ST3_tibble_2 <- ST3_tibble_inspect |>
  filter(best_fragment == "★") |>
  filter((Tie == "Y" & name_sub_lib_id %in% tiebreak_ST3_fragments) | Tie == "")
#- 16.3.5: Determine if any outright failed
  ST3_fail_check <- validation_check_files_unfiltered |>
    select(id, cas, short_name, quality, state, modification) |>
    arrange(quality) |>
    arrange(short_name) |>
    filter(cas %in% ST3_tibble_2$CAS)
#- 16.3.6: Create a vector of failed chemicals
failed_ST3 <- c(
  "50-32-8", # Benzo[a]pyrene
  "207-08-9", #Benzo(k)fluoranthene
  # "218-01-9" # Chrysene, was going to include but validation showed failure in part 2 standard and part 3 may still be viable
  # "91-20-3", # Naphthalene, was originally going to include but part 3 was failure plot and not sure if part 2 ever inspected so going to include the part 2 one
  "92-87-5" # Benzidine
)
#- 16.3.7: Finalize ST3 tibble with hierarchical structure
ST3_tibble_flat <- ST3_tibble_2 |>
  filter(!CAS %in% failed_ST3) |>
  left_join(literature_ST3 |> select(CAS, `Usage Class`, AT_manuscript, AT_ref, plasma_manuscript, plasma_ref, urine_manuscript, urine_ref), by = "CAS") |>
  mutate(
    Name = gsub("\\$\\^\\\\dagger\\$|\\*", "", Name),
    `Range (PPB)` = range,
    `Mean Non-Cancer Thyroid Concentration (PPB)` = mean_ctrl,
    `Mean Tumor Concentration (PPB)` = mean_tumor,
    # Combine literature values with superscripted reference numbers
    # Reference mapping: mlyczynska2023=11, maier2022=12, riffelmann1995=13, cdc2024=14, wang2010=16
    `Adipose Tissue (PPB)†` = case_when(
      !is.na(AT_manuscript) & !is.na(AT_ref) ~ paste0(AT_manuscript, "\\textsuperscript{", 
        case_when(
          AT_ref == "mlyczynska2023" ~ "11",
          AT_ref == "maier2022" ~ "12",
          AT_ref == "riffelmann1995" ~ "13",
          AT_ref == "cdc2024" ~ "14",
          AT_ref == "wang2010" ~ "16",
          TRUE ~ AT_ref
        ), "}"),
      !is.na(AT_manuscript) ~ AT_manuscript,
      TRUE ~ NA_character_
    ),
    `Urine (PPB)‡` = case_when(
      !is.na(urine_manuscript) & !is.na(urine_ref) ~ paste0(urine_manuscript, "\\textsuperscript{", 
        case_when(
          urine_ref == "mlyczynska2023" ~ "11",
          urine_ref == "maier2022" ~ "12",
          urine_ref == "riffelmann1995" ~ "13",
          urine_ref == "cdc2024" ~ "14",
          urine_ref == "wang2010" ~ "16",
          TRUE ~ urine_ref
        ), "}"),
      !is.na(urine_manuscript) ~ urine_manuscript,
      TRUE ~ NA_character_
    ),
    `Serum/Plasma (PPB)‡` = case_when(
      !is.na(plasma_manuscript) & !is.na(plasma_ref) ~ paste0(plasma_manuscript, "\\textsuperscript{", 
        case_when(
          plasma_ref == "mlyczynska2023" ~ "11",
          plasma_ref == "maier2022" ~ "12",
          plasma_ref == "riffelmann1995" ~ "13",
          plasma_ref == "cdc2024" ~ "14",
          plasma_ref == "wang2010" ~ "16",
          TRUE ~ plasma_ref
        ), "}"),
      !is.na(plasma_manuscript) ~ plasma_manuscript,
      TRUE ~ NA_character_
    )
  ) |>
  select(
    `Usage Class`,
    Name,
    CAS,
    `IARC Group` = IARC_Group,
    `Mean Non-Cancer Thyroid Concentration (PPB)` = `Mean Non-Cancer Thyroid Concentration (PPB)`,
    `Mean Tumor Concentration (PPB)` = `Mean Tumor Concentration (PPB)`,
    `Range (PPB)*` = `Range (PPB)`,
    `Adipose Tissue (PPB)†` = `Adipose Tissue (PPB)†`,
    `Urine (PPB)‡` = `Urine (PPB)‡`,
    `Serum/Plasma (PPB)‡` = `Serum/Plasma (PPB)‡`
  ) |>
  arrange(`Usage Class`, Name)
#- 16.3.8: Build hierarchical structure with Usage Class headers
ST3_list <- list()
all_usage_classes <- unique(ST3_tibble_flat$`Usage Class`)
# Sort usage classes alphabetically
all_usage_classes <- sort(all_usage_classes)
for (usage_idx in seq_along(all_usage_classes)) {
  usage_class <- all_usage_classes[usage_idx]
  # Add Usage Class header row (in Name column)
  ST3_list[[length(ST3_list) + 1]] <- tibble(
    Name = usage_class,
    CAS = NA_character_,
    `IARC Group` = NA_character_,
    `Mean Non-Cancer Thyroid Concentration (PPB)` = NA_character_,
    `Mean Tumor Concentration (PPB)` = NA_character_,
    `Range (PPB)*` = NA_character_,
    `Adipose Tissue (PPB)†` = NA_character_,
    `Urine (PPB)‡` = NA_character_,
    `Serum/Plasma (PPB)‡` = NA_character_
  )
  # Add chemicals for this usage class with 2 spaces indent (already sorted by Name)
  usage_data <- ST3_tibble_flat |> 
    filter(`Usage Class` == usage_class) |>
    mutate(Name = paste0("  ", Name)) |>  # Add 2-space indent for chemicals
    select(-`Usage Class`)  # Remove Usage Class column
  ST3_list[[length(ST3_list) + 1]] <- usage_data
  # Add blank spacer row after each usage class (except the last)
  if (usage_idx < length(all_usage_classes)) {
    ST3_list[[length(ST3_list) + 1]] <- tibble(
      Name = "",
      CAS = "",
      `IARC Group` = "",
      `Mean Non-Cancer Thyroid Concentration (PPB)` = "",
      `Mean Tumor Concentration (PPB)` = "",
      `Range (PPB)*` = "",
      `Adipose Tissue (PPB)†` = "",
      `Urine (PPB)‡` = "",
      `Serum/Plasma (PPB)‡` = ""
    )
  }
}
#- 16.3.9: Combine into final ST3 tibble
ST3_tibble <- bind_rows(ST3_list) |>
  mutate(
    # Add en-dash for missing literature values in chemical rows (not headers or spacers)
    across(
      c(`Adipose Tissue (PPB)†`, `Urine (PPB)‡`, `Serum/Plasma (PPB)‡`),
      ~ case_when(
        Name != "" & !is.na(CAS) & (is.na(.) | . == "") ~ "–",
        TRUE ~ .
      )
    )
  )
#- 16.3.10: Build and format ST3 gt table
gt_ST3 <- build_ST3(ST3_tibble)
#- 16.3.11: Save ST3 as LaTeX (without table wrapper) to Supplementary/Components/Tables
latex_code <- gt::as_latex(gt_ST3) |> as.character()
latex_code <- fix_ST3_latex(latex_code)
# Remove table wrapper for direct inclusion in supplementary
latex_lines <- strsplit(latex_code, "\n")[[1]]
latex_lines <- latex_lines[-c(1, length(latex_lines))]  # Remove \begin{table} and \end{table}
latex_code <- paste(latex_lines, collapse = "\n")
writeLines(latex_code, "Supplementary/Components/Tables/ST3.tex")
#+ 16.4: ST4
#! Pending
#+ 16.5 Abbreviations Dictionary
#- 11.5.1: Build abbreviations list
abbrev_list <- ST_abbrevs |>
  arrange(formatted)
#- 11.5.2: Convert to LaTeX itemize list with reduced spacing
abbrev_latex <- c(
  "\\begin{itemize}",
  "\\setlength{\\itemsep}{2pt}",
  "\\setlength{\\parskip}{0pt}",
  "\\setlength{\\parsep}{0pt}",
  paste0("  \\item ", abbrev_list$formatted),
  "\\end{itemize}"
)
#- 11.5.3: Save abbreviations to Supplementary/Components/Sections
writeLines(abbrev_latex, "Supplementary/Components/Sections/abbreviations.tex")