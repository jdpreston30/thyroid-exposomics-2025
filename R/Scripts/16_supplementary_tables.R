#* 16: Supplementary Tables
#+ 16.1: ST1: Chemical Library (pivoted subid)
#- 16.1.1: Build and format ST1 gt table
gt_ST1 <- build_ST1(ST1_tibble)
#- 16.1.2: Save ST1 as LaTeX (without table wrapper) to Supplementary/Components/Tables
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
  mutate(`Superclass: Class` = if_else(is.na(Superclass), Class, paste(Superclass, Class, sep = ": "))) |>
  select(GROUP = Table_Header, Class = Table_Class, Subclass = Table_Subclass, CAS = cas, `Potential EDC` = Potential_EDC, `IARC Group` = IARC_Group, `Superclass: Class`) |>
  mutate(GROUP = toupper(GROUP)) |>
  left_join(ST1_tibble |> select(Name, CAS), by = "CAS") |>
  mutate(Name = gsub("\u2021", "", Name)) |>
  mutate(Name = gsub("\\*", "\u2020", Name)) |>
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
#+ 16.3: Abbreviations Dictionary
#- 16.3.1: Build abbreviations list
abbrev_list <- ST_abbrevs |>
  arrange(formatted)
#- 16.3.2: Convert to LaTeX itemize list with reduced spacing
abbrev_latex <- c(
  "\\begin{itemize}",
  "\\setlength{\\itemsep}{2pt}",
  "\\setlength{\\parskip}{0pt}",
  "\\setlength{\\parsep}{0pt}",
  paste0("  \\item ", abbrev_list$formatted),
  "\\end{itemize}"
)
#- 16.3.3: Save abbreviations to Supplementary/Components/Sections
writeLines(abbrev_latex, "Supplementary/Components/Sections/abbreviations.tex")
#+ 16.4: ST3: Cohort demographics (DTC tumors vs cadaver controls)
#- 16.4.1: Assemble combined demographic frame (age, sex, sample collection timing)
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
#- 16.4.2: Summarize with TernTables (returns display tibble; also writes Word version)
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
#- 16.4.3: Build ST3 LaTeX (self-contained tabular) and save to Supplementary/Components/Tables
st3_latex <- build_ST3(ST3_tern)
writeLines(st3_latex, "Supplementary/Components/Tables/ST3.tex")
