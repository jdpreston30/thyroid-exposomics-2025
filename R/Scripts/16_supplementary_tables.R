#* 16: Supplementary Tables
#+ 16.1: ST1: Chemical library (pivoted subid)
#- 16.1.0: Get list of Y/N expanded fragments
expanded_chemicals <- expanded_validation |>
  select(id, expanded) |>
  filter(expanded == "Y") |>
  pull(id)
#- 11.1.1: Prepare ST1 tibble
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
  # Final cleanup - no # column
  select(id, name, cas, monoisotopic, trt, starts_with("mz")) |>
  rename(
    `Library ID` = id,
    Name = name,
    CAS = cas,
    `Monoisotopic Mass` = monoisotopic,
    `Target RT (min.)` = trt
  )
#- 11.1.2: Build and format ST1 gt table
gt_ST1 <- build_ST1(ST1_tibble)
#- 11.1.3: Save ST1 as LaTeX
save_table_latex(gt_ST1, "Outputs/Tables/ST1.tex")
#+ 11.2: ST2
#+ Abbreviations Table
#- Build abbreviations table
abbrev_table <- ST1_abbrevs