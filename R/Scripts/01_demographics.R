#* 1: Demographics
#+ 1.1: Compute Mean age
age <- demographics |>
  group_by(variant) |>
  summarise(Mean_Stdev = paste0(round(mean(age_collection), 1), " ± ", round(sd(age_collection), 1))) |>
  pivot_wider(names_from = variant, values_from = Mean_Stdev) |>
  mutate(Variable = "Age") |>
  bind_cols(
    tibble(Total = paste0(round(mean(demographics$age_collection), 1), " ± ", round(sd(demographics$age_collection), 1)))
  ) |>
  select(Variable, Follicular, FVPTC, Papillary, Total)
#+ 1.2: Summarize Sex
sex <- demographics |>
  group_by(variant) |>
  summarise(n_female = sum(sex == "Female"), total = n(), .groups = "drop") |>
  mutate(
    percent = round((n_female / total) * 100, 1),
    count_percent = paste0(n_female, " (", percent, "%)")
  ) |>
  select(variant, count_percent) |>
  pivot_wider(names_from = variant, values_from = count_percent) |>
  mutate(Variable = "Sex (Female)") |>
  bind_cols(
    demographics |>
      summarise(n_female = sum(sex == "Female"), total = n()) |>
      mutate(Total = paste0(n_female, " (", round((n_female / total) * 100, 1), "%)")) |>
      select(Total)
  ) |>
  select(Variable, Follicular, FVPTC, Papillary, Total)
#+ 1.3: Determine Collection bins
bins <- demographics |>
  group_by(variant, year_bin) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(variant) |>
  mutate(percent = round((n / sum(n)) * 100, 1)) |>
  ungroup() |>
  mutate(count_percent = paste0(n, " (", percent, "%)")) |>
  select(year_bin, variant, count_percent) |>
  pivot_wider(names_from = variant, values_from = count_percent) |>
  mutate(Variable = year_bin) |>
  left_join(
    demographics |>
      count(year_bin) |>
      mutate(Total = paste0(n, " (", round((n / sum(n)) * 100, 1), "%)")) |>
      select(year_bin, Total),
    by = "year_bin"
  ) |>
  select(Variable, Follicular, FVPTC, Papillary, Total) |>
  arrange(Variable)
#+ 1.4: Combine and Export
demographic_table <- rbind(age, sex, bins)
write.csv(demographic_table, "Outputs/Tables/T1.csv")
