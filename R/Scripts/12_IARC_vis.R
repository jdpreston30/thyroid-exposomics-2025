#* 12: Tumor v Control IARC Plots
#+ 12.1: IARC Stats (ttest on log transformed)
IARC_ttests <- full_joiner |>
  mutate(across(where(is.numeric), ~ log2(.))) |>
  mutate(group = if_else(str_starts(sample_ID, "T00"), "Control", "Tumor")) |>
  pivot_longer(
    cols = where(is.numeric) & !any_of("group"),
    names_to = "chemical", values_to = "value"
  ) |>
  group_by(chemical) |>
  summarise(
    p_value = t.test(value ~ group)$p.value,
    mean_control = mean(value[group == "Control"]),
    mean_tumor = mean(value[group == "Tumor"]),
    .groups = "drop"
  ) |>
  arrange(p_value)
#+ 12.2: Graph tumor v cadaver IARC1
#- 12.2.4: o-Toluidine_0_BP3.GC2_CP3017
# P-value
toluidine_p <- IARC_ttests |> 
  filter(chemical == "o-Toluidine_0_BP3.GC2_CP3017") |> 
  pull(p_value)
# Data
toluidine_data <- full_joiner %>%
  select(tumor_vs_ctrl, `o-Toluidine_0_BP3.GC2_CP3017`) %>%
  rename(concentration = `o-Toluidine_0_BP3.GC2_CP3017`)
p3E <- plot_iarc(toluidine_data, chemical_name = "o-Toluidine", p_value = toluidine_p)
#- 12.2.2: 4-aminobiphenyl_0_BP3.GC2_CP3002
#! This was the most detected and top frag
# P-value
aminobiphenyl_0_p <- IARC_ttests |> 
  filter(chemical == "4-aminobiphenyl_0_BP3.GC2_CP3002") |> 
  pull(p_value)
# Data
aminobiphenyl_0_data <- full_joiner %>%
  select(tumor_vs_ctrl, `4-aminobiphenyl_0_BP3.GC2_CP3002`) %>%
  rename(concentration = `4-aminobiphenyl_0_BP3.GC2_CP3002`)
p3F <- plot_iarc(aminobiphenyl_0_data, chemical_name = "4-Aminobiphenyl", p_value = aminobiphenyl_0_p)
#+ 12.3: Advanced Carcinogen Classificaiton
#- 12.3.1: Advanced carcinogen classification run for variants
carc_by_variant <- MTi |>
  filter(cas %in% MT_final_cas_list) |>
  # Expand rows for double counting when multiple highest groups
  mutate(
    highest = strsplit(highest, ", ")
  ) |>
  unnest(highest) |>
  # Map the group names to the variant labels
  mutate(
    Variant = recode(highest,
      "FTC" = "Follicular",
      "FV_PTC" = "FV-PTC",
      "PTC" = "Papillary"
    )
  ) |>
  # Count occurrences for each Variant and Carcinogenicity
  count(Variant, Carcinogenicity) |>
  # Pivot to wide format
  pivot_wider(
    names_from = Carcinogenicity,
    values_from = n,
    values_fill = 0 # Fill missing counts with 0
  ) |>
  # Reorder columns for readability and exclude Unclassified
  select(
    Variant,
    "Known Carcinogen",
    "Likely Carcinogen",
    "Possible Carcinogen",
    "Uncertain Risk"
  ) |>
  filter(Variant != "Equal")
#- 12.3.2: Create carcinogen classification stacked bar plot
p3D <- plot_carcinogen_stacked(carc_by_variant)
#- 13.3.3: Summary of updated carcinogen classification
carc_summary <- MTi |>
  filter(cas %in% MT_final_cas_list) |>
  filter(Carcinogenicity != "Unclassified") |>
  # Expand rows for double counting when multiple highest groups
  mutate(
    highest = strsplit(highest, ", ")
  ) |>
  unnest(highest) |>
  select(short_name, Carcinogenicity) |>
  unique()
  # Map the group names to the variant labels
  mutate(
    Variant = recode(highest,
      "FTC" = "Follicular",
      "FV_PTC" = "FV-PTC",
      "PTC" = "Papillary"
    )
  ) |>
  # Count occurrences for each Variant and Carcinogenicity
  count(Variant, Carcinogenicity)
#+ 12.5: IARC detection heatmap
#! Leaving out for now but may revisit if individual check of top fragments proceeds