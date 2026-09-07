#* 12: Tumor v Control IARC Plots
#+ 12.1: IARC Stats (ttest on Log transformed)
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
#+ 12.2: Graph Tumor v Cadaver IARC1
#! FRAGMENT CHOICE -- both panels hardcode the _0 fragment because it is the one the source data flags as the quantifier: lib.subject.qsummary carries bestFlag = 1 AND bestQuantFlag = 1 on o-Toluidine_0 and 4-aminobiphenyl_0, and the cadaver sheet agrees. The criterion is detection then intensity, NOT p-value -- 4-aminobiphenyl_0 is detected in 60/60 tumors and 8/8 cadavers (detectFrac 1.00) versus _1 at 55/60 and _2 at 52/60. So _0 is chosen even though _1 has the smaller P (7.18e-17 vs 1.33e-11); picking on P would be selecting the fragment by its own result.
#! This mirrors the rule still visible at 06_tumor_cadaver.R:131, `arrange(desc(pct_det_tumor), desc(pct_det_ctrl), desc(iMean_tumors))`. That arrange was once followed by `group_by(cas) |> slice_head(n = 1)`, which auto-selected the ideal fragment; it was removed in c7ecd9c (2025-12-16) so full_joiner could retain every fragment. These hardcoded names are that removed selector's output, written out by hand -- 8dc27ec, the commit that introduced them, is from the same day.
#- 12.2.1: o-Toluidine_0_BP3.GC2_CP3017
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
      "FV_PTC" = "IEFVPTC",
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
#- 12.3.3: Summary of updated carcinogen classification
carc_summary <- MTi |>
  filter(cas %in% MT_final_cas_list) |>
  filter(Carcinogenicity != "Unclassified") |>
  # Expand rows for double counting when multiple highest groups
  mutate(
    highest = strsplit(highest, ", ")
  ) |>
  unnest(highest) |>
  select(short_name, Carcinogenicity, highest) |>
  unique() |>
  # Map the group names to the variant labels
  mutate(
    Variant = recode(highest,
      "FTC" = "Follicular",
      "FV_PTC" = "IEFVPTC",
      "PTC" = "Papillary"
    )
  ) |>
  # Count occurrences for each Variant and Carcinogenicity
  count(Variant, Carcinogenicity)
#+ 12.4: IARC Detection Heatmap
#! Leaving out for now but may revisit if individual check of top fragments proceeds
