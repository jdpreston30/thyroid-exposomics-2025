#* 12: Tumor v Control IARC Plots
#+ 12.1: IARC Stats (ttest on log transformed)
IARC_ttests <- full_joiner |>
  mutate(across(where(is.numeric), ~ log2(. + 1e-6))) |>
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
p3D <- plot_iarc(toluidine_data, chemical_name = "o-Toluidine", p_value = toluidine_p)
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
p3E <- plot_iarc(aminobiphenyl_0_data, chemical_name = "4-Aminobiphenyl", p_value = aminobiphenyl_0_p)
#+ 12.5: IARC detection heatmap
#! Leaving out for now but may revisit if individual check of top fragments proceeds