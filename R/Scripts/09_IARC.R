#* 9: Tumor v Control IARC Analysis
#+ 9.1: IARC Stats (ttest on log transformed)
IARC_ttests <- full_joiner |>
  select(sample_ID, any_of(IARC_tumors_ctrl_filtered)) |>
  mutate(across(
    .cols = -c(sample_ID),
    .fns = ~ {
      col_min <- min(.x, na.rm = TRUE)
      replace(.x, is.na(.x), 0.5 * col_min)
    }
  )) |>
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
#+ 9.2: Graph tumor v cadaver IARC1
#- 9.2.1: Get p-value from IARC_ttests
toluidine_p <- IARC_ttests |> 
  filter(chemical == "o-Toluidine_0_BP3.GC2_CP3017") |> 
  pull(p_value)
#- 9.2.2: Subset plot data plot
p3D_data <- full_joiner %>%
  select(tumor_vs_ctrl, `o-Toluidine_0_BP3.GC2_CP3017`) %>%
  rename(concentration = `o-Toluidine_0_BP3.GC2_CP3017`)
#- 9.2.3: Create plot
{
source("R/Utilities/Visualization/plot_iarc.R")
p3D <- plot_iarc(p3D_data, chemical_name = "o-Toluidine", p_value = toluidine_p)
source("temp3.R")
}
