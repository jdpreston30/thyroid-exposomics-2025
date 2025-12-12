#+ 5.5: Create all variant graphs
#- 5.5.1: Bar plot for top 5 quant
p3A <- plot_top5_quant(tumors_quant_sig_top5)
#- 5.5.2: Heatmap for qualitative features
p3B <- plot_qualitative_heatmap(qual_i_reordered)
#- 5.5.3: Balloon plot for usage class vs variant
p3C <- plot_balloon(balloon_data_graph)
#+ 5.4: Prepare Balloon Plot Data
#- 5.4.1: Organize data
balloon_data <- MTii |>
  filter(!is.na(highest)) |>
  select(cas, usage_class, highest) |>
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
  # Remove 'Equal' and NA values
  filter(Variant != "Equal", !is.na(Variant)) |>
  # Count occurrences for each variant within each usage class
  count(usage_class, Variant) |>
  mutate(
    usage_class = str_replace_all(usage_class, "[†‡]", "")
  )
#- 5.4.2: Determine the correct order for usage_class
usage_class_order <- balloon_data |>
  pivot_wider(names_from = Variant, values_from = n, values_fill = 0) |>
  arrange(Papillary, `FV-PTC`, Follicular) |>
  pull(usage_class)
#- 5.4.3: Apply the correct order for y-axis
balloon_data_graph <- balloon_data |>
  mutate(
    Variant = factor(Variant, levels = c("Follicular", "FV-PTC", "Papillary")),
    usage_class = factor(usage_class, levels = usage_class_order)
  )