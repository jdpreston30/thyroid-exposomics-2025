#* 11: Visualization of Variant Differences Post-Validation
#+ 11.1: Top 5 Quant Graph
#- 11.1.1: Pull fragments for top 5 in order
top_5_ordered <- MT_final |>
  filter(mode == "quantitative") |>
  arrange(p_value) |>
  slice_head(n = 5) |>
  pull(name_sub_lib_id)
#- 11.1.2: Subset data to those with correct column order
tumors_quant_sig_top5 <- tumors_quant_sig |>
  select(variant, all_of(top_5_ordered))
#- 11.1.3: Create plot
p3A <- plot_top5_quant(tumors_quant_sig_top5, add_cld = TRUE)
#+ 11.2: Qualitative Features Heatmap

#+ 11.3: Balloon Plot
#- 11.3.1: Organize data
balloon_data <- MT_final |>
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
    usage_class = str_replace_all(usage_class, "[†‡ᵃᵇ]", "")
  )
#- 11.3.2: Get summary
balloon_data %>%
  group_by(Variant) %>%
  summarise(total = sum(n)) %>%
  arrange(desc(total))
#- 11.3.3: Determine the correct order for usage_class
usage_class_order <- balloon_data |>
  pivot_wider(names_from = Variant, values_from = n, values_fill = 0) |>
  arrange(Papillary, `FV-PTC`, Follicular) |>
  pull(usage_class)
#- 11.3.4: Apply the correct order for y-axis
balloon_data_graph <- balloon_data |>
  mutate(
    Variant = factor(Variant, levels = c("Follicular", "FV-PTC", "Papillary")),
    usage_class = factor(usage_class, levels = usage_class_order)
  )
#- 11.3.5: Make Balloon Plot
p3C <- plot_balloon(balloon_data_graph, show_x_labels = TRUE)
#+ 11.4: Qualitative Features Heatmap
#- 11.4.1: Get qualitative features data
qualitative_final_features <- MT_final |>
  filter(mode == "qualitative") |>
  select(short_name, FTC, FV_PTC, PTC)
#- 11.4.2: Prepare data for heatmap (convert % to numeric)
qualitative_heatmap_data <- qualitative_final_features |>
  pivot_longer(cols = c(FTC, FV_PTC, PTC), 
               names_to = "variant", 
               values_to = "pct_detection") |>
  mutate(
    pct_detection = as.numeric(gsub("%", "", pct_detection)),
    variant = factor(variant, 
                    levels = c("FTC", "FV_PTC", "PTC"),
                    labels = c("Follicular", "FV-PTC", "Papillary"))
  )
#- 11.4.3: Order compounds from least to greatest by Papillary, then FV-PTC, then Follicular
compound_order <- qualitative_final_features |>
  mutate(across(c(FTC, FV_PTC, PTC), ~as.numeric(gsub("%", "", .)))) |>
  arrange(PTC, FV_PTC, FTC) |>
  pull(short_name)
#- 11.4.4: Apply order
qualitative_heatmap_data <- qualitative_heatmap_data |>
  mutate(short_name = factor(short_name, levels = compound_order))
#- 11.4.5: Create heatmap
p3B <- plot_qualitative_heatmap(qualitative_heatmap_data)
