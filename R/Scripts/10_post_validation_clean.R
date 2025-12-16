#* 10: Post Validation Cleaning
#+ 10.1: Subset to validated compounds
#- 10.1.1: Validated compounds
validated_compounds <- validation_check_files |>
  select(id, quality, cas, source, short_name) |>
  filter(quality %in% c(1, 2))
#- 10.1.2: Validated Variant Differences list
validated_variant <- validated_compounds |>
  filter(source %in% c("VD", "IARC and VD")) |>
  pull(id)
#- 10.1.2: Validated list
validated_iarc <- validated_compounds |>
  filter(source %in% c("IARC", "IARC and VD")) |>
  pull(id)
#+ 10.2: Modify MT_final per validation results
MT_final <- MT_final_i |>
  filter(id %in% validated_list) |>
  inner_join(validated_variant, by = "id") |>
  relocate(quality, .before = annot_ident)
#+ 10.3: Modify cadaver/tumor comparisons by validation results
#- 10.3.1: Subset controls to validated IARCs
IARC_controls <- IARC_controls_i |>
  filter(id %in% validated_iarc) |>
  left_join(validated_compounds |>
    select(id, short_name), by = "id") |>
  select(pct_NA_ctrl, id, short_name, name_sub_lib_id, iMean_ctrl, subid)
#- 10.3.2: Subset tumors to validated IARCs
IARC_tumors <- IARC_tumors_i |>
  filter(id %in% validated_iarc) |>
  left_join(validated_compounds |>
    select(id, short_name), by = "id") |>
  select(pct_NA_tumor = pct_NA, id, short_name, name_sub_lib_id, iMean_tumor = iMean, subid)
#- 10.3.3: Do a full join of the detection percentages of both
IARC_combined <- IARC_controls |>
  full_join(IARC_tumors, by = c("name_sub_lib_id", "id", "subid", "short_name"), suffix = c("_ctrl", "_tumor")) |>
    mutate(
    pct_NA_tumor = replace_na(pct_NA_tumor, 1),
    pct_NA_ctrl = replace_na(pct_NA_ctrl, 1)
  ) |>
  arrange(short_name, pct_NA_tumor) |>
  select(short_name, id, subid, pct_NA_tumor, pct_NA_ctrl, iMean_ctrl, iMean_tumor)
#- 10.3.2: full_joiner_i
full_joiner <- full_joiner_i |>
  select(sample_ID, variant, tumor_vs_ctrl, any_of(cadaver_iarc_keep))
#+ 10.4: IARC detection heatmap
#- 10.4.1: Prepare data for heatmap
IARC_heatmap_data <- IARC_combined |>
  mutate(compound_var = paste0(short_name, "_", subid)) |>
  mutate(compound_var = factor(compound_var, levels = unique(compound_var))) |>
  pivot_longer(cols = c(pct_NA_tumor, pct_NA_ctrl), 
               names_to = "tissue_type", 
               values_to = "pct_NA") |>
  mutate(tissue_type = case_when(
    tissue_type == "pct_NA_tumor" ~ "Tumor",
    tissue_type == "pct_NA_ctrl" ~ "Control"
  ))
#- 10.4.2: Create heatmap
IARC_heatmap <- ggplot(IARC_heatmap_data, aes(x = tissue_type, y = compound_var, fill = pct_NA * 100)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 50, limits = c(0, 100),
                       name = "% Non-detection") +
  labs(x = NULL, y = NULL, title = "IARC Group 1 Chemical Detection") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())
#+ 10.5: 
