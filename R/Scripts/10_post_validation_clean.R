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
#- 10.1.2: Validated iarc
validated_iarc <- validated_compounds |>
  filter(source %in% c("IARC", "IARC and VD")) |>
  pull(id)
#- 10.1.2: Validated variant joiner tibble
validated_list <- validated_compounds |>
  select(id, quality)
#+ 10.2: Modify MT_final per validation results
MT_final <- MT_final_i |>
  filter(id %in% validated_variant) |>
  inner_join(validated_list, by = "id") |>
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
#- 10.3.4: Do a full join of the detection percentages of both
IARC_combined <- IARC_controls |>
  full_join(IARC_tumors, by = c("name_sub_lib_id", "id", "subid", "short_name"), suffix = c("_ctrl", "_tumor")) |>
  mutate(
    pct_detected_tumor = 1 - replace_na(pct_NA_tumor, 1),
    pct_detected_ctrl = 1 - replace_na(pct_NA_ctrl, 1)
  ) |>
  filter(pct_detected_tumor >= 0.7 & pct_detected_ctrl >= 0.7) |>
  arrange(short_name, desc(pct_detected_tumor)) |>
  select(name_sub_lib_id, short_name, id, subid, pct_detected_tumor, pct_detected_ctrl, iMean_ctrl, iMean_tumor) |>
  arrange(desc(pct_detected_tumor), desc(pct_detected_ctrl))
#- 10.3.5: Create list of those remaining IARCs
IARC_namesub_pull <- IARC_combined |>
  pull(name_sub_lib_id)
#- 10.3.2: Use these to complete full joiner (ppm ppb estimates)
full_joiner <- full_joiner_i |>
  select(sample_ID, variant, tumor_vs_ctrl, any_of(IARC_namesub_pull))



#!!!!!!!!
#- 6.5.4: Determine the matches, filter to 30% NA quant, pull features
IARC_tumors_ctrl_filtered_id_frag <- IARC_tumors_i |>
  left_join(IARC_controls_i, by = "name_sub_lib_id") |>
  select(cas, name, name_sub_lib_id, id = id.x, pct_NA, pct_NA_ctrl, iMean, iMean_ctrl) |>
  filter(!is.na(pct_NA_ctrl)) |>
  group_by(cas) |>
  arrange(pct_NA, pct_NA_ctrl, desc(iMean)) |>
  slice(1) |>
  ungroup() |>
  filter(pct_NA < 0.3, pct_NA_ctrl < 0.3)
#- 6.5.5: Pull the name_sub_lib_id
IARC_tumors_ctrl_filtered <- IARC_tumors_ctrl_filtered_id_frag |>
  pull(name_sub_lib_id)
#- 6.5.5: Pull ID
IARC_tumors_ctrl_filtered_id <- IARC_tumors_ctrl_filtered_id_frag |>
  pull(id)
#!!!!!!!!!!

#+ 10.4: Full Joiner (Quant)

#!!!!!!!!!!!!!!!!!!!!
#+ 10.5: IARC detection heatmap
#- 10.4.0: Consolidate percentages per id
IARC_consolidated <- IARC_combined_top |>
  group_by(short_name, id) |>
  summarise(
    pct_NA_tumor = mean(pct_NA_tumor, na.rm = TRUE),
    pct_NA_ctrl = mean(pct_NA_ctrl, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(pct_NA_tumor))
#- 10.4.1: Prepare data for heatmap
IARC_heatmap_data <- IARC_consolidated |>
  mutate(short_name = factor(short_name, levels = unique(short_name))) |>
  pivot_longer(cols = c(pct_NA_tumor, pct_NA_ctrl), 
               names_to = "tissue_type", 
               values_to = "pct_NA") |>
  mutate(
    pct_detection = (1 - pct_NA) * 100,
    tissue_type = case_when(
      tissue_type == "pct_NA_tumor" ~ "Tumor",
      tissue_type == "pct_NA_ctrl" ~ "Control"
    )
  )
#- 10.4.2: Create heatmap
IARC_heatmap <- ggplot(IARC_heatmap_data, aes(x = tissue_type, y = short_name, fill = pct_detection)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#0B5DA4", mid = "#dfdbd7ff", high = "#BF2E39", 
                       midpoint = 50, limits = c(0, 100),
                       name = "% Detection") +
  labs(x = NULL, y = NULL, title = "IARC Group 1 Chemical Detection") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())
#+ 10.5: 
