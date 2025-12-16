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
IARC_tumors_i |>
  filter(id == "CP2535")
IARC_tumors <- IARC_tumors_i |>
  filter(id %in% validated_iarc) |>
  left_join(validated_compounds |>
    select(id, short_name), by = "id") |>
  select(pct_NA_tumor = pct_NA, id, short_name, name_sub_lib_id, iMean_tumor = iMean, subid)
#- 10.3.4: Do a full join of the detection percentages of both
IARC_combined <- IARC_controls |>
  full_join(IARC_tumors, by = c("name_sub_lib_id", "id", "subid", "short_name"), suffix = c("_ctrl", "_tumor")) |>
    mutate(
    pct_NA_tumor = replace_na(pct_NA_tumor, 1),
    pct_NA_ctrl = replace_na(pct_NA_ctrl, 1)
  ) |>
  arrange(short_name, pct_NA_tumor) |>
  select(short_name, id, subid, pct_NA_tumor, pct_NA_ctrl, iMean_ctrl, iMean_tumor)
#- 10.3.5: Subset to top subids of the variants
validation_check_files_subids <- validation_check_files |>
  filter(!is.na(top_frag_subid)) |>
  rename(subid = top_frag_subid) |>
  select(id, subid, short_name)
#- 10.3.6: Now filter the combined table to just those subids
IARC_combined_top <- validation_check_files_subids |>
  left_join(IARC_combined, by = c("id", "subid"))
#! Pentachlorophenol excluded as the mz4 used for the "top" was not in original algorithm list
#+ 10.4: Full Joiner (Quant)
#- 10.3.2: full_joiner_i
#!!!!!!!!!!!!!!!!!!!!!
full_joiner <- full_joiner_i |>
  select(sample_ID, variant, tumor_vs_ctrl, any_of(cadaver_iarc_keep))
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
