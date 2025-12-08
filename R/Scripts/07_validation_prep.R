#* 7: Manual Spectral Validation QC
#+ 7.1: Pull variant comparison features that need validation
#- 7.1.0: Pull the id_subid for all variant features from MT_final_i
fragements_variant_pull <- MT_final_i |> pull(id_subid)
#- 7.1.1: Subset to relevant features (all in MT_final_i)
variant_validate_ids <- MT_final_i |>
  pull(id)
#- 7.1.2: Determine order of figures
variant_validate_order <- MT_final_i |>
  mutate(order = row_number()) |>
  select(id, order)
#- 7.1.2: Pull RT ranges and mzs for validation features
variant_validate_metadata_long <- combined_peakwalk_tumor |>
  filter(id %in% variant_validate_ids) |>
  mutate(
    rt_lower = round(pmax(trtlower - 5/60, 0), 2),
    rt_upper = round(trtupper + 5/60, 2),
    id_subid = paste0(id, "_", subid),
    asterisk = if_else(id_subid %in% fragements_variant_pull, "Y", "N")
  ) |>
  select(id, subid, id_subid, asterisk, tmz, monoisotopic, rt_lower, rt_upper)
#- 7.1.3: Calculate expanded RT range for each ID (min rt_lower, max rt_upper across all mzs)
expanded_rt_ranges <- variant_validate_metadata_long |>
  group_by(id) |>
  summarize(
    rt_lower_expanded = min(rt_lower, na.rm = TRUE),
    rt_upper_expanded = max(rt_upper, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(tumor_rt_range = paste0("c(", rt_lower_expanded, ", ", rt_upper_expanded, ")"))
#- 7.1.4: Pivot to wide format with mz columns; add on standard file names
variant_validate_metadata_wide_i <- variant_validate_metadata_long |>
  select(id, subid, id_subid, asterisk, tmz, monoisotopic) |>
  # Get unique monoisotopic per ID (same for all mzs)
  group_by(id) |>
  mutate(monoisotopic = first(monoisotopic)) |>
  ungroup() |>
  # Pivot tmz and asterisk to wide
  pivot_wider(
    id_cols = c(id, monoisotopic),
    names_from = subid,
    values_from = c(tmz, asterisk),
    names_glue = "{.value}{subid}",
    values_fill = list(tmz = NA_real_, asterisk = "N")
  ) |>
  # Rename tmz columns to mz columns
  rename_with(~gsub("^tmz", "mz", .x), starts_with("tmz")) |>
  # Join the expanded RT range
  left_join(expanded_rt_ranges |> select(id, tumor_rt_range), by = "id") |>
  # Join the order
  left_join(variant_validate_order, by = "id") |>
  # Reorder columns: order, id, monoisotopic, tumor_rt_range, then all mz columns
  relocate(order, id, monoisotopic, tumor_rt_range) |>
  arrange(order) |>
  # Add standards column based on CP prefix
  mutate(
    bp_prefix = case_when(
      grepl("^CP1", id) ~ "BP1",
      grepl("^CP2", id) ~ "BP2",
      grepl("^CP3", id) ~ "BP3",
      TRUE ~ NA_character_
    )
  ) |>
  rowwise() |>
  mutate(
    standards = {
      if (!is.na(bp_prefix)) {
        standard_cols <- colnames(combined_peakwalk_tumor)[grepl(paste0("^", bp_prefix), colnames(combined_peakwalk_tumor))]
        paste(standard_cols, collapse = ", ")
      } else {
        NA_character_
      }
    }
  ) |>
  ungroup() |>
  select(-bp_prefix)
#- 7.1.5: Pull top 6 sample files for each id (separate columns)
top_sample_files_i <- combined_peakwalk_tumor |>
  filter(id_subid %in% fragements_variant_pull) |>
  # Group by id first to get one row per id
  group_by(id) |>
  slice(1) |>
  ungroup() |>
  select(id, starts_with("BL_")) |>
  rowwise() |>
  mutate(
    # Get all BL_ values for this row (excluding id)
    bl_values = list(c_across(-id)),
    # Get the column names
    bl_cols = list(names(combined_peakwalk_tumor)[grepl("^BL_", names(combined_peakwalk_tumor))]),
    # Create named vector and sort by intensity
    sorted_files = list({
      vals <- bl_values[[1]]
      names(vals) <- bl_cols[[1]]
      top_files <- names(sort(vals, decreasing = TRUE, na.last = TRUE))
      top_files[!is.na(top_files)][1:6]  # Get top 6, may have NAs if fewer than 6
    }),
    # Create individual file columns
    file1 = sorted_files[[1]][1],
    file2 = sorted_files[[1]][2],
    file3 = sorted_files[[1]][3],
    file4 = sorted_files[[1]][4],
    file5 = sorted_files[[1]][5],
    file6 = sorted_files[[1]][6]
  ) |>
  ungroup() |>
  select(id, file1, file2, file3, file4, file5, file6)
#- 7.1.6: Combine metadata and top sample files
vv_wide_ii <- variant_validate_metadata_wide_i |>
  left_join(top_sample_files_i, by = "id") |>
  left_join(MT_final_i |> select(id, short_name), by = "id") |>
  relocate(short_name, .after = id)
#- 7.1.7: Consolidate asterisk columns into single column
vv_wide_i <- vv_wide_ii |>
  rowwise() |>
  mutate(
    asterisk = {
      # Get all asterisk columns and their values
      asterisk_vals <- c_across(matches("^asterisk[0-9]+$"))
      asterisk_col_names <- names(vv_wide_ii)[grepl("^asterisk[0-9]+$", names(vv_wide_ii))]
      
      # Find which ones are "Y" and build mz list
      marked_mzs <- character(0)
      for (i in seq_along(asterisk_vals)) {
        if (!is.na(asterisk_vals[i]) && asterisk_vals[i] == "Y") {
          subid <- sub("asterisk", "", asterisk_col_names[i])
          marked_mzs <- c(marked_mzs, paste0("mz", subid))
        }
      }
      
      if (length(marked_mzs) > 0) {
        paste(marked_mzs, collapse = ", ")
      } else {
        NA_character_
      }
    }
  ) |>
  ungroup() |>
  # Remove individual asterisk columns
  select(-matches("^asterisk[0-9]+$")) |>
  distinct()
#- 7.1.8: Add file-specific RT range vectors
vv_wide <- vv_wide_i |>
  rowwise() |>
  mutate(
    # Parse standards (comma-separated, only 1-2 standards)
    standard_files_vec = list(strsplit(standards, ", ")[[1]]),
    # Sample RT range vectors (f1-f6) - use direct column references
    f1_rt = get_rt_range(id, file1),
    f2_rt = get_rt_range(id, file2),
    f3_rt = get_rt_range(id, file3),
    f4_rt = get_rt_range(id, file4),
    f5_rt = get_rt_range(id, file5),
    f6_rt = get_rt_range(id, file6),
    # Standard RT range vectors (s1-s2)
    s1_rt = if (length(standard_files_vec[[1]]) >= 1) get_rt_range(id, standard_files_vec[[1]][1]) else NA_character_,
    s2_rt = if (length(standard_files_vec[[1]]) >= 2) get_rt_range(id, standard_files_vec[[1]][2]) else NA_character_
  ) |>
  ungroup() |>
  select(-standard_files_vec)
#+ 7.2: Pull tumor/cadaver comparison features that need validation
vv_wide |> colnames()

# Test first row to see file columns and RT ranges
vv_wide |> filter(order <= 1) |> 
  select(id, file1, file2, file3, file4, file5, file6, f1_rt, f2_rt, f3_rt, f4_rt, f5_rt, f6_rt)

# Check that RT ranges are populated
vv_wide |> 
  summarize(
    f1_populated = sum(!is.na(f1_rt)),
    f2_populated = sum(!is.na(f2_rt)),
    f3_populated = sum(!is.na(f3_rt)),
    f4_populated = sum(!is.na(f4_rt)),
    f5_populated = sum(!is.na(f5_rt)),
    f6_populated = sum(!is.na(f6_rt))
  )

# Full structure check
vv_wide |> glimpse()
#!!!!!!!!!!!!!!!
# #- 7.1.1: Subset to relevant features (all in MT_final_i)
# top_5_quant_cas <- MT_final_i |>
#   filter(mode == "quantitative") |>
#   arrange(p_value) |>
#   slice_head(n = 5) |>
#   pull(cas)
# #- 7.1.2: Subset to relevant IARCs (group 1)
# iarc_1_validate <- tumor_raw |>
#   filter(name_sub_lib_id %in% IARC_tumors_ctrl_filtered) |>
#   pull(cas)
# #+ 7.2: Create tumor validation table
# tumor_top5_iarc <- tumor_raw |>
#   filter(cas %in% top_5_quant_cas | cas %in% iarc_1_validate) |>
#   mutate(
#     validate_for = case_when(
#       cas %in% top_5_quant_cas & cas %in% iarc_1_validate ~ "quant and iarc",
#       cas %in% top_5_quant_cas ~ "quant",
#       cas %in% iarc_1_validate ~ "iarc",
#       TRUE ~ NA_character_
#     )
#   ) |>
#   select(cas, short_display_name, lib, mzMed, rtMin, rtMed, rtMax, tmz, trt, id, subid, id_subid, name_sub_lib_id, validate_for) |>
#   bind_cols(
#     tumor_raw |>
#       filter(cas %in% top_5_quant_cas | cas %in% iarc_1_validate) |>
#       select(matches("^F[0-9]+$")) |>
#       apply(1, function(row) {
#         if (all(is.na(row))) {
#           tibble(check1 = NA_character_, check2 = NA_character_, check3 = NA_character_,
#                  check4 = NA_character_, check5 = NA_character_)
#         } else {
#           top_indices <- order(row, decreasing = TRUE, na.last = TRUE)
#           top_5 <- head(names(row)[top_indices][!is.na(row[top_indices])], 5)
#           top_5 <- c(top_5, rep(NA_character_, 5 - length(top_5)))
#           tibble(check1 = top_5[1], check2 = top_5[2], check3 = top_5[3],
#                  check4 = top_5[4], check5 = top_5[5])
#         }
#       }) |>
#       bind_rows()
#   ) |>
#   add_file_references(file_list) |>
#   arrange(cas, id_subid) |>
#   arrange(desc(validate_for))
# #+ 7.3: Create controls validation table
# cadaver_iarc <- IARC_controls_i |>
#   filter(cas %in% iarc_1_validate) |>
#   left_join(feature_metadata |> select(cas, Short_display_name), by = "cas") |>
#   mutate(across(matches("^T[0-9]+$"), as.numeric)) |>
#   mutate(validate_for = "iarc") |>
#   select(cas, short_display_name = Short_display_name, lib, mzMed, rtMin, rtMed, rtMax, tmz, trt, id, subid, id_subid, name_sub_lib_id, validate_for) |>
#   bind_cols(
#     IARC_controls_i |>
#       filter(cas %in% iarc_1_validate) |>
#       mutate(across(matches("^T[0-9]+$"), as.numeric)) |>
#       select(matches("^T[0-9]+$")) |>
#       apply(1, function(row) {
#         if (all(is.na(row))) {
#           tibble(check1 = NA_character_, check2 = NA_character_, check3 = NA_character_,
#                  check4 = NA_character_, check5 = NA_character_)
#         } else {
#           top_indices <- order(row, decreasing = TRUE, na.last = TRUE)
#           top_5 <- head(names(row)[top_indices][!is.na(row[top_indices])], 5)
#           top_5 <- c(top_5, rep(NA_character_, 5 - length(top_5)))
#           tibble(check1 = top_5[1], check2 = top_5[2], check3 = top_5[3],
#                  check4 = top_5[4], check5 = top_5[5])
#         }
#       }) |>
#       bind_rows()
#   ) |>
#   add_file_references(file_list) |>
#   arrange(cas, id_subid) |>
#   arrange(desc(validate_for))
# #+ 7.4: Create mz reference table
# #- 7.4.1: Create base reference table 
# mz_reference_table_i <- GC2_features |>
#   left_join(
#     feature_metadata |> select(cas, short_display_name = Short_display_name),
#     by = "cas"
#   ) |>
#   filter(cas %in% top_5_quant_cas | cas %in% iarc_1_validate) |>
#   select(cas, name, short_display_name, id, trt, subid, tmz)|>
#   # Pivot to create mz0, mz1, mz2, etc columns
#   pivot_wider(
#     id_cols = c(cas, name, short_display_name, id, trt),
#     names_from = subid,
#     values_from = tmz,
#     names_prefix = "mz",
#     values_fill = NA_real_,
#     names_sort = TRUE
#   ) |>
#   arrange(short_display_name, id)
# #- 7.4.2: Extract RT ranges from tumor data
# tumor_rt_ranges <- tumor_top5_iarc |>
#   distinct(cas, .keep_all = TRUE) |>
#   select(cas, rtMin, rtMax) |>
#   mutate(
#     rtMin = round(rtMin, 2),
#     rtMax = round(rtMax, 2),
#     tumor_rt_range = paste0("c(", rtMin, ", ", rtMax, ")")
#   ) |>
#   select(cas, tumor_rt_range)
# #- 7.4.3: Extract RT ranges from cadaver data
# cadaver_rt_ranges <- cadaver_iarc |>
#   distinct(cas, .keep_all = TRUE) |>
#   select(cas, rtMin, rtMax) |>
#   mutate(
#     rtMin = round(rtMin, 2),
#     rtMax = round(rtMax, 2),
#     cadaver_rt_range = paste0("c(", rtMin, ", ", rtMax, ")")
#   ) |>
#   select(cas, cadaver_rt_range)
# #- 7.4.4: Join RT ranges to mz_reference_table
# mz_reference_table <- mz_reference_table_i |>
#   left_join(tumor_rt_ranges, by = "cas") |>
#   left_join(cadaver_rt_ranges, by = "cas")
# #+ 7.5: Store appropriate metadata
# #- 7.5.1: Pull list of chosen IARC fragments
# IARC_fragments <- IARC_tumors_ctrl_filtered
# #- 7.5.2: Pull list of analyzed fragments (variant)
# quant_5_chosen_fragments <- quant_qual_results |>
#   filter(mode == "quantitative") |>
#   arrange(p_value) |>
#   slice_head(n = 5) |>
#   pull(name_sub_lib_id)
# #- 7.5.3: Parse IARC fragments
# iarc_parsed <- tibble(name_sub_lib_id = IARC_fragments) |>
#   mutate(
#     id = sub(".*_(CP[0-9]+)$", "\\1", name_sub_lib_id),
#     subid = as.numeric(sub(".*_([0-9]+)_BP.*", "\\1", name_sub_lib_id)),
#     source = "IARC"
#   ) |>
#   select(id, subid, source)
# #- 7.5.4: Parse quantitative fragments
# quant_parsed <- tibble(name_sub_lib_id = quant_5_chosen_fragments) |>
#   mutate(
#     id = sub(".*_(CP[0-9]+)$", "\\1", name_sub_lib_id),
#     subid = as.numeric(sub(".*_([0-9]+)_BP.*", "\\1", name_sub_lib_id)),
#     source = "quant"
#   ) |>
#   select(id, subid, source)
# #- 7.5.5: Combine and tag sources (handle overlaps)
# parsed_fragments <- bind_rows(iarc_parsed, quant_parsed) |>
#   group_by(id, subid) |>
#   summarize(
#     source = if_else(n() > 1, "both", first(source)),
#     .groups = "drop"
#   )
# #- 7.5.6: Filter GC2_features to only selected fragments
# selected_gc2_features <- GC2_features |>
#   inner_join(parsed_fragments, by = c("id", "subid")) |>
#   select(source, id, subid, name, tmz, trt, cas) |>
#   arrange(name, subid) |>
#   arrange(source)
# #- 7.5.7: Create short versions of references to visualize 
# CS <- cadaver_iarc |>
#   select(cas, short_display_name, id, file1, file2, file3)
# TS <- tumor_top5_iarc |>
  # select(cas, short_display_name, id, file1, file2, file3)

