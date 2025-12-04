#* 7: Manual Spectral Validation QC
#+ 7.1: Pull features that need validation
#- 7.1.1: Subset to relevant features (top 5 quant)
top_5_quant_cas <- MT_final_i |>
  filter(mode == "quantitative") |>
  arrange(p_value) |>
  slice_head(n = 5) |>
  pull(cas)
#- 7.1.2: Subset to relevant IARCs (group 1)
iarc_1_validate <- tumor_raw |>
  filter(name_sub_lib_id %in% IARC_tumors_ctrl_filtered) |>
  pull(cas)
#+ 7.2: Create tumor validation table
tumor_top5_iarc <- tumor_raw |>
  filter(cas %in% top_5_quant_cas | cas %in% iarc_1_validate) |>
  mutate(
    validate_for = case_when(
      cas %in% top_5_quant_cas & cas %in% iarc_1_validate ~ "quant and iarc",
      cas %in% top_5_quant_cas ~ "quant",
      cas %in% iarc_1_validate ~ "iarc",
      TRUE ~ NA_character_
    )
  ) |>
  select(cas, short_display_name, lib, mzMed, rtMin, rtMed, rtMax, tmz, trt, id, subid, id_subid, name_sub_lib_id, validate_for) |>
  bind_cols(
    tumor_raw |>
      filter(cas %in% top_5_quant_cas | cas %in% iarc_1_validate) |>
      select(matches("^F[0-9]+$")) |>
      apply(1, function(row) {
        if (all(is.na(row))) {
          tibble(check1 = NA_character_, check2 = NA_character_, check3 = NA_character_,
                 check4 = NA_character_, check5 = NA_character_)
        } else {
          top_indices <- order(row, decreasing = TRUE, na.last = TRUE)
          top_5 <- head(names(row)[top_indices][!is.na(row[top_indices])], 5)
          top_5 <- c(top_5, rep(NA_character_, 5 - length(top_5)))
          tibble(check1 = top_5[1], check2 = top_5[2], check3 = top_5[3],
                 check4 = top_5[4], check5 = top_5[5])
        }
      }) |>
      bind_rows()
  ) |>
  add_file_references(file_list) |>
  arrange(cas, id_subid) |>
  arrange(desc(validate_for))
#+ 7.3: Create controls validation table
cadaver_iarc <- IARC_controls_i |>
  filter(cas %in% iarc_1_validate) |>
  left_join(feature_metadata |> select(cas, Short_display_name), by = "cas") |>
  mutate(across(matches("^T[0-9]+$"), as.numeric)) |>
  mutate(validate_for = "iarc") |>
  select(cas, short_display_name = Short_display_name, lib, mzMed, rtMin, rtMed, rtMax, tmz, trt, id, subid, id_subid, name_sub_lib_id, validate_for) |>
  bind_cols(
    IARC_controls_i |>
      filter(cas %in% iarc_1_validate) |>
      mutate(across(matches("^T[0-9]+$"), as.numeric)) |>
      select(matches("^T[0-9]+$")) |>
      apply(1, function(row) {
        if (all(is.na(row))) {
          tibble(check1 = NA_character_, check2 = NA_character_, check3 = NA_character_,
                 check4 = NA_character_, check5 = NA_character_)
        } else {
          top_indices <- order(row, decreasing = TRUE, na.last = TRUE)
          top_5 <- head(names(row)[top_indices][!is.na(row[top_indices])], 5)
          top_5 <- c(top_5, rep(NA_character_, 5 - length(top_5)))
          tibble(check1 = top_5[1], check2 = top_5[2], check3 = top_5[3],
                 check4 = top_5[4], check5 = top_5[5])
        }
      }) |>
      bind_rows()
  ) |>
  add_file_references(file_list) |>
  arrange(cas, id_subid) |>
  arrange(desc(validate_for))
#+ 7.4: Create mz reference table
#- 7.4.1: Create base reference table 
mz_reference_table_i <- GC2_features |>
  left_join(
    feature_metadata |> select(cas, short_display_name = Short_display_name),
    by = "cas"
  ) |>
  filter(cas %in% top_5_quant_cas | cas %in% iarc_1_validate) |>
  select(cas, name, short_display_name, id, trt, subid, tmz)|>
  # Pivot to create mz0, mz1, mz2, etc columns
  pivot_wider(
    id_cols = c(cas, name, short_display_name, id, trt),
    names_from = subid,
    values_from = tmz,
    names_prefix = "mz",
    values_fill = NA_real_,
    names_sort = TRUE
  ) |>
  arrange(short_display_name, id)
#- 7.4.2: Extract RT ranges from tumor data
tumor_rt_ranges <- tumor_top5_iarc |>
  distinct(cas, .keep_all = TRUE) |>
  select(cas, rtMin, rtMax) |>
  mutate(
    rtMin = round(rtMin, 2),
    rtMax = round(rtMax, 2),
    tumor_rt_range = paste0("c(", rtMin, ", ", rtMax, ")")
  ) |>
  select(cas, tumor_rt_range)
#- 7.4.3: Extract RT ranges from cadaver data
cadaver_rt_ranges <- cadaver_iarc |>
  distinct(cas, .keep_all = TRUE) |>
  select(cas, rtMin, rtMax) |>
  mutate(
    rtMin = round(rtMin, 2),
    rtMax = round(rtMax, 2),
    cadaver_rt_range = paste0("c(", rtMin, ", ", rtMax, ")")
  ) |>
  select(cas, cadaver_rt_range)
#- 7.4.4: Join RT ranges to mz_reference_table
mz_reference_table <- mz_reference_table_i |>
  left_join(tumor_rt_ranges, by = "cas") |>
  left_join(cadaver_rt_ranges, by = "cas")
#+ 7.5: Store appropriate metadata
#- 7.5.1: Pull list of chosen IARC fragments
IARC_fragments <- IARC_tumors_ctrl_filtered
#- 7.5.2: Pull list of analyzed fragments (variant)
quant_5_chosen_fragments <- quant_qual_results |>
  filter(mode == "quantitative") |>
  arrange(p_value) |>
  slice_head(n = 5) |>
  pull(name_sub_lib_id)
#- 7.5.3: Parse IARC fragments
iarc_parsed <- tibble(name_sub_lib_id = IARC_fragments) |>
  mutate(
    id = sub(".*_(CP[0-9]+)$", "\\1", name_sub_lib_id),
    subid = as.numeric(sub(".*_([0-9]+)_BP.*", "\\1", name_sub_lib_id)),
    source = "IARC"
  ) |>
  select(id, subid, source)
#- 7.5.4: Parse quantitative fragments
quant_parsed <- tibble(name_sub_lib_id = quant_5_chosen_fragments) |>
  mutate(
    id = sub(".*_(CP[0-9]+)$", "\\1", name_sub_lib_id),
    subid = as.numeric(sub(".*_([0-9]+)_BP.*", "\\1", name_sub_lib_id)),
    source = "quant"
  ) |>
  select(id, subid, source)
#- 7.5.5: Combine and tag sources (handle overlaps)
parsed_fragments <- bind_rows(iarc_parsed, quant_parsed) |>
  group_by(id, subid) |>
  summarize(
    source = if_else(n() > 1, "both", first(source)),
    .groups = "drop"
  )
#- 7.5.6: Filter GC2_features to only selected fragments
selected_gc2_features <- GC2_features |>
  inner_join(parsed_fragments, by = c("id", "subid")) |>
  select(source, id, subid, name, tmz, trt, cas) |>
  arrange(name, subid) |>
  arrange(source)
#- 7.5.7: Create short versions of references to visualize 
CS <- cadaver_iarc |>
  select(cas, short_display_name, id, file1, file2, file3)
TS <- tumor_top5_iarc |>
  select(cas, short_display_name, id, file1, file2, file3)
