build_validation_table <- function(validate_ids, 
                                    validate_order = NULL, 
                                    source_label = "validation",
                                    short_name_join = NULL,
                                    peakwalk_data = combined_peakwalk_tumor,
                                    rt_data = tumor_rt_long) {
  #' Build Validation Table
  #' 
  #' Creates a validation table with RT ranges, mz values, top sample files, and file-specific RT ranges
  #' 
  #' @param validate_ids Character vector of IDs to include in validation
  #' @param validate_order Optional tibble with columns: id, order (for ordering the final table). If NULL, no ordering applied.
  #' @param source_label Character string to label the source column (default: "validation")
  #' @param short_name_join Optional tibble with columns: id, short_name (for adding short names)
  #' @param peakwalk_data Tibble with peakwalk data (default: combined_peakwalk_tumor, can use combined_peakwalk_cadaver)
  #' @param rt_data Tibble with RT data in long format (default: tumor_rt_long, can use cadaver_rt_long)
  #' 
  #' @return Tibble with validation metadata including file-specific RT ranges
  
  # Pull mzs for validation features
  validate_metadata_long <- peakwalk_data |>
    filter(id %in% validate_ids) |>
    mutate(id_subid = paste0(id, "_", subid)) |>
    select(id, subid, id_subid, tmz, monoisotopic)
  
  # Calculate expanded RT range for each ID from actual RT measurements
  expanded_rt_ranges <- rt_data |>
    mutate(id = sub("_[0-9]+$", "", id_subid)) |>
    filter(id %in% validate_ids) |>
    group_by(id) |>
    summarize(
      rt_min = min(rt, na.rm = TRUE),
      rt_max = max(rt, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      rt_lower_expanded = round(pmax(rt_min - 5/60, 0), 2),
      rt_upper_expanded = round(rt_max + 5/60, 2),
      compound_rt_range = paste0("c(", rt_lower_expanded, ", ", rt_upper_expanded, ")")
    ) |>
    select(id, compound_rt_range)
  
  # Pivot to wide format with mz columns
  validate_metadata_wide_i <- validate_metadata_long |>
    select(id, subid, id_subid, tmz, monoisotopic) |>
    group_by(id) |>
    mutate(monoisotopic = first(monoisotopic)) |>
    ungroup() |>
    pivot_wider(
      id_cols = c(id, monoisotopic),
      names_from = subid,
      values_from = tmz,
      names_prefix = "mz",
      values_fill = NA_real_
    ) |>
    left_join(expanded_rt_ranges |> select(id, compound_rt_range), by = "id")
  
  # Apply ordering if validate_order is provided
  if (!is.null(validate_order)) {
    validate_metadata_wide_i <- validate_metadata_wide_i |>
      left_join(validate_order, by = "id") |>
      relocate(order, id, monoisotopic, compound_rt_range) |>
      arrange(order)
  } else {
    validate_metadata_wide_i <- validate_metadata_wide_i |>
      relocate(id, monoisotopic, compound_rt_range)
  }
  
  validate_metadata_wide_i <- validate_metadata_wide_i |>
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
          standard_cols <- colnames(peakwalk_data)[grepl(paste0("^", bp_prefix), colnames(peakwalk_data))]
          paste(standard_cols, collapse = ", ")
        } else {
          NA_character_
        }
      }
    ) |>
    ungroup() |>
    select(-bp_prefix)
  
  # Pull top 6 sample files for each compound
  top_sample_files <- peakwalk_data |>
    filter(id %in% validate_ids) |>
    group_by(id) |>
    slice(1) |>
    ungroup() |>
    select(id, starts_with("BL_")) |>
    pivot_longer(cols = starts_with("BL_"), names_to = "file", values_to = "intensity") |>
    filter(!is.na(intensity)) |>
    group_by(id) |>
    arrange(id, desc(intensity)) |>
    slice(1:6) |>
    mutate(file_rank = row_number()) |>
    ungroup() |>
    pivot_wider(
      id_cols = id,
      names_from = file_rank,
      values_from = file,
      names_prefix = "file"
    )
  
  # Combine metadata and top sample files
  validate_wide_ii <- validate_metadata_wide_i |>
    left_join(top_sample_files, by = "id")
  
  # Add short names if provided
  if (!is.null(short_name_join)) {
    validate_wide_ii <- validate_wide_ii |>
      left_join(short_name_join |> select(id, short_name), by = "id") |>
      relocate(short_name, .after = id)
  }
  
  # Add file-specific RT values and RT range vectors
  buffer <- 10/60  # Default buffer for rt_window
  validate_wide <- validate_wide_ii |>
    rowwise() |>
    mutate(
      standard_files_vec = list(strsplit(standards, ", ")[[1]]),
      f1_rt = get_rt_range(id, file1, rt_data = rt_data),
      f2_rt = get_rt_range(id, file2, rt_data = rt_data),
      f3_rt = get_rt_range(id, file3, rt_data = rt_data),
      f4_rt = get_rt_range(id, file4, rt_data = rt_data),
      f5_rt = get_rt_range(id, file5, rt_data = rt_data),
      f6_rt = get_rt_range(id, file6, rt_data = rt_data),
      s1_rt = if (length(standard_files_vec[[1]]) >= 1) get_rt_range(id, standard_files_vec[[1]][1], rt_data = rt_data) else NA_real_,
      s2_rt = if (length(standard_files_vec[[1]]) >= 2) get_rt_range(id, standard_files_vec[[1]][2], rt_data = rt_data) else NA_real_,
      f1_rt_range = if (!is.na(f1_rt)) sprintf("c(%.2f, %.2f)", f1_rt - buffer, f1_rt + buffer) else NA_character_,
      f2_rt_range = if (!is.na(f2_rt)) sprintf("c(%.2f, %.2f)", f2_rt - buffer, f2_rt + buffer) else NA_character_,
      f3_rt_range = if (!is.na(f3_rt)) sprintf("c(%.2f, %.2f)", f3_rt - buffer, f3_rt + buffer) else NA_character_,
      f4_rt_range = if (!is.na(f4_rt)) sprintf("c(%.2f, %.2f)", f4_rt - buffer, f4_rt + buffer) else NA_character_,
      f5_rt_range = if (!is.na(f5_rt)) sprintf("c(%.2f, %.2f)", f5_rt - buffer, f5_rt + buffer) else NA_character_,
      f6_rt_range = if (!is.na(f6_rt)) sprintf("c(%.2f, %.2f)", f6_rt - buffer, f6_rt + buffer) else NA_character_,
      s1_rt_range = if (!is.na(s1_rt)) sprintf("c(%.2f, %.2f)", s1_rt - buffer, s1_rt + buffer) else NA_character_,
      s2_rt_range = if (!is.na(s2_rt)) sprintf("c(%.2f, %.2f)", s2_rt - buffer, s2_rt + buffer) else NA_character_
    ) |>
    ungroup() |>
    select(-standard_files_vec) |>
    mutate(source = source_label)
  
  return(validate_wide)
}
