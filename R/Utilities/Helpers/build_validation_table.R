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
  
  # Pull ALL sample files for each compound (not just top 6)
  all_sample_files <- peakwalk_data |>
    filter(id %in% validate_ids) |>
    group_by(id) |>
    slice(1) |>
    ungroup() |>
    select(id, starts_with("BL_")) |>
    pivot_longer(cols = starts_with("BL_"), names_to = "file", values_to = "intensity") |>
    filter(!is.na(intensity)) |>
    group_by(id) |>
    arrange(id, desc(intensity)) |>
    mutate(
      file_rank = row_number(),
      all_files = paste(file, collapse = ", ")
    ) |>
    ungroup()
  
  # Pivot to wide format with ALL file columns (file1, file2, ..., fileN)
  top_sample_files <- all_sample_files |>
    select(id, file_rank, file, all_files) |>
    pivot_wider(
      id_cols = c(id, all_files),
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
  
  # Add file-specific RT values and RT range vectors dynamically for ALL files
  buffer <- 10/60  # Default buffer for rt_window
  
  # Get all file column names (file1, file2, ..., fileN)
  file_cols <- names(validate_wide_ii)[grepl("^file[0-9]+$", names(validate_wide_ii))]
  
  # Add RT values and ranges dynamically
  validate_wide <- validate_wide_ii |>
    rowwise() |>
    mutate(standard_files_vec = list(strsplit(standards, ", ")[[1]]))
  
  # Add f*_rt columns dynamically
  for (file_col in file_cols) {
    file_num <- gsub("file", "", file_col)
    rt_col <- paste0("f", file_num, "_rt")
    rt_range_col <- paste0("f", file_num, "_rt_range")
    
    validate_wide <- validate_wide |>
      mutate(
        !!rt_col := get_rt_range(id, .data[[file_col]], rt_data = rt_data),
        !!rt_range_col := if_else(
          !is.na(.data[[rt_col]]),
          sprintf("c(%.2f, %.2f)", .data[[rt_col]] - buffer, .data[[rt_col]] + buffer),
          NA_character_
        )
      )
  }
  
  # Add standard RT columns (s1_rt, s2_rt, etc.)
  validate_wide <- validate_wide |>
    mutate(
      s1_rt = if (length(standard_files_vec[[1]]) >= 1) get_rt_range(id, standard_files_vec[[1]][1], rt_data = rt_data) else NA_real_,
      s2_rt = if (length(standard_files_vec[[1]]) >= 2) get_rt_range(id, standard_files_vec[[1]][2], rt_data = rt_data) else NA_real_,
      s1_rt_range = if_else(!is.na(s1_rt), sprintf("c(%.2f, %.2f)", s1_rt - buffer, s1_rt + buffer), NA_character_),
      s2_rt_range = if_else(!is.na(s2_rt), sprintf("c(%.2f, %.2f)", s2_rt - buffer, s2_rt + buffer), NA_character_)
    ) |>
    ungroup() |>
    select(-standard_files_vec) |>
    mutate(source = source_label)
  
  return(validate_wide)
}
