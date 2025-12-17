#* 7: Prep Manual Spectral Validation QC
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
#- 7.1.3: Build validation table using helper function
vv_wide_base <- build_validation_table(
  validate_ids = variant_validate_ids,
  validate_order = variant_validate_order,
  source_label = "diff_variant",
  short_name_join = MT_final_i
)
#- 7.1.4: Add asterisk marking for significant fragments (using actual m/z values)
vv_wide_i <- vv_wide_base |>
  # Create id_subid combinations for each mz column to check against fragements_variant_pull
  rowwise() |>
  mutate(
    # Find which mz columns correspond to significant fragments
    asterisk = {
      # Get all mz column names (mz0, mz1, etc.)
      mz_cols <- names(vv_wide_base)[grepl("^mz[0-9]+$", names(vv_wide_base))]
      marked_mz_values <- numeric(0)
      
      # Get current row as a list to access columns by name
      current_row <- pick(everything())
      
      for (mz_col in mz_cols) {
        # Extract subid from column name (mz0 -> 0, mz1 -> 1, etc.)
        subid <- sub("^mz", "", mz_col)
        # Build id_subid
        id_subid_check <- paste0(id, "_", subid)
        # Check if this fragment is significant
        if (id_subid_check %in% fragements_variant_pull) {
          # Get the actual m/z value from this column
          mz_value <- current_row[[mz_col]]
          if (!is.na(mz_value)) {
            marked_mz_values <- c(marked_mz_values, mz_value)
          }
        }
      }
      if (length(marked_mz_values) > 0) {
        paste(marked_mz_values, collapse = ", ")
      } else {
        NA_character_
      }
    }
  ) |>
  ungroup()
#+ 7.2: Pull IARC Group 1 compounds from tumors for validation
#- 7.2.0: Pull the id_subid for all IARC fragments (cas match)
fragments_iarc_pull <- combined_peakwalk_tumor |>
  filter(cas %in% iarc_1) |>
  pull(id_subid)
#- 7.2.1: Subset to IDs for IARC compounds
iarc_validate_ids <- combined_peakwalk_tumor |>
  filter(cas %in% iarc_1) |>
  distinct(id) |>
  pull(id)
#- 7.2.3: Prepare short_name join table
iarc_short_names <- combined_peakwalk_tumor |>
  filter(cas %in% iarc_1) |>
  left_join(feature_metadata |> select(cas, Short_display_name), by = "cas") |>
  distinct(id, Short_display_name) |>
  rename(short_name = Short_display_name)
#- 7.2.4: Build validation table using helper function (no asterisk marking needed)
iv_wide_i <- build_validation_table(
  validate_ids = iarc_validate_ids,
  source_label = "iarc_variant",
  short_name_join = iarc_short_names
) |>
  mutate(asterisk = NA_character_) |>
  # Filtered out all PCBs
  filter(!grepl("^PCB", short_name))
#+ 7.3: Pull IARC Group 1 compounds from cadavers for validation
#- 7.3.1: Build cadaver validation table (reuse iarc_validate_ids and iarc_short_names from 7.2)
ic_wide_i <- build_validation_table(
  validate_ids = iarc_validate_ids,
  source_label = "iarc_cadaver",
  short_name_join = iarc_short_names,
  peakwalk_data = combined_peakwalk_cadaver,
  rt_data = cadaver_rt_long
) |>
  mutate(asterisk = NA_character_) |>
  # Filtered out all PCBs
  filter(!grepl("^PCB", short_name))
#+ 7.4: Create expanded validation combined table
#- 7.4.1: Merge all wides and remove duplicates
validation_combined <- bind_rows(
  vv_wide_i |> select(id, short_name, monoisotopic, mz0, mz1, mz2, mz3),
  iv_wide_i |> select(id, short_name, monoisotopic, mz0, mz1, mz2, mz3),
  ic_wide_i |> select(id, short_name, monoisotopic, mz0, mz1, mz2, mz3)
) |>
  distinct()
#- 7.4.2: Get validation IDs
val_ids <- validation_combined |>
  pull(id)
#- 7.4.3: Read expanded library
expanded_lib_features <- expanded_lib |>
  select(id = ID, emz1:emz8) |>
  filter(id %in% val_ids)
#- 7.4.4: Convert validation_combined mz columns to long format
mz_long <- validation_combined |>
  select(id, mz0:mz3) |>
  pivot_longer(cols = mz0:mz3, names_to = "fragment", values_to = "mz") |>
  filter(!is.na(mz))
#- 7.4.5: Convert expanded library to long format
emz_long <- expanded_lib_features |>
  pivot_longer(cols = emz1:emz8, names_to = "fragment", values_to = "mz") |>
  filter(!is.na(mz))
#- 7.4.6: Merge: Keep original mz fragments, add unique emz fragments
merged_fragments <- bind_rows(
  # Original fragments (priority)
  mz_long,
  # Add unique emz fragments (those not already in mz_long for each id)
  emz_long |>
    anti_join(mz_long, by = c("id", "mz"))
) |>
  distinct(id, mz, .keep_all = TRUE) |>
  arrange(id, mz) |>
  # Add potential_duplicate flag (within 10 ppm mass error)
  group_by(id) |>
  mutate(
    potential_duplicate = {
      # For each m/z, check if any other m/z is within 10 ppm
      ppm_diffs <- outer(mz, mz, function(x, y) abs(x - y) / x * 1e6)
      # Set diagonal to Inf (don't compare with self)
      diag(ppm_diffs) <- Inf
      # Check if any value in each row is <= 10 ppm
      ifelse(apply(ppm_diffs, 1, function(row) any(row <= 10)), "Y", "N")
    }
  ) |>
  ungroup() |>
  # Filter out emz fragments when potential_duplicate is "Y", keep original mz fragments
  filter(!(potential_duplicate == "Y" & grepl("^emz", fragment)))
#- 7.4.7: Convert back to wide format with proper numbering
expanded_validation_i <- merged_fragments |>
  group_by(id) |>
  # Sort: original mz fragments first (by name), then emz fragments
  arrange(id, !grepl("^mz", fragment), fragment) |>
  mutate(
    is_original = grepl("^mz", fragment),
    max_mz_num = max(as.numeric(gsub("mz", "", fragment[grepl("^mz", fragment)])), na.rm = TRUE),
    emz_counter = cumsum(!is_original)
  ) |>
  mutate(
    fragment_num = if_else(
      is_original,
      fragment,  # Keep original mz0, mz1, mz2, mz3 names
      paste0("mz", max_mz_num + emz_counter)  # Number emz starting after last original mz
    )
  ) |>
  ungroup() |>
  select(id, fragment_num, mz) |>
  # Pivot to wide format (will create as many columns as needed)
  pivot_wider(names_from = fragment_num, values_from = mz) |>
  # Add metadata back
  left_join(validation_combined |> select(id, short_name, monoisotopic), by = "id") |>
  # Reorder columns: id, short_name, monoisotopic, then all mz columns
  select(id, short_name, monoisotopic, starts_with("mz"))
#- 7.4.8: Reorganize fragments to ensure mz0 is always monoisotopic
expanded_validation <- expanded_validation_i |>
  rowwise() |>
  mutate(
    # Get all mz values and check if monoisotopic matches any within 20 ppm
    mz_values = list(c_across(starts_with("mz"))),
    ppm_diffs = list(abs(mz_values - monoisotopic) / monoisotopic * 1e6),
    min_ppm = min(ppm_diffs, na.rm = TRUE),
    matches_fragment = min_ppm <= 20,
    match_position = if_else(matches_fragment, which.min(ppm_diffs), NA_integer_)
  ) |>
  ungroup() |>
  # Now reorganize based on match results
  mutate(
    # Get mz column names dynamically
    across(starts_with("mz"), ~., .names = "original_{.col}")
  ) |>
  rowwise() |>
  mutate(
    # Reorganize based on match scenario
    reorganized_mz = list({
      current_mz <- c_across(starts_with("original_mz"))
      
      if (matches_fragment && match_position == 1) {
        # Scenario 1: Monoisotopic already matches mz0, no action needed
        current_mz
      } else if (matches_fragment && match_position > 1) {
        # Scenario 2: Monoisotopic matches a different column (e.g., mz4)
        # Move matched value to position 1, shift others up, fill gap
        matched_value <- current_mz[match_position]
        # Remove the matched position and shift everything down to fill gap
        before_match <- current_mz[1:(match_position - 1)]
        after_match <- if (match_position < length(current_mz)) {
          current_mz[(match_position + 1):length(current_mz)]
        } else {
          numeric(0)
        }
        # New order: matched value first, then all others (gap filled)
        c(matched_value, before_match, after_match)
      } else {
        # Scenario 3: Monoisotopic doesn't match any fragment
        # Add monoisotopic as new mz0, shift everything up
        c(monoisotopic, current_mz)
      }
    })
  ) |>
  ungroup() |>
  # Extract reorganized values back into individual mz columns
  mutate(
    mz0 = map_dbl(reorganized_mz, ~if_else(length(.x) >= 1, .x[1], NA_real_)),
    mz1 = map_dbl(reorganized_mz, ~if_else(length(.x) >= 2, .x[2], NA_real_)),
    mz2 = map_dbl(reorganized_mz, ~if_else(length(.x) >= 3, .x[3], NA_real_)),
    mz3 = map_dbl(reorganized_mz, ~if_else(length(.x) >= 4, .x[4], NA_real_)),
    mz4 = map_dbl(reorganized_mz, ~if_else(length(.x) >= 5, .x[5], NA_real_)),
    mz5 = map_dbl(reorganized_mz, ~if_else(length(.x) >= 6, .x[6], NA_real_)),
    mz6 = map_dbl(reorganized_mz, ~if_else(length(.x) >= 7, .x[7], NA_real_)),
    mz7 = map_dbl(reorganized_mz, ~if_else(length(.x) >= 8, .x[8], NA_real_)),
    mz8 = map_dbl(reorganized_mz, ~if_else(length(.x) >= 9, .x[9], NA_real_)),
    mz9 = map_dbl(reorganized_mz, ~if_else(length(.x) >= 10, .x[10], NA_real_))
  ) |>
  # Clean up temporary columns, keeping id, short_name, monoisotopic
  select(id, short_name, monoisotopic, mz0:mz9, -starts_with("original_"), -mz_values, -ppm_diffs, -min_ppm, -matches_fragment, -match_position, -reorganized_mz) |>
  # Special case: Remove 105.0699 from CP3017 and shift fragments down
  #! Removed CP3017 specific fragment due to adding too much noise
  rowwise() |>
  mutate(
    # For CP3017, identify and remove the 105.0699 fragment
    fragments_to_keep = if (id == "CP3017") {
      list({
        all_mz <- c_across(starts_with("mz"))
        # Remove 105.0699 (within small tolerance)
        all_mz[is.na(all_mz) | abs(all_mz - 105.0699) > 0.0001]
      })
    } else {
      list(c_across(starts_with("mz")))
    }
  ) |>
  ungroup() |>
  # Reassign mz columns with filtered/shifted values
  mutate(
    mz0 = map_dbl(fragments_to_keep, ~if_else(length(.x) >= 1, .x[1], NA_real_)),
    mz1 = map_dbl(fragments_to_keep, ~if_else(length(.x) >= 2, .x[2], NA_real_)),
    mz2 = map_dbl(fragments_to_keep, ~if_else(length(.x) >= 3, .x[3], NA_real_)),
    mz3 = map_dbl(fragments_to_keep, ~if_else(length(.x) >= 4, .x[4], NA_real_)),
    mz4 = map_dbl(fragments_to_keep, ~if_else(length(.x) >= 5, .x[5], NA_real_)),
    mz5 = map_dbl(fragments_to_keep, ~if_else(length(.x) >= 6, .x[6], NA_real_)),
    mz6 = map_dbl(fragments_to_keep, ~if_else(length(.x) >= 7, .x[7], NA_real_)),
    mz7 = map_dbl(fragments_to_keep, ~if_else(length(.x) >= 8, .x[8], NA_real_)),
    mz8 = map_dbl(fragments_to_keep, ~if_else(length(.x) >= 9, .x[9], NA_real_)),
    mz9 = map_dbl(fragments_to_keep, ~if_else(length(.x) >= 10, .x[10], NA_real_))
  ) |>
  # Clean up temporary columns, keeping id, short_name, monoisotopic and all mz columns
  select(id, short_name, monoisotopic, mz0:mz9, -fragments_to_keep)
#- 7.4.9: Verification: Check that mz0 now matches monoisotopic within 20 ppm for all compounds
monoisotopic_verification <- expanded_validation |>
  mutate(
    mz0_ppm_diff = abs(mz0 - monoisotopic) / monoisotopic * 1e6,
    mz0_matches_monoisotopic = mz0_ppm_diff <= 20
  ) |>
  select(id, short_name, monoisotopic, mz0, mz0_ppm_diff, mz0_matches_monoisotopic)
#+ 7.5: Update validation tables with expanded fragments
#- 7.5.1: Update variant validation table
vv_wide <- vv_wide_i |>
  select(-starts_with("mz"), -asterisk) |> # Remove old mz columns and asterisk
  left_join(
    expanded_validation |> select(id, starts_with("mz")),
    by = "id"
  ) |>
  # Convert asterisk from m/z values back to column names
  rowwise() |>
  mutate(
    asterisk = if (!is.na(vv_wide_i$asterisk[vv_wide_i$id == id])) {
      # Get the original asterisk values (m/z values as string)
      original_asterisk <- vv_wide_i$asterisk[vv_wide_i$id == id]
      # Split into individual m/z values
      asterisk_mz_values <- as.numeric(strsplit(original_asterisk, ", ")[[1]])
      # Find which columns match these m/z values
      current_mz <- c_across(starts_with("mz"))
      matched_cols <- character(0)
      for (mz_val in asterisk_mz_values) {
        # Find which column position matches this m/z (within small tolerance)
        match_pos <- which(abs(current_mz - mz_val) < 0.0001)
        if (length(match_pos) > 0) {
          matched_cols <- c(matched_cols, paste0("mz", match_pos[1] - 1))
        }
      }
      if (length(matched_cols) > 0) {
        paste(matched_cols, collapse = ", ")
      } else {
        NA_character_
      }
    } else {
      NA_character_
    }
  ) |>
  ungroup()
#- 7.5.2: Update IARC tumor validation table
iv_wide <- iv_wide_i |>
  select(-starts_with("mz")) |> # Remove old mz columns
  left_join(
    expanded_validation |> select(id, starts_with("mz")),
    by = "id"
  )
#- 7.5.3: Update IARC cadaver validation table
ic_wide <- ic_wide_i |>
  select(-starts_with("mz")) |> # Remove old mz columns
  left_join(
    expanded_validation |> select(id, starts_with("mz")),
    by = "id"
  )
#+ 7.6: Create subsetted version based on validated IARC (post hoc in step 9) with tag-ordered files
#- 7.6.0: Get list of top frags 
top_frags_validated_iarc <- validation_check_files |>
  filter(!is.na(top_frag)) |>
  select(id, top_frag)
#- 7.6.1: Get list of validated IARC1 chemicals
validation_iarcs <- top_frags_validated_iarc |>
  pull(id)
#- 7.6.1: Build IARC tumor validation table with tag ordering (reuse iarc_short_names from 7.2.3)
iv_wide_iarc_validated_i <- build_validation_table(
  validate_ids = validation_iarcs,
  source_label = "iarc_tumor_validated",
  short_name_join = iarc_short_names,
  order_by = "tag",
  buffer = (10/60)
) |>
  mutate(asterisk = NA_character_) |>
  filter(!grepl("^PCB", short_name))
#- 7.6.2: Build IARC cadaver validation table with tag ordering (reuse iarc_short_names from 7.2.3, peakwalk/rt data from 7.3)
ic_wide_iarc_validated_i <- build_validation_table(
  validate_ids = validation_iarcs,
  source_label = "iarc_cadaver_validated",
  short_name_join = iarc_short_names,
  peakwalk_data = combined_peakwalk_cadaver,
  rt_data = cadaver_rt_long,
  order_by = "tag",
  buffer = (10/60)
) |>
  mutate(asterisk = NA_character_) |>
  filter(!grepl("^PCB", short_name))
#- 7.6.3: Update with expanded fragments (reuse expanded_validation from 7.4.8)
# IARC tumor
iv_wide_iarc_validated <- iv_wide_iarc_validated_i |>
  select(-starts_with("mz")) |>
  left_join(
    expanded_validation |> select(id, starts_with("mz")),
    by = "id"
  ) |>
  left_join(top_frags_validated_iarc, by = "id") |>
  select(short_name, top_frag, everything())
# IARC cadaver
ic_wide_iarc_validated <- ic_wide_iarc_validated_i |>
  select(-starts_with("mz")) |>
  left_join(
    expanded_validation |> select(id, starts_with("mz")),
    by = "id"
  ) |>
  left_join(top_frags_validated_iarc, by = "id") |>
  select(short_name, top_frag, everything())
