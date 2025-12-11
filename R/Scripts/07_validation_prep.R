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
#- 7.1.3: Build validation table using helper function; add asterisk
vv_wide_base <- build_validation_table(
  validate_ids = variant_validate_ids,
  validate_order = variant_validate_order,
  source_label = "diff_variant",
  short_name_join = MT_final_i
)
#- 7.1.4: Add asterisk marking for significant fragments
vv_wide_i <- vv_wide_base |>
  # Create id_subid combinations for each mz column to check against fragements_variant_pull
  rowwise() |>
  mutate(
    # Find which mz columns correspond to significant fragments
    asterisk = {
      # Get all mz column names (mz0, mz1, etc.)
      mz_cols <- names(vv_wide_base)[grepl("^mz[0-9]+$", names(vv_wide_base))]
      marked_mzs <- character(0)
      
      for (mz_col in mz_cols) {
        # Extract subid from column name (mz0 -> 0, mz1 -> 1, etc.)
        subid <- sub("^mz", "", mz_col)
        # Build id_subid
        id_subid_check <- paste0(id, "_", subid)
        # Check if this fragment is significant
        if (id_subid_check %in% fragements_variant_pull) {
          marked_mzs <- c(marked_mzs, mz_col)
        }
      }
      
      if (length(marked_mzs) > 0) {
        paste(marked_mzs, collapse = ", ")
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
expanded_validation <- merged_fragments |>
  group_by(id) |>
  # Arrange by original fragment priority (mz before emz) and then by mz value
  arrange(id, !grepl("^mz", fragment), mz) |>
  # Renumber fragments sequentially starting from 0
  mutate(fragment_num = paste0("mz", row_number() - 1)) |>
  ungroup() |>
  select(id, fragment_num, mz) |>
  # Pivot to wide format (will create as many columns as needed)
  pivot_wider(names_from = fragment_num, values_from = mz) |>
  # Add metadata back
  left_join(validation_combined |> select(id, short_name, monoisotopic), by = "id") |>
  # Reorder columns: id, short_name, monoisotopic, then all mz columns
  select(id, short_name, monoisotopic, starts_with("mz"))
#+ 7.5: Create final versions with expanded fragments
#+ 7.5: Update validation tables with expanded fragments
#- 7.5.1: Update variant validation table
vv_wide <- vv_wide_i |>
  select(-starts_with("mz")) |> # Remove old mz columns
  left_join(
    expanded_validation |> select(id, starts_with("mz")),
    by = "id"
  )
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
