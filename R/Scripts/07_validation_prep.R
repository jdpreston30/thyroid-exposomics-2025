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
vv_wide <- vv_wide_base |>
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
iv_wide <- build_validation_table(
  validate_ids = iarc_validate_ids,
  source_label = "iarc_variant",
  short_name_join = iarc_short_names
) |>
  mutate(asterisk = NA_character_) |>
  # Filtered out all PCBs
  filter(!grepl("^PCB", short_name))
#+ 7.3: Pull IARC Group 1 compounds from cadavers for validation
#- 7.3.1: Build cadaver validation table (reuse iarc_validate_ids and iarc_short_names from 7.2)
ic_wide <- build_validation_table(
  validate_ids = iarc_validate_ids,
  source_label = "iarc_cadaver",
  short_name_join = iarc_short_names,
  peakwalk_data = combined_peakwalk_cadaver,
  rt_data = cadaver_rt_long
) |>
  mutate(asterisk = NA_character_) |>
  # Filtered out all PCBs
  filter(!grepl("^PCB", short_name))
