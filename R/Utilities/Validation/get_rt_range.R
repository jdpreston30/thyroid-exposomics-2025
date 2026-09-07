get_rt_range <- function(id_val, file_name, rt_data = tumor_rt_long) {
  if (is.na(file_name)) return(NA_real_)
  # Look up all id_subid values that match this compound (e.g., CP3017_0, CP3017_1, etc.)
#! arrange() before first(): the grepl matches EVERY subid for the compound, each a separate fragment with its own retention time, so first() on an unordered result silently depended on upstream row order. Both matches carried the same rt when this was found, so it never bit, but any change to the PeakWalk row set could have swapped which fragment was used. vp.R Step 6.5 orders identically so the two agree.
  rt_val <- rt_data |>
    filter(grepl(paste0("^", id_val, "_"), id_subid), file == file_name) |>
    arrange(id_subid) |>
    pull(rt) |>
    first()
  if (is.na(rt_val) || length(rt_val) == 0) return(NA_real_)
  return(rt_val)
}

