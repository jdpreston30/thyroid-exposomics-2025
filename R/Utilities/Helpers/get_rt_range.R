get_rt_range <- function(id_val, file_name, rt_data = tumor_rt_long) {
  if (is.na(file_name)) return(NA_real_)
  # Look up all id_subid values that match this compound (e.g., CP3017_0, CP3017_1, etc.)
  rt_val <- rt_data |>
    filter(grepl(paste0("^", id_val, "_"), id_subid), file == file_name) |>
    pull(rt) |>
    first()
  if (is.na(rt_val) || length(rt_val) == 0) return(NA_real_)
  return(rt_val)
}

