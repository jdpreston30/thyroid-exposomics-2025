subtitle_rt_of <- function(plot_obj, fallback_range = NULL) {
  #' Resolve the retention time to print in a rebuilt subtitle
  #'
  #' @param plot_obj A validation plot object from process_single_compound().
  #' @param fallback_range Optional numeric RT range used only when nothing better exists.
  #'
  #' @return Numeric retention time in minutes.
#! Never recompute mean(rt_range) here. rt_range comes from build_validation_table.R:169, which writes its endpoints through sprintf("%.2f"), so averaging them loses up to 0.005 min and flips the printed digit on exact .xx5 ties (o-Toluidine at 8.085). Prefer subtitle_rt, which process_single_compound.R carries from the unrounded fN_rt column; then the RT already rendered into the subtitle; and only then the lossy midpoint, for grobs generated before subtitle_rt existed.
  if (!is.null(plot_obj$subtitle_rt) && !is.na(plot_obj$subtitle_rt)) {
    return(as.numeric(plot_obj$subtitle_rt))
  }
  st <- plot_obj$plot$labels$subtitle
  if (!is.null(st) && !is.na(st) && grepl("RT = [0-9.]+ min", st)) {
    return(as.numeric(sub(".*RT = ([0-9.]+) min.*", "\\1", st)))
  }
  rng <- if (!is.null(fallback_range)) fallback_range else plot_obj$rt_range
  if (!is.null(rng) && length(rng) >= 2) return(mean(rng))
  NA_real_
}
