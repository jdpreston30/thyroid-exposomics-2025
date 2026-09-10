#* rtx plot parts -- the pieces of a validation plot that must be built OUTSIDE process_single_compound()
#! aes() quosures and the scale breaks/labels closures capture the environment they are created in. Inside process_single_compound() that is the worker's frame: the chromatogram tables, the mzML slices and every earlier plot of the compound. Under ggplot2 4.0 each saved grob therefore serialized to 60-170 MB (5-7 MB gzipped), the disk-read compile churned through ~1.5 GB per compound and 9 GB went to swap on 2026-09-10. Built here, the captured frame is these arguments only. Adding labs()/theme()/coord_cartesian() to the returned plot back in the big frame is fine: those carry values, not closures (verified: 16.5 MB helper-built vs 16.8 MB after + labs + theme in a 40 MB frame). The ~16 MB floor per plot is ggplot2 4.0's S7 class overhead, not ours.
#+ Sample-only chromatogram (rtx mode 1)
rtx_sample_plot <- function(sample_chrom, stick) {
  p <- ggplot(sample_chrom, aes(x = rt, y = intensity, color = mz_label))
  if (stick) p + geom_segment(aes(xend = rt, yend = 0), linewidth = 0.4) else p + geom_line(linewidth = 0.4)
}
rtx_sample_nodata_plot <- function(sample_rt_range, y_limit) {
  ggplot() +
    annotate("text", x = mean(sample_rt_range), y = y_limit/2,
             label = "NO DATA", size = 6, color = "gray50", fontface = "bold") +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, y_limit),
      n.breaks = 8,
      labels = scales::label_scientific(digits = 2)
    )
}
rtx_sample_scale_x <- function(use_hard_limits, sample_rt_range) {
  brk <- function(limits) {
    start <- ceiling(limits[1] * 20) / 20
    end <- floor(limits[2] * 20) / 20
    if (start < end) seq(start, end, by = 0.05) else seq(start, end, by = -0.05)
  }
  mbrk <- function(limits) {
    start <- ceiling(limits[1] * 40) / 40
    end <- floor(limits[2] * 40) / 40
    if (start < end) seq(start, end, by = 0.025) else seq(start, end, by = -0.025)
  }
  if (use_hard_limits) {
    scale_x_continuous(expand = c(0, 0), limits = sample_rt_range, breaks = brk, minor_breaks = mbrk)
  } else {
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05), add = 0), breaks = brk, minor_breaks = mbrk)
  }
}
#+ Sample-over-standard mirror chromatogram (rtx mode 2)
rtx_mirror_plot <- function(combined_data, stick) {
  p <- ggplot(combined_data, aes(x = rt, y = plot_intensity, color = mz_label, group = interaction(mz_label, type)))
  p <- if (stick) p + geom_segment(aes(xend = rt, yend = 0), linewidth = 0.4) else p + geom_line(linewidth = 0.4)
  p + geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4)
}
rtx_mirror_y_labels <- function(x) scales::label_scientific(digits = 2)(abs(x))
rtx_mirror_scale_y <- function(y_limit) {
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(-y_limit, y_limit),
    labels = rtx_mirror_y_labels,
    n.breaks = 8
  )
}
rtx_mirror_nodata_plot <- function(sample_rt_range, y_limit) {
  ggplot() +
    annotate("text", x = mean(sample_rt_range), y = 0,
             label = "NO DATA", size = 6, color = "gray50", fontface = "bold") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
    rtx_mirror_scale_y(y_limit)
}
rtx_mirror_scale_x <- function(use_hard_limits, sample_rt_range) {
  brk  <- function(limits) seq(ceiling(limits[1] * 20) / 20, floor(limits[2] * 20) / 20, by = 0.05)
  mbrk <- function(limits) seq(ceiling(limits[1] * 40) / 40, floor(limits[2] * 40) / 40, by = 0.025)
  if (use_hard_limits) {
    scale_x_continuous(expand = c(0, 0), limits = sample_rt_range, breaks = brk, minor_breaks = mbrk)
  } else {
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05), add = 0), breaks = brk, minor_breaks = mbrk)
  }
}
