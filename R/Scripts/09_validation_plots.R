#* 9: Validation Plots Adjustment and Manual Review
#+ 9.1: Load compiled validation plots
validation_plots <- load_validation_plots(
  onedrive_base_path = config$paths$validation_plot_directory_onedrive
)
#+ 9.2: Adjust x-axis RT ranges for each plot
validation_plots_adjusted <- adjust_validation_plot_ranges(
  validation_plots = validation_plots,
  validation_curated = validation_check_files
)
#+ 9.4: Manual Adjustment of Specific Plots
#- 8.4.1: Set vector of plots that need adjustment
adjust_ids <- c("CP2382", "CP3007", "CP2486", "CP2212", "CP1090", "CP3113", "CP3193", "CP3182", "CP3148", "CP2365", "CP3174", "CP3017", "CP3021", "CP2487", "CP1016", "CP2107", "CP3066")
#- 8.3.4: Pull all adjust plots for each ID into separate lists
aps <- list()
for (id in adjust_ids) {
  id_plots <- validation_plots_adjusted[grepl(paste0("_", id, "$"), names(validation_plots_adjusted))]
  if (length(id_plots) > 0) {
    aps[[id]] <- id_plots
    cat(sprintf("Pulled %d plots for %s\n", length(id_plots), id))
  }
}
#+ 8.4: Plots that failed after review
#- 8.4.1: Molinate (CP2486)
F3_S1_CP2486_revised <- remove_standard(aps$CP2486$F1_S1_CP2486, xl = 5.1, xu = 5.175)
#- 8.4.2: Atrazine (CP3113)
write_small(aps$CP3113)
#- 8.4.16: Bupirimate (CP2107)
write_small(aps$CP2107)
#- 8.4.10: Resmethrin (CP3174)
write_small(aps$CP3174)
#+ 8.6: Plots that required revision
#- 8.4.1: MDA (CP3007)
#! removed standard, zoomed Y axis, calling level 2
F3_S1_CP3007_revised <- remove_standard(aps$CP3007$F3_S1_CP3007, xl = 3.8, xu = 4)
#- 8.4.5: N-MeFOSAA (CP3193)
#- 8.4.3: MEHP
#- 8.4.7: TEEP (CP3182)
#- 8.4.8: Menthone (CP3148)
#- 8.4.8: Menthone (CP3148)
#- 8.4.9 : Prosulfuron (CP2365)
#- 8.4.10: Resmethrin (CP3174)
#- 8.4.11: o-Toluidine (CP3017)
#- 8.4.12: o-Anisidine (CP3021)
#- 8.4.13: Vernolate (CP2487)
#- 8.4.16: Bupirimate (CP2107)
#- 8.4.17: o-Cresol (CP3066)
#+ 8.5 Save 
