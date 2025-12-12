#* 9: Validation Plots Adjustment and Manual Review
#+ 9.1: Load compiled validation plots
{
  modify_curated_path <- file.path(config$paths$validation_plot_directory_onedrive, "curated", "modify_curated.rds")
  final_curated_path <- file.path(config$paths$validation_plot_directory_onedrive, "curated", "final_curated.rds")
  modify_curated <- readRDS(modify_curated_path)
  final_curated <- readRDS(final_curated_path)
  validation_plots <- c(modify_curated, final_curated)
  cat(sprintf("\n✓ Loaded %d plots from modify_curated and %d plots from final_curated\n", 
              length(modify_curated), length(final_curated)))
}
#+ 9.2: Global adjustments
#!!!!!!!!!!!!! 
validation_plots <- modify_curated
#!!!!!!!!!!!!!
#- 9.2.1: Adjust x-axis ranges based on all plots
VPA_colored <- adjust_VP_ranges(
  validation_plots = validation_plots,
  validation_curated = validation_check_files
)
#- 9.2.2: Adjust colors
VPA_colored <- adjust_VP_colors(
  validation_plots = validation_plots_adjusted_x
)
#- 9.2.3: Adjust y-axis titles
VPA <- adjust_VP_yaxis_title(VPA_colored, y_title = "← Standard | Tumor →")
#+ 9.3: Plots that failed after review
#- 9.3.1: Molinate (CP2486)
remove_standard(VPA$F1_S1_CP2486, xl = 5.1, xu = 5.175, subfolder = "failed")
#- 9.3.2: Atrazine (CP3113)
write_small(VPA$F6_S1_CP3113, subfolder = "failed")
#- 9.3.3: Bupirimate (CP2107)
write_small(VPA$F1_S1_CP2107, subfolder = "failed")
#- 9.3.4: Resmethrin (CP3174)
write_small(VPA$F2_S1_CP3174, subfolder = "failed")
#+ 9.6: Revise Plots (X-axis adjust and/or fragment zoom)
#- 9.6.1: MEHP (CP2382)
# quant frag tiny but there, zooming x axis, calling level 1
F2_S1_CP2382_R <- zoom_x(VPA$F2_S1_CP2382, xl = 9.75, xu = 10.2)
#!!!!!!!!!!!!! ZOOM PLOT?
#- 9.6.2: TEEP (CP3182)
# quant frag tiny but there, zooming x axis, calling level 1
F6_S1_CP3182_R <- zoom_x(VPA$F6_S1_CP3182, xl = 19.36, xu = 19.66)
#!!!!!!!!!!!!! ZOOM PLOT?
#- 9.6.3: N-MeFOSAA (CP3193)
# Removed mz3 to show smaller fragments
F3_S1_CP3193_R <- zoom_fragment(VPA$F3_S1_CP3193, -3)
#!!!!!!!!!!1 mz3 zoom plot?
#!!!! show original too though 
#+ 9.7: Revise Plots (Remove standard peak -> level 2)
#- 9.7.1: MDA (CP3007)
# removed standard, zoomed Y axis, calling level 2
F3_S1_CP3007_R <- remove_standard(VPA$F3_S1_CP3007, xl = 3.8, xu = 4)
#- 9.4.8: Menthone (CP3148)
F3_S1_CP3007_R <- remove_standard(VPA$F3_S1_CP3007, xl = 3.8, xu = 4)


#- 9.4.9 : Prosulfuron (CP2365)
#- 9.4.10: Resmethrin (CP3174)
#- 9.4.11: o-Toluidine (CP3017)
#- 9.4.12: o-Anisidine (CP3021)
#- 9.4.13: Vernolate (CP2487)
#- 9.4.16: Bupirimate (CP2107)
#- 9.4.17: o-Cresol (CP3066)
#+ 9.8: Revise Plots (Other)
#- 9.8.1: o-Toluidine (CP3017)
#! Already modified in earlier steps with removal of bad fragment

#+ 9.5 Save 
#!!! SAVE PLOTS ONLY LOCALLY TO REPO