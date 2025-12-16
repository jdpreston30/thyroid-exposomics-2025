#- Copy

{
CP <- "CP3021"
source <- "variant_rtx"
source_dir <- file.path(config$paths$validation_plot_directory_onedrive, source)
dest_dir <- file.path(config$paths$validation_plot_directory_onedrive, "curated", "original")
files <- list.files(source_dir, pattern = paste0(".*", CP, ".*\\.rds"), full.names = TRUE)
file.copy(files, dest_dir, overwrite = TRUE)
cat(sprintf("✓ Copied %d files matching %s\n", length(files), CP))
}
#- Check other files
#_ o-Tolidine (CP2475)
F1_S1_CP2475_R <- vp(F1_S1_CP2475, xl = 7.85, xu = 8,  yl = -1.3e6, yu = 2e6)
mz_fragment = -3)
F2_S1_CP2475_R <- vp(F2_S1_CP2475)
F3_S1_CP2475_R <- vp(F3_S1_CP2475)
F4_S1_CP2475_R <- vp(F4_S1_CP2475)
F5_S1_CP2475_R <- vp(F5_S1_CP2475)
F6_S1_CP2475_R <- vp(F6_S1_CP2475)
#_Methylparaben (CP2252)
F1_S1_CP2252_R <- vp(F1_S1_CP2252, xl = 4.3, xu = 4.50)

#+ check
#_2-Naphthylamine (CP2535)
F5_S1_CP2535_R <- vp(F5_S1_CP2535, xl = 7.4, xu = 7.6)
#_Dibutyl phthalate (CP3047)
F1_S1_CP3047_R <- vp(F1_S1_CP3047, xl = NULL, xu = NULL)
#_Napropamide (CP3135)
F1_S1_CP3135_R <- vp(F1_S1_CP3135, xl = NULL, xu = NULL)
#_OD-PABA (CP2331)
F4_S2_CP2331_R <- vp(F4_S2_CP2331, xl = NULL, xu = NULL)
#_Flucythrinate (CP3166)
F3_S1_CP3166_R <- vp(F3_S1_CP3166, xl = NULL, xu = NULL)


