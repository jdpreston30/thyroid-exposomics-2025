#* 8: Manual Spectral Validation
#+ 8.1: Convert raw files to mzML
#- 8.1.1: Setup base paths
#- 8.1.2: Convert all .raw files to mzML format
file_inventory <- convert_raw_to_mzml(
  file_list = file_list,
  tumor_raw_dir = config$paths$tumor_raw_dir,
  cadaver_raw_dir = config$paths$cadaver_raw_dir
) #! Function will skip if already completed
#+ 8.2: Manual Validation Plots (on 4-fragment versions)
#- 8.2.1: IARC Tumor
iarc_tumor_rtx <- rtx(
  validation_list = iv_wide,
  iterate_through = 6,
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "iarc_tumor_rtx",
  overwrite_rds = TRUE,
  use_parallel = TRUE,
  n_cores = 8
)
#- 8.2.2: IARC Cadaver
iarc_cadaver_rtx <- rtx(
  validation_list = ic_wide,
  study = "cadaver",
  iterate_through = 7,
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "iarc_cadaver_rtx",
  overwrite_rds = TRUE,
  use_parallel = TRUE,
  n_cores = 8
)
#- 8.2.3: Variant Differences Chemicals
variant_rtx <- rtx(
  validation_list = vv_wide,
  iterate_through = 6,
  output_dir = "Outputs/Validation/initial_compile/",
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "variant_rtx",
  overwrite_rds = TRUE,
  use_parallel = FALSE
)
#+ 8.3: Manually copy over files, read in, adjust x ranges
#- 8.3.1: IARC Tumor
compile_validation_pdf(
  compound_plots = iarc_tumor_rtx,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = "iarc_tumor_rtx.pdf",
  add_plot_tags = TRUE
)
#- 8.3.2: IARC Cadaver
compile_validation_pdf(
  compound_plots = iarc_cadaver_rtx,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = "iarc_cadaver_rtx.pdf",
  add_plot_tags = TRUE
)
#- 8.3.3: Variant Differences Chemicals
compile_validation_pdf(
  compound_plots = variant_rtx,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = "variant_rtx.pdf",
  add_plot_tags = TRUE
)
#+ 8.4: Manually copy over files, read in, adjust x ranges
#! Rerun 8.3 once run is done
#- 8.3.0: Read in manual validation results
validation_check_files <- validation_check |>
  filter(!state %in% c("failed", "not used")) |>
  mutate(rt_range = (rtu-rtl)/2) |>
  select(-c(modification, note, rtl, rtu)) |>
  arrange(order)
#- 8.3.1: Pull all those files to repo
copy_raw_validation_plots(
  validation_curated = validation_check_files,
  config = config,
  output_dir = config$paths$validation_plots_raw
)
#- 8.3.2: Read all copied RDS files into a single object
rds_dir <- config$paths$validation_plots_raw
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
validation_plots <- list()
for (rds_file in rds_files) {
  # Extract plot_tag from filename
  plot_tag <- tools::file_path_sans_ext(basename(rds_file))
  # Read RDS file with error handling
  tryCatch({
    plot_data <- readRDS(rds_file)
    # Store in list with plot_tag as key
    validation_plots[[plot_tag]] <- plot_data
    cat(sprintf("  Loaded: %s\n", plot_tag))
  }, error = function(e) {
    warning(sprintf("  ⚠️  Failed to load %s: %s\n  File may be corrupted. Skipping.", 
                    plot_tag, e$message))
  })
}
#- 8.3.3: Adjust x-axis RT ranges for each plot
validation_plots_adjusted <- adjust_validation_plot_ranges(
  validation_plots = validation_plots,
  validation_curated = validation_check_files
)
#+ 8.4: Manual Adjustment of Specific Plots
source("R/Utilities/Helpers/remove_standard.R")
source("R/Utilities/Helpers/write_small.R")
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
