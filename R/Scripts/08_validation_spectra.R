#* 8: Manual Spectral Validation
#+ 8.1: Convert raw files to mzML
#- 8.1.1: Setup base paths
#- 8.1.2: Convert all .raw files to mzML format
file_inventory <- convert_raw_to_mzml(
  file_list = file_list,
  tumor_raw_dir = config$paths$tumor_raw_dir,
  cadaver_raw_dir = config$paths$cadaver_raw_dir
)
#! Function will skip if already completed
#+ 8.2: Manual Validation (on 4-fragment versions)
#- 8.2.1: IARC Tumor
rtx(
  validation_list = iv_wide,
  iterate_through = 6,
  output_dir = "Outputs/Validation",
  pdf_name = "iarc_tumor_rtx.pdf",
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "iarc_tumor_rtx",
  overwrite_rds = TRUE,
  use_parallel = TRUE,
  n_cores = 10
)
#- 8.2.2: IARC Cadaver
rtx(
  validation_list = ic_wide,
  study = "cadaver",
  iterate_through = 7,
  output_dir = "Outputs/Validation",
  pdf_name = "iarc_cadaver_rtx.pdf",
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "iarc_cadaver_rtx",
  overwrite_rds = TRUE,
  use_parallel = TRUE,
  n_cores = 10
)
#- 8.2.3: Variant Differences Chemicals
rtx(
  validation_list = vv_wide,
  iterate_through = 6,
  output_dir = "Outputs/Validation/",
  pdf_name = "variant_rtx.pdf",
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "variant_rtx",
  overwrite_rds = TRUE,
  use_parallel = TRUE,
  n_cores = 6
)
#+ 8.3: Manually copy over files, read in, adjust x ranges
#- 8.3.0: Read in manual validation results
validation_check_files <- validation_check |>
  filter(!state %in% c("failed", "not used")) |>
  mutate(rt_range = (rtu-rtl)/2) |>
  select(-c(to_do, note, rtl, rtu)) |>
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
#- 8.4.1: Set vector of plots that need adjustment
adjust_ids <- c("CP2382", "CP3007", "CP2486", "CP2212", "CP1090", "CP3113", "CP3193", "CP3182", "CP3148", "CP2365", "CP3174", "CP3017", "CP3021", "CP2487", "CP1016", "CP2107", "CP3066")