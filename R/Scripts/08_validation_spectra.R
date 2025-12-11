validation_check <- read_xlsx(config$paths$variant_validation, sheet = "validation")
#! Delete when done
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
#- 8.2.1: Variant Differences Chemicals
rtx(
  validation_list = vv_wide,
  iterate_through = 6,
  output_dir = "Outputs/Validation/",
  pdf_name = "variant_rtx.pdf",
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "variant_rtx",
  overwrite_rds = FALSE,
  use_parallel = TRUE,
  n_cores = 8
)
#- 8.2.2: IARC Tumor
rtx(
  validation_list = iv_wide,
  iterate_through = 6,
  output_dir = "Outputs/Validation",
  pdf_name = "iarc_tumor_rtx.pdf",
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "iarc_tumor_rtx",
  overwrite_rds = FALSE,
  use_parallel = TRUE,
  n_cores = 8
)
#- 8.2.3: IARC Cadaver
rtx(
  validation_list = ic_wide,
  study = "cadaver",
  iterate_through = 7,
  output_dir = "Outputs/Validation",
  pdf_name = "iarc_cadaver_rtx.pdf",
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "iarc_cadaver_rtx",
  overwrite_rds = FALSE,
  use_parallel = TRUE,
  n_cores = 8
)
#+ 8.3: Finalize curated validation plots
#- 8.3.0: Read in manual validation results
# validation_check_curate <- validation_check |>
#   filter(state == "final") |>
#   mutate(rt_range = (rtu-rtl)/2) |>
#   select(-c(to_do, note, rtl, rtu)) |>
#   arrange(order)
# #- 8.3.1: Pull all those files to repo
# copy_raw_validation_plots(
#   validation_curated = validation_check_curate,
#   config = config,
#   output_dir = config$paths$validation_plots_raw
# )
# #- 8.3.2: Read all copied RDS files into a single object
# rds_dir <- config$paths$validation_plots_raw
# rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
# validation_plots <- list()
# for (rds_file in rds_files) {
#   # Extract plot_tag from filename
#   plot_tag <- tools::file_path_sans_ext(basename(rds_file))
#   # Read RDS file
#   plot_data <- readRDS(rds_file)
#   # Store in list with plot_tag as key
#   validation_plots[[plot_tag]] <- plot_data
#   cat(sprintf("  Loaded: %s\n", plot_tag))
# }
# #- 8.3.3: Adjust x-axis RT ranges for each plot
# source("R/Utilities/Validation/adjust_validation_plot_ranges.R")
# validation_plots_adjusted <- adjust_validation_plot_ranges(
#   validation_plots = validation_plots,
#   validation_curated = validation_check_curate
# )

#!!!!!!!!!!!!!
#+ 8.5: List of Validated Chemicals
#- 8.5.0: Read in manual validation results
#- 8.5.1: Construct list of validated chemicals
#- 8.5.2: Construct list of validation graphs to dipslay
#+ 8.6: List of Validated Chemicals (IARC)
#- 8.6.0: Read in manual validation results
#- 8.6.1: Construct list of validated chemicals
#- 8.6.2: Construct list of validation graphs to dipslay