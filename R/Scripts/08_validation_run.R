#* 8: Manual Spectral Validation
if (config$run_validation_step) {
#+ 8.1: Convert raw files to mzML
file_inventory <- convert_raw_to_mzml(
  file_list = file_list,
  tumor_raw_dir = config$paths$tumor_raw_dir,
  cadaver_raw_dir = config$paths$cadaver_raw_dir
)
#+ 8.2: Manual Validation Plots Creation
#- 8.2.1: IARC Tumor
iarc_tumor_rtx <- rtx(
  validation_list = iv_wide,
  iterate_through = 6,
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "iarc_tumor_rtx",
  overwrite_rds = TRUE,
  save_compiled_rds = TRUE,
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
  save_compiled_rds = TRUE,
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
  save_compiled_rds = TRUE,
  use_parallel = FALSE
)
#+ 8.3: Compile PDFs of all Validation Plots
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
} else {
  cat("⏭️  Skipping validation step (config$run_validation_step = FALSE)\n")
}