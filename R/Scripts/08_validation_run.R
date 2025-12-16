#* 8: Manual Spectral Validation
if (config$analysis$run_validation_step) {
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
#- 8.2.3: Variant Differences Chemicals (Part 1)
vv_wide_pt1 <- vv_wide |>
  slice(1:20)
variant_rtx_pt1 <- rtx(
  validation_list = vv_wide_pt1,
  iterate_through = 6,
  output_dir = "Outputs/Validation/initial_compile/",
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "variant_rtx",
  overwrite_rds = TRUE,
  save_compiled_rds = TRUE,
  use_parallel = FALSE
)
#!!!!
variant_rtx_pt1 <- variant_rtx
#!!!
#- 8.2.4: Variant Differences Chemicals (Part 2)
vv_wide_pt2 <- vv_wide |>
  slice(21:40)
variant_rtx_pt2<- rtx(
  validation_list = vv_wide_pt2,
  iterate_through = 6,
  output_dir = "Outputs/Validation/initial_compile/",
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "variant_rtx",
  overwrite_rds = TRUE,
  save_compiled_rds = TRUE,
  use_parallel = FALSE
)
#- 8.2.5: Variant Differences Chemicals (Part 3)
vv_wide_pt3 <- vv_wide |>
  slice(41:n())
variant_rtx_pt3 <- rtx(
  validation_list = vv_wide_pt3,
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
#- 8.3.3: Variant Differences Chemicals (Part 1)
compile_validation_pdf(
  compound_plots = variant_rtx_pt1,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = "variant_rtx_pt1.pdf",
  add_plot_tags = TRUE
)
#- 8.3.4: Variant Differences Chemicals (Part 2)
compile_validation_pdf(
  compound_plots = variant_rtx_pt2,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = "variant_rtx_pt2.pdf",
  add_plot_tags = TRUE
)
#- 8.3.5: Variant Differences Chemicals (Part 3)
compile_validation_pdf(
  compound_plots = variant_rtx_pt3,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = "variant_rtx_pt3.pdf",
  add_plot_tags = TRUE
)
#+ 8.4: Skip entire section if YAML specifies
} else {
  cat("⏭️  Skipping validation step (config$analysis$run_validation_step = FALSE)\n")
}