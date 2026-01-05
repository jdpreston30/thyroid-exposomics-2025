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
# Run rtx
iarc_tumor_rtx <- rtx(
  validation_list = iv_wide,
  iterate_through = 6,
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "iarc_tumor_rtx",
  overwrite_rds = TRUE,
  use_parallel = TRUE,
  n_cores = 9
)
# Create compiled PDF
compile_validation_pdf(
  compound_plots = iarc_tumor_rtx,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = "iarc_tumor_rtx.pdf",
  add_plot_tags = TRUE,
  external_subfolder = "iarc_tumor_rtx"
)
#- 8.2.2: IARC Cadaver
# Run rtx
iarc_cadaver_rtx <- rtx(
  validation_list = ic_wide,
  study = "cadaver",
  iterate_through = 7,
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "iarc_cadaver_rtx",
  overwrite_rds = TRUE,
  use_parallel = TRUE,
  n_cores = 9
)
# Create compiled PDF
compile_validation_pdf(
  compound_plots = iarc_cadaver_rtx,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = "iarc_cadaver_rtx.pdf",
  add_plot_tags = TRUE,
  external_subfolder = "iarc_cadaver_rtx"
)
#- 8.2.3: Variant Differences Chemicals (Part 1)
# Subset to part 1
vv_wide_pt1 <- vv_wide |>
  slice(1:20)
# Run rtx
variant_rtx_pt1 <- rtx(
  validation_list = vv_wide_pt1,
  iterate_through = 6,
  output_dir = "Outputs/Validation/initial_compile/",
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "variant_rtx",
  overwrite_rds = TRUE,
  use_parallel = TRUE,
  n_cores = 8
)
# Write compiled PDF
compile_validation_pdf(
  compound_plots = variant_rtx_pt1,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = "variant_rtx_pt1.pdf",
  add_plot_tags = TRUE,
  external_subfolder = "variant_rtx"
)
#- 8.2.4: Variant Differences Chemicals (Part 2)
# Subset to part 2
vv_wide_pt2 <- vv_wide |>
  slice(21:40)
# Run rtx
variant_rtx_pt2<- rtx(
  validation_list = vv_wide_pt2,
  iterate_through = 6,
  output_dir = "Outputs/Validation/initial_compile/",
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "variant_rtx",
  overwrite_rds = TRUE,
  use_parallel = TRUE,
  n_cores = 8
)
# Write compiled PDF
compile_validation_pdf(
  compound_plots = variant_rtx_pt2,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = "variant_rtx_pt2.pdf",
  add_plot_tags = TRUE,
  external_subfolder = "variant_rtx"
)
#- 8.2.5: Variant Differences Chemicals (Part 3)
# Subset to part 3
vv_wide_pt3 <- vv_wide |>
  slice(41:n())
# Run rtx
variant_rtx_pt3 <- rtx(
  validation_list = vv_wide_pt3,
  iterate_through = 6,
  output_dir = "Outputs/Validation/initial_compile/",
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "variant_rtx",
  overwrite_rds = TRUE,
  use_parallel = TRUE,
  n_cores = 8
)
# Write compiled PDF
compile_validation_pdf(
  compound_plots = variant_rtx_pt3,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = "variant_rtx_pt3.pdf",
  add_plot_tags = TRUE,
  external_subfolder = "variant_rtx"
)
#+ 8.4: Iterate through all validated IARC1 (Post-hoc per step 9)
#- 8.4.1: IARC Tumor - Pentachlorophenol
iarc_tumor_rtx_validated_pt1 <- rtx(
  validation_list = iv_wide_iarc_validated |> slice(1),
  iterate_through = 120,
  rt_lookup = "window",
  window = 10/60,
  save_rds = TRUE,
  rds_save_folder = "iarc_tumor_rtx_validated_check",
  overwrite_rds = TRUE,
  use_parallel = FALSE,
  run_standard = FALSE,
  fragment_pare = FALSE,
  force_plot = TRUE
)
compile_validation_pdf(
  compound_plots = iarc_tumor_rtx_validated_pt1,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = paste0("iarc_tumor_", gsub("[^A-Za-z0-9]", "_", iv_wide_iarc_validated$short_name[1]), ".pdf"),
  add_plot_tags = TRUE,
  external_subfolder = "iarc_tumor_rtx_validated_check"
)
#- 8.4.2: IARC Tumor - γ-BHC
iarc_tumor_rtx_validated_pt2 <- rtx(
  validation_list = iv_wide_iarc_validated |> slice(2),
  iterate_through = 120,
  rt_lookup = "window",
  window = 10/60,
  save_rds = TRUE,
  rds_save_folder = "iarc_tumor_rtx_validated_check",
  overwrite_rds = TRUE,
  use_parallel = FALSE,
  run_standard = FALSE,
  fragment_pare = FALSE,
  force_plot = TRUE
)
compile_validation_pdf(
  compound_plots = iarc_tumor_rtx_validated_pt2,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = paste0("iarc_tumor_", gsub("[^A-Za-z0-9]", "_", iv_wide_iarc_validated$short_name[2]), ".pdf"),
  add_plot_tags = TRUE,
  external_subfolder = "iarc_tumor_rtx_validated_check"
)
#- 8.4.3: IARC Tumor - 2-Naphthylamine
iarc_tumor_rtx_validated_pt3 <- rtx(
  validation_list = iv_wide_iarc_validated |> slice(3),
  iterate_through = 120,
  rt_lookup = "window",
  window = 10/60,
  save_rds = TRUE,
  rds_save_folder = "iarc_tumor_rtx_validated_check",
  overwrite_rds = TRUE,
  use_parallel = FALSE,
  run_standard = FALSE,
  fragment_pare = FALSE,
  force_plot = TRUE
)
compile_validation_pdf(
  compound_plots = iarc_tumor_rtx_validated_pt3,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = paste0("iarc_tumor_", gsub("[^A-Za-z0-9]", "_", iv_wide_iarc_validated$short_name[3]), ".pdf"),
  add_plot_tags = TRUE,
  external_subfolder = "iarc_tumor_rtx_validated_check"
)
#- 8.4.4: IARC Tumor - Phenacetin
iarc_tumor_rtx_validated_pt4 <- rtx(
  validation_list = iv_wide_iarc_validated |> slice(4),
  iterate_through = 120,
  rt_lookup = "window",
  window = 10/60,
  save_rds = TRUE,
  rds_save_folder = "iarc_tumor_rtx_validated_check",
  overwrite_rds = TRUE,
  use_parallel = FALSE,
  run_standard = FALSE,
  fragment_pare = FALSE,
  force_plot = TRUE
)
compile_validation_pdf(
  compound_plots = iarc_tumor_rtx_validated_pt4,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = paste0("iarc_tumor_", gsub("[^A-Za-z0-9]", "_", iv_wide_iarc_validated$short_name[4]), ".pdf"),
  add_plot_tags = TRUE,
  external_subfolder = "iarc_tumor_rtx_validated_check"
)
#- 8.4.5: IARC Tumor - 4-ABP
iarc_tumor_rtx_validated_pt5 <- rtx(
  validation_list = iv_wide_iarc_validated |> slice(5),
  iterate_through = 120,
  rt_lookup = "window",
  window = 10/60,
  save_rds = TRUE,
  rds_save_folder = "iarc_tumor_rtx_validated_check",
  overwrite_rds = TRUE,
  use_parallel = FALSE,
  run_standard = FALSE,
  fragment_pare = FALSE,
  force_plot = TRUE
)
compile_validation_pdf(
  compound_plots = iarc_tumor_rtx_validated_pt5,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = paste0("iarc_tumor_", gsub("[^A-Za-z0-9]", "_", iv_wide_iarc_validated$short_name[5]), ".pdf"),
  add_plot_tags = TRUE,
  external_subfolder = "iarc_tumor_rtx_validated_check"
)
#- 8.4.6: IARC Tumor - MOCA
iarc_tumor_rtx_validated_pt6 <- rtx(
  validation_list = iv_wide_iarc_validated |> slice(6),
  iterate_through = 120,
  rt_lookup = "window",
  window = 10/60,
  save_rds = TRUE,
  rds_save_folder = "iarc_tumor_rtx_validated_check",
  overwrite_rds = TRUE,
  use_parallel = FALSE,
  run_standard = FALSE,
  fragment_pare = FALSE,
  force_plot = TRUE
)
compile_validation_pdf(
  compound_plots = iarc_tumor_rtx_validated_pt6,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = paste0("iarc_tumor_", gsub("[^A-Za-z0-9]", "_", iv_wide_iarc_validated$short_name[6]), ".pdf"),
  add_plot_tags = TRUE,
  external_subfolder = "iarc_tumor_rtx_validated_check"
)
#- 8.4.7: IARC Tumor - o-Toluidine
iarc_tumor_rtx_validated_pt7 <- rtx(
  validation_list = iv_wide_iarc_validated |> slice(7),
  iterate_through = 120,
  rt_lookup = "window",
  window = 10/60,
  save_rds = TRUE,
  rds_save_folder = "iarc_tumor_rtx_validated_check",
  overwrite_rds = TRUE,
  use_parallel = FALSE,
  run_standard = FALSE,
  fragment_pare = FALSE,
  force_plot = TRUE
)
compile_validation_pdf(
  compound_plots = iarc_tumor_rtx_validated_pt7,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = paste0("iarc_tumor_", gsub("[^A-Za-z0-9]", "_", iv_wide_iarc_validated$short_name[7]), ".pdf"),
  add_plot_tags = TRUE,
  external_subfolder = "iarc_tumor_rtx_validated_check"
)
#- 8.4.8: IARC Tumor - 2-ABP
iarc_tumor_rtx_validated_pt8 <- rtx(
  validation_list = iv_wide_iarc_validated |> slice(8),
  iterate_through = 120,
  rt_lookup = "window",
  window = 10/60,
  save_rds = TRUE,
  rds_save_folder = "iarc_tumor_rtx_validated_check",
  overwrite_rds = TRUE,
  use_parallel = FALSE,
  run_standard = FALSE,
  fragment_pare = FALSE,
  force_plot = TRUE
)
compile_validation_pdf(
  compound_plots = iarc_tumor_rtx_validated_pt8,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = paste0("iarc_tumor_", gsub("[^A-Za-z0-9]", "_", iv_wide_iarc_validated$short_name[8]), ".pdf"),
  add_plot_tags = TRUE,
  external_subfolder = "iarc_tumor_rtx_validated_check"
)
#- 8.4.9: Combine all IARC Tumor PDFs
cat("\n📚 Combining all IARC Tumor PDFs into single document...\n")
iarc_tumor_pdf_files <- list.files(
  path = "Outputs/Validation/initial_compile/",
  pattern = "^iarc_tumor_.*\\.pdf$",
  full.names = TRUE
)
pdftools::pdf_combine(
  input = iarc_tumor_pdf_files,
  output = "Outputs/Validation/initial_compile/iarc_tumor_rtx_validated_COMBINED.pdf"
)
cat(sprintf("✓ Combined %d PDFs into iarc_tumor_rtx_validated_COMBINED.pdf\n", length(iarc_tumor_pdf_files)))
#- 8.4.10: IARC Cadaver
iarc_cadaver_rtx_validated <- rtx(
  validation_list = ic_wide_iarc_validated,
  study = "cadaver",
  iterate_through = 16,  # 16
  rt_lookup = "window",
  window = 10/60,
  save_rds = TRUE,
  rds_save_folder = "iarc_cadaver_rtx_validated_check",
  overwrite_rds = TRUE,
  use_parallel = FALSE,
  run_standard = FALSE,
  fragment_pare = FALSE, 
  force_plot = TRUE,
  debug = TRUE
)
# Generate compiled PDF
compile_validation_pdf(
  compound_plots = iarc_cadaver_rtx_validated,
  output_dir = "Outputs/Validation/initial_compile/",
  pdf_name = "iarc_cadaver_rtx_validated.pdf",
  add_plot_tags = TRUE,
  external_subfolder = "iarc_cadaver_rtx_validated_check"
)
#+ 8.4: Skip entire section if YAML specifies
} else {
  cat("⏭️  Skipping validation step (config$analysis$run_validation_step = FALSE)\n")
}