#* 8: Manual Spectral Validation
if (config$analysis$run_validation_step) {
#+ 8.1: Convert Raw Files to mzML
file_inventory <- convert_raw_to_mzml(
  file_list = file_list,
  tumor_raw_dir = config$paths$tumor_raw_dir,
  cadaver_raw_dir = config$paths$cadaver_raw_dir
)
#+ 8.2: Manual Validation Plots Creation
#! 8.2.1 and 8.2.2 are skippable when their grobs are already current on OneDrive. vp() falls back to reading validation_plots/{iarc_tumor_rtx,iarc_cadaver_rtx} when an object is absent from .GlobalEnv (vp.R:50-83), and script 09 needs only 16 of those objects. Nothing downstream consumes the iarc_tumor_rtx / iarc_cadaver_rtx variables -- they are used solely by the compile_validation_pdf() calls in this block.
#! ONLY skip after confirming the OneDrive folders hold the CURRENT grobs. If they are stale, vp() loads the old ones silently and the figures come out with the previous chemical names -- no error, no warning.
if (isTRUE(config$analysis$rebuild_iarc_grobs)) {
#- 8.2.1: IARC tumor
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
#- 8.2.2: IARC cadaver
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
} else {
  cat("⏭️  Skipping 8.2.1/8.2.2 IARC grob rebuild (config$analysis$rebuild_iarc_grobs = FALSE);\n")
  cat("   script 09 will load those grobs from the OneDrive backup via vp().\n")
}
#- 8.2.3: Variant differences chemicals (Part 1)
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
#- 8.2.4: Variant differences chemicals (Part 2)
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
#- 8.2.5: Variant differences chemicals (Part 3)
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
#+ 8.3: Iterate Through All Validated IARC1 (Post-hoc per Step 9)
#! Diagnostic only, and OFF by default. Nothing here reaches the supplement: the figure_order sheet points exclusively at Outputs/Validation/revised/grobs/, script 09 never references a *_validated / *_check object, and its section 9.10 ("IARC 1 Top Fragments") is an empty stub. These are nine SERIAL single-compound runs at iterate_through = 120 with use_parallel = FALSE, so they cost hours while producing PDFs for eyeballing fragments. Flip run_iarc_fragment_check to true in the yaml only when you actually want that inspection.
if (isTRUE(config$analysis$run_iarc_fragment_check)) {
#- 8.3.1: IARC tumor - Pentachlorophenol
{
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
}
#- 8.3.2: IARC tumor - γ-BHC
{
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
}
#- 8.3.3: IARC tumor - 2-Naphthylamine
{
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
}
#- 8.3.4: IARC tumor - Phenacetin
{
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
}
#- 8.3.5: IARC tumor - 4-ABP
{
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
}
#- 8.3.6: IARC tumor - MOCA
{
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
}
#- 8.3.7: IARC tumor - o-Toluidine
{
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
}
#- 8.3.8: IARC tumor - 2-ABP
{
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
}
#- 8.3.9: IARC cadaver
{
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
}
} else {
  cat("⏭️  Skipping 8.3 IARC fragment checks (config$analysis$run_iarc_fragment_check = FALSE)\n")
}
#+ 8.4: Skip Entire Section If YAML Specifies
} else {
  cat("⏭️  Skipping validation step (config$analysis$run_validation_step = FALSE)\n")
}
