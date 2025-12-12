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
moca <- iv_wide |>
  filter(short_name == "MOCA")
iarc_tumor_rtx <- rtx(
  validation_list = moca,
  iterate_through = 6,
  rt_lookup = "sample",
  save_rds = TRUE,
  rds_save_folder = "moca",
  overwrite_rds = TRUE,
  save_compiled_rds = FALSE,
  use_parallel = FALSE,
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
#+ 8.4: Read validation plots, compile, adjust x ranges
#- 8.4.1: Derive a list of all unique plots to read in
variant_plot_list <- validation_check_files %>%
  pull(plot) %>%
  str_split(",\\s*") %>%
  unlist() %>%
  unique()
#- 8.4.2: Copy all selected plots to curated/original folder
{
  variant_rtx_dir <- file.path(config$paths$validation_plot_directory_onedrive, "variant_rtx")
  iarc_tumor_dir <- file.path(config$paths$validation_plot_directory_onedrive, "iarc_tumor_rtx")
  curated_original_dir <- file.path(config$paths$validation_plot_directory_onedrive, "curated", "original")
  dir.create(curated_original_dir, showWarnings = FALSE, recursive = TRUE)
  copied_count <- 0
  for (plot_name in variant_plot_list) {
    # Check variant_rtx first
    variant_file <- file.path(variant_rtx_dir, paste0(plot_name, ".rds"))
    if (file.exists(variant_file)) {
      file.copy(variant_file, file.path(curated_original_dir, paste0(plot_name, ".rds")), overwrite = TRUE)
      copied_count <- copied_count + 1
      next
    }
    # Check iarc_tumor_rtx
    iarc_file <- file.path(iarc_tumor_dir, paste0(plot_name, ".rds"))
    if (file.exists(iarc_file)) {
      file.copy(iarc_file, file.path(curated_original_dir, paste0(plot_name, ".rds")), overwrite = TRUE)
      copied_count <- copied_count + 1
    }
  }
  cat(sprintf("✓ Copied %d plots to curated/original/\n\n", copied_count))
}
#- 8.4.3: Separate into modify vs final batches
{
  modify_plots <- validation_check_files %>%
    filter(state != "final") %>%
    pull(plot) %>%
    str_split(",\\s*") %>%
    unlist() %>%
    unique()
  final_plots <- validation_check_files %>%
    filter(state == "final") %>%
    pull(plot) %>%
    str_split(",\\s*") %>%
    unlist() %>%
    unique()
}
#- 8.4.4: Read modify_curated batch
{
  modify_curated <- list()
  for (plot_name in modify_plots) {
    file_path <- file.path(curated_original_dir, paste0(plot_name, ".rds"))
    if (file.exists(file_path)) {
      modify_curated[[plot_name]] <- readRDS(file_path)
      cat(sprintf("  ✓ Loaded: %s\n", plot_name))
    }
  }
  # Save modify_curated as RDS
  modify_rds_path <- file.path(config$paths$validation_plot_directory_onedrive, "curated", "modify_curated.rds")
  saveRDS(modify_curated, modify_rds_path)
  cat(sprintf("\n✓ Saved modify_curated RDS with %d plots: %s\n", length(modify_curated), modify_rds_path))
} 
#- 8.4.5: Read final_curated batch
{
  cat(sprintf("\nReading %d final plots...\n", length(final_plots)))
  final_curated <- list()
  for (plot_name in final_plots) {
    file_path <- file.path(curated_original_dir, paste0(plot_name, ".rds"))
    if (file.exists(file_path)) {
      final_curated[[plot_name]] <- readRDS(file_path)
      cat(sprintf("  ✓ Loaded: %s\n", plot_name))
    }
  }
  # Save final_curated as RDS
  final_rds_path <- file.path(config$paths$validation_plot_directory_onedrive, "curated", "final_curated.rds")
  saveRDS(final_curated, final_rds_path)
  cat(sprintf("\n✓ Saved final_curated RDS with %d plots: %s\n", length(final_curated), final_rds_path))
}
#- 8.4.6: Skip entire section if YAML specifies
} else {
  cat("⏭️  Skipping validation step (config$run_validation_step = FALSE)\n")
}