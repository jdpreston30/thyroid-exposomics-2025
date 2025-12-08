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
#+ 8.2: Manual Validation (Variant Differences)
#- 8.2.1: Run combined RTX + MZX
variant_validation <- rtx_mzx(
  validation_list = vv_wide,
  study = "tumor",
  iterate_through = 6,
  output_dir = "Outputs/Validation",
  pdf_name = "tumor_validation.pdf",
  ppm_tolerance = 5,
  rt_tolerance = 10/60,
  rt_lookup = "sample",
  display_mode = "max",
  stick = FALSE,
  max_i = FALSE,
  plot_width = 3.9,
  plot_height = 3.25
)
#+ 8.3: IARC Chemicals Validation (Tumor)
#- 8.3.1: Run combined RTX + MZX
iarc_tumor_validation <- rtx_mzx(
  validation_list = iv_wide,
  study = "tumor",
  iterate_through = 6,
  output_dir = "Outputs/Validation",
  pdf_name = "iarc_tumor_validation.pdf",
  ppm_tolerance = 5,
  rt_tolerance = 10/60,
  rt_lookup = "sample",
  display_mode = "max",
  stick = FALSE,
  max_i = FALSE,
  plot_width = 3.9,
  plot_height = 3.25
)
#+ 8.4: IARC Chemicals Validation (Cadaver)
#- 8.4.1: Run combined RTX + MZX
iarc_cadaver_validation <- rtx_mzx(
  validation_list = ic_wide,
  study = "cadaver",
  iterate_through = 6,
  output_dir = "Outputs/Validation",
  pdf_name = "iarc_cadaver_validation.pdf",
  ppm_tolerance = 5,
  rt_tolerance = 10/60,
  rt_lookup = "sample",
  display_mode = "max",
  stick = FALSE,
  max_i = FALSE,
  plot_width = 3.9,
  plot_height = 3.25
)
#+ 8.5: List of Validated Chemicals
#- 8.5.0: Read in manual validation results
variant_diff_valid <- read_xlsx("metadata_files/manual_validation.csv", sheet = "variant_diff")
#- 8.5.1: Construct list of validated chemicals
#- 8.5.2: Construct list of validation graphs to dipslay
#+ 8.6: List of Validated Chemicals (IARC)
#- 8.6.0: Read in manual validation results
iarc_valid <- read_xlsx("metadata_files/manual_validation.csv", sheet = "iarc")
#- 8.6.1: Construct list of validated chemicals
#- 8.6.2: Construct list of validation graphs to dipslay