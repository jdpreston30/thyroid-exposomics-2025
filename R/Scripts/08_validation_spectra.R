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
variant_validation <- rtx(
  validation_list = vv_wide,
  study = "tumor",
  iterate_through = 6,
  output_dir = "Outputs/Validation",
  pdf_name = "tumor_spectra.pdf",
  ppm_tolerance = 5,
  stick = FALSE,
  max_i = FALSE,
  plot_width = 3.9,
  plot_height = 3.25,
  rt_lookup = "sample"
)
#+ 8.3: IARC chemicals validation
#- 8.3.1: Run for tumors
iarc_tumor_validation <- rtx(
  validation_list = iv_wide,
  study = "tumor",
  iterate_through = 6,
  output_dir = "Outputs/Validation",
  pdf_name = "iarc_tumor_spectra.pdf",
  ppm_tolerance = 5,
  stick = FALSE,
  max_i = FALSE,
  plot_width = 3.9,
  plot_height = 3.25,
  rt_lookup = "sample"
)
#- 8.3.2: Run for cadaver
iarc_cadaver_validation <- rtx(
  validation_list = ic_wide,
  study = "cadaver",
  iterate_through = 6,
  output_dir = "Outputs/Validation",
  pdf_name = "iarc_cadaver_spectra.pdf",
  ppm_tolerance = 5,
  stick = FALSE,
  max_i = FALSE,
  plot_width = 3.9,
  plot_height = 3.25,
  rt_lookup = "sample"
)
#+ 8.4: List of Validated Chemicals (Variant)
#- 8.4.1: Construct list of validated chemicals

#- 8.4.2: Construct list of validation graphs to dipslay
#+ 8.5: List of Validated Chemicals (IARC)
#- 8.5.1: Construct list of validated chemicals
#- 8.5.2: Construct list of validation graphs to dipslay