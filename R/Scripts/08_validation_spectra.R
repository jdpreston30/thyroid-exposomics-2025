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
variant_rtx <- rtx(
  validation_list = vv_wide,
  iterate_through = 6,
  output_dir = "Outputs/Validation/",
  pdf_name = "variant_rtx.pdf",
  rt_lookup = "sample"
)
#- 8.2.2: IARC Tumor
iarc_tumor_rtx <- rtx(
  validation_list = iv_wide,
  iterate_through = 6,
  output_dir = "Outputs/Validation",
  pdf_name = "iarc_tumor_rtx.pdf",
  rt_lookup = "sample"
)
#- 8.2.3: IARC Cadaver
iarc_cadaver_rtx <- rtx(
  validation_list = ic_wide,
  study = "cadaver",
  iterate_through = 7,
  output_dir = "Outputs/Validation",
  pdf_name = "iarc_cadaver_rtx.pdf",
  rt_lookup = "sample"
)
#+ 8.3: 


#!!!!!!!!!!!!!
#+ 8.5: List of Validated Chemicals
#- 8.5.0: Read in manual validation results
#- 8.5.1: Construct list of validated chemicals
#- 8.5.2: Construct list of validation graphs to dipslay
#+ 8.6: List of Validated Chemicals (IARC)
#- 8.6.0: Read in manual validation results
#- 8.6.1: Construct list of validated chemicals
#- 8.6.2: Construct list of validation graphs to dipslay