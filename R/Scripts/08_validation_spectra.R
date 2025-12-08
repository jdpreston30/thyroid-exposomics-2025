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
vv_wide_temp <- vv_wide |>
  filter(order <= 1)
vv_wide_temp |> select(id, f1_rt, f2_rt, f3_rt)
#+ 8.2: Manual Validation (Variant Differences)
source("R/Utilities/Validation/rtx.R")
variant_validation <- rtx(
  validation_list = vv_wide_temp,
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

# #- 8.2.1: 2,6-DCP-4'-NPE
# s2.1.1 <- pvc_rtx("CP2458", "BL_12082022_004", "BP2-1_1", rtr = c(10.02, 10.78), ppm_tolerance = 5, png_name = "quant_2_6-DCP-4'-NPE")
# #- 8.2.2: Cyfluthrin
# s2.2.1 <- pvc_rtx("CP3153", "BL_12082022_021", "BP3-1_1", rtr = c(8.01, 8.77), ppm_tolerance = 5, png_name = "quant_Cyfluthrin")
# #- 8.2.3: Guthion
# s2.3.1 <- pvc_rtx("CP2002", "BL_12082022_001", "BP2-1_2", rtr = c(14.72, 15.48), ppm_tolerance = 5, png_name = "quant_Guthion")
# #- 8.2.4: o-Toluidine
# s2.4.1 <- pvc_rtx("CP3017", "BL_12082022_093", source = "quant", "BP3-1_2", rtr = c(7.49, 8.25), ppm_tolerance = 5, png_name = "quant_o-Toluidine")
# #- 8.2.5: DNOP
# s2.5.1 <- pvc_rtx("CP2187", "BL_12082022_011", "BP2-1_1", rtr = c(17.42, 18.18), ppm_tolerance = 5, png_name = "quant_DNOP")
# #+ 8.3: Manual spectral extraction (rtx IARC)
# #- 8.3.2: o-Toluidine
# #_Cadaver
# s3.2.1 <- pvc_rtx("CP3017", "BL_08222024_ThyroidTssue_003", "BP3-1_2", study = "cadaver", source = "IARC", rtr = c(7.49, 8.25), ppm_tolerance = 5, png_name = "iarc_cadaver_toluidine")
# #_Tumor
# s3.1.1 <- pvc_rtx("CP3017", "BL_12082022_093", "BP3-1_2", study = "tumor", source = "IARC", rtr = c(7.49, 8.25), ppm_tolerance = 5, png_name = "iarc_tumor_toluidine")
# #+ 8.6: Modify MT_final per validation results
# MT_final <- MT_final_i
# #+ 8.7: Modify cadaver/tumor comparisons by validation results
# cadaver_iarc_keep <- cadaver_iarc |>
#   filter(cas == "95-53-4") |>
#   pull(name_sub_lib_id)
# #- 8.6.2: full_joiner_i
# full_joiner <- full_joiner_i |>
#   select(sample_ID, variant, tumor_vs_ctrl, any_of(cadaver_iarc_keep))
  