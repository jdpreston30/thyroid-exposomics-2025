#* 8: Manual Spectral Validation
#+ 8.1: Convert raw files to mzML
#- 8.1.1: Setup base paths
tumor_raw_dir <- config$paths$tumor_raw_dir
cadaver_raw_dir <- config$paths$cadaver_raw_dir
#- 8.1.2: Convert all .raw files to mzML format
file_inventory <- convert_raw_to_mzml(
  file_list = file_list,
  tumor_raw_dir = tumor_raw_dir,
  cadaver_raw_dir = cadaver_raw_dir
)
#! Function will skip if already completed
#+ 8.2: Manual spectral extraction (rtx Quant)
#- 8.2.1: 2,6-DCP-4'-NPE
s2.1.1 <- pvc_rtx("CP2458", "BL_12082022_004", "BP2-1_1", rtr = c(10.2, 10.5), ppm_tolerance = 1000, png_name = "quant_2_6-DCP-4'-NPE")
#- 8.2.2: Cyfluthrin
s2.2.1 <- pvc_rtx("CP3153", "BL_12082022_021", "BP3-1_1", rtr = c(8.05, 8.25), ppm_tolerance = 1000, png_name = "quant_Cyfluthrin")
#- 8.2.3: Guthion
s2.3.1 <- pvc_rtx("CP2002", "BL_12082022_001", "BP2-1_2", rtr = c(15.05, 15.2), ppm_tolerance = 1000, png_name = "quant_Guthion")
#- 8.2.4: o-Toluidine
s2.4.1 <- pvc_rtx("CP3017", "BL_12082022_093", source = "quant", "BP3-1_2", rtr = c(8.05, 8.25), ppm_tolerance = 1000, png_name = "quant_o-Toluidine")
#- 8.2.5: DNOP
s2.5.1 <- pvc_rtx("CP2187", "BL_12082022_011", "BP2-1_1", rtr = c(17.95, 18.15), ppm_tolerance = 100, png_name = "quant_DNOP")
#+ 8.3: Manual spectral extraction (rtx IARC)
#- 8.3.1: 4-aminobiphenyl
#_Cadaver
pvc_rtx("CP3002", "BL_08222024_ThyroidTssue_007", "BP3-1_2", study = "cadaver", source = "IARC", rtr = c(5, 5.6), ppm_tolerance = 1000, png_name = "iarc_cadaver_4-ABP")
#_Tumor
pvc_rtx("CP3002", "BL_12082022_004", "BP3-1_2", study = "tumor", source = "IARC", rtr = c(4.9, 5.7), ppm_tolerance = 1000, png_name = "iarc_tumor_4-ABP")
#- 8.3.2: o-Toluidine
#_Cadaver
s3.2.1 <- pvc_rtx("CP3017", "BL_08222024_ThyroidTssue_003", "BP3-1_2", study = "cadaver", source = "IARC", rtr = c(7.7, 8), ppm_tolerance = 1000, png_name = "iarc_cadaver_toluidine")
#_Tumor
s3.1.1 <- pvc_rtx("CP3017", "BL_12082022_093", "BP3-1_2", study = "tumor", source = "IARC", rtr = c(8.05, 8.25), ppm_tolerance = 1000, png_name = "iarc_tumor_toluidine")
#- 8.3.3: 2-Naphthylamine
#_Cadaver
pvc_rtx("CP3014", "BL_08222024_ThyroidTssue_003", "BP3-1_1", study = "cadaver", source = "IARC", rtr = c(7.1, 7.7), ppm_tolerance = 1000, png_name = "iarc_cadaver_2-Naphthylamine")
#_Tumor
pvc_rtx("CP3014", "BL_12082022_115", "BP3-1_2", study = "tumor", source = "IARC", rtr = c(6.9, 7.6), ppm_tolerance = 1000, png_name = "iarc_tumor_2-Naphthylamine")
#- 8.3.4: Benz[a]pyrene
#_Cadaver
pvc_rtx("CP3028", "BL_08222024_ThyroidTssue_003", "BP3-1_2", study = "cadaver", source = "IARC", rtr = c(16.35, 16.8), ppm_tolerance = 1000, png_name = "iarc_cadaver_Benz[a]pyrene")
#_Tumor
pvc_rtx("CP3028", "BL_12082022_023", "BP3-1_1", study = "tumor", source = "IARC", rtr = c(16.2, 17.1), ppm_tolerance = 1000, png_name = "iarc_tumor_Benz[a]pyrene")
#- 8.3.5: Benzidine
#_Cadaver
pvc_rtx("CP2215", "BL_08222024_ThyroidTssue_003", "BP2-1_2", study = "cadaver", source = "IARC", rtr = c(9.7, 10.3), ppm_tolerance = 1000, png_name = "iarc_cadaver_benzidine")
#_Tumor
pvc_rtx("CP2215", "BL_12082022_003", "BP2-1_2", study = "tumor", source = "IARC", rtr = c(9.7, 10.3), ppm_tolerance = 1000, png_name = "iarc_tumor_benzidine")
#+ 8.4: Manual spectral extraction (mzx Quant)
#- 8.4.1: 2,6-DCP-4'-NPE
s2.1.2 <- pvc_mzx("CP2458", "BL_12082022_004", "BP2-1_1", rt = 10.35, rt_window = 0.15, ppm_filter = 1000, png_name = "quant_2_6-DCP-4'-NPE")
#- 8.4.2: Cyfluthrin
s2.2.2 <- pvc_mzx("CP3153", "BL_12082022_021", "BP3-1_1", rt = 8.15, rt_window = 0.10, ppm_filter = 1000, block_label = "top", png_name = "quant_Cyfluthrin")
#- 8.4.3: Guthion
s2.3.2 <- pvc_mzx("CP2002", "BL_12082022_001", "BP2-1_2", rt = 15.125, rt_window = 0.075, ppm_filter = 1000, png_name = "quant_Guthion")
#- 8.4.4: o-Toluidine
s2.4.2 <- pvc_mzx("CP3017", "BL_12082022_093", "BP3-1_2", source = "quant",
        rt = 8.2, rt_window = 0.1, 
        ppm_filter = 1000,
        block_label = "top",
        png_name = "quant_o-Toluidine")
#- 8.4.5: DNOP
s2.5.2 <- pvc_mzx("CP2187", "BL_12082022_011", "BP2-1_1", rt = 18.05, rt_window = 0.20, ppm_filter = 100, block_label = "top", png_name = "quant_DNOP")
#+ 8.5: Manual spectral extraction (mzx IARC)
#- 8.5.1: o-Toluidine
#_Cadaver
s3.2.2 <- pvc_mzx("CP3017", "BL_08222024_ThyroidTssue_003", "BP3-1_2", source = "IARC", study = "cadaver", rt = 7.85, rt_window = 0.15, ppm_filter = 1000, png_name = "iarc_cadaver_toluidine_mz")
#_Tumor
s3.1.2 <- pvc_mzx("CP3017", "BL_12082022_093", "BP3-1_2", source = "IARC", study = "tumor", rt = 8.2, rt_window = 0.1, ppm_filter = 1000, block_label = "top", png_name = "iarc_tumor_toluidine")
#+ 8.6: Modify MT_final per validation results
MT_final <- MT_final_i
#+ 8.7: Modify cadaver/tumor comparisons by validation results
cadaver_iarc_keep <- cadaver_iarc |>
  filter(cas == "95-53-4") |>
  pull(name_sub_lib_id)
#- 8.6.2: full_joiner_i
full_joiner <- full_joiner_i |>
  select(sample_ID, variant, tumor_vs_ctrl, any_of(cadaver_iarc_keep))
  