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
#- 8.2.0: Pull RT ranges of top 5
#_ Pulling in data
bp1 <- read_csv("/Users/jdp2019/Library/CloudStorage/OneDrive-Emory/Research/Manuscripts and Projects/Active Projects/Thyroid Exposomics Paper/MANUSCRIPT SUBMISSIONS/Thyroid Lancet/Other/Sami_data_12_5_24/main/BP1.GC2/feature.sample.i.csv")
bp2 <- read_csv("/Users/jdp2019/Library/CloudStorage/OneDrive-Emory/Research/Manuscripts and Projects/Active Projects/Thyroid Exposomics Paper/MANUSCRIPT SUBMISSIONS/Thyroid Lancet/Other/Sami_data_12_5_24/main/BP2.GC2/feature.sample.i.csv")
bp3 <- read_csv("/Users/jdp2019/Library/CloudStorage/OneDrive-Emory/Research/Manuscripts and Projects/Active Projects/Thyroid Exposomics Paper/MANUSCRIPT SUBMISSIONS/Thyroid Lancet/Other/Sami_data_12_5_24/main/BP3.GC2/feature.sample.i.csv")
#_Binding
bp2_bp2 <- rbind(
  bp2 |> filter(id %in% c("CP2002", "CP2187", "CP2458")),
  bp3 |> filter(id %in% c("CP3153", "CP3017"))) |>
  select(id, name, trtlower, trtupper)
#- 8.2.1: 2,6-DCP-4'-NPE
s2.1.1 <- pvc_rtx("CP2458", "BL_12082022_004", "BP2-1_1", rtr = c(10.02, 10.78), ppm_tolerance = 5, png_name = "quant_2_6-DCP-4'-NPE")
#- 8.2.2: Cyfluthrin
s2.2.1 <- pvc_rtx("CP3153", "BL_12082022_021", "BP3-1_1", rtr = c(8.01, 8.77), ppm_tolerance = 5, png_name = "quant_Cyfluthrin")
#- 8.2.3: Guthion
s2.3.1 <- pvc_rtx("CP2002", "BL_12082022_001", "BP2-1_2", rtr = c(14.72, 15.48), ppm_tolerance = 5, png_name = "quant_Guthion")
#- 8.2.4: o-Toluidine
s2.4.1 <- pvc_rtx("CP3017", "BL_12082022_093", source = "quant", "BP3-1_2", rtr = c(7.49, 8.25), ppm_tolerance = 5, png_name = "quant_o-Toluidine")
#- 8.2.5: DNOP
s2.5.1 <- pvc_rtx("CP2187", "BL_12082022_011", "BP2-1_1", rtr = c(17.42, 18.18), ppm_tolerance = 5, png_name = "quant_DNOP")
#+ 8.3: Manual spectral extraction (rtx IARC)
#- 8.3.2: o-Toluidine
#_Cadaver
s3.2.1 <- pvc_rtx("CP3017", "BL_08222024_ThyroidTssue_003", "BP3-1_2", study = "cadaver", source = "IARC", rtr = c(7.49, 8.25), ppm_tolerance = 5, png_name = "iarc_cadaver_toluidine")
#_Tumor
s3.1.1 <- pvc_rtx("CP3017", "BL_12082022_093", "BP3-1_2", study = "tumor", source = "IARC", rtr = c(7.49, 8.25), ppm_tolerance = 5, png_name = "iarc_tumor_toluidine")
#+ 8.4: Manual spectral extraction (mzx Quant)
#- 8.4.1: 2,6-DCP-4'-NPE
s2.1.2 <- pvc_mzx("CP2458", "BL_12082022_004", "BP2-1_1", rt = 10.4, rt_window = 0.38, ppm_filter = 5, png_name = "quant_2_6-DCP-4'-NPE")
#- 8.4.2: Cyfluthrin
s2.2.2 <- pvc_mzx("CP3153", "BL_12082022_021", "BP3-1_1", rt = 8.39, rt_window = 0.38, ppm_filter = 5, block_label = "top", png_name = "quant_Cyfluthrin")
#- 8.4.3: Guthion
s2.3.2 <- pvc_mzx("CP2002", "BL_12082022_001", "BP2-1_2", rt = 15.1, rt_window = 0.38, ppm_filter = 5, png_name = "quant_Guthion")
#- 8.4.4: o-Toluidine
s2.4.2 <- pvc_mzx("CP3017", "BL_12082022_093", "BP3-1_2", source = "quant",
        rt = 7.87, rt_window = 0.38, 
        ppm_filter = 5,
        block_label = "top",
        png_name = "quant_o-Toluidine")
#- 8.4.5: DNOP
s2.5.2 <- pvc_mzx("CP2187", "BL_12082022_011", "BP2-1_1", rt = 17.8, rt_window = 0.38, ppm_filter = 5, block_label = "top", png_name = "quant_DNOP")
#+ 8.5: Manual spectral extraction (mzx IARC)
#- 8.5.1: o-Toluidine
#_Cadaver
s3.2.2 <- pvc_mzx("CP3017", "BL_08222024_ThyroidTssue_003", "BP3-1_2", source = "IARC", study = "cadaver", rt = 7.87, rt_window = 0.38, ppm_filter = 5, png_name = "iarc_cadaver_toluidine_mz")
#_Tumor
s3.1.2 <- pvc_mzx("CP3017", "BL_12082022_093", "BP3-1_2", source = "IARC", study = "tumor", rt = 7.87, rt_window = 0.38, ppm_filter = 5, block_label = "top", png_name = "iarc_tumor_toluidine")
#+ 8.6: Modify MT_final per validation results
MT_final <- MT_final_i
#+ 8.7: Modify cadaver/tumor comparisons by validation results
cadaver_iarc_keep <- cadaver_iarc |>
  filter(cas == "95-53-4") |>
  pull(name_sub_lib_id)
#- 8.6.2: full_joiner_i
full_joiner <- full_joiner_i |>
  select(sample_ID, variant, tumor_vs_ctrl, any_of(cadaver_iarc_keep))
  