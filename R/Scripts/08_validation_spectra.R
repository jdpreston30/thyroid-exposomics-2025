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
#+ 8.2: Manual spectral extraction (Quant)
#- 8.2.1: 2,6-DCP-4'-NPE
pvcd("CP2458", "BL_12082022_004", "BP2-1_1", rtr = c(10.2, 10.5), ppm_tolerance = 1000, png_name = "quant_2_6-DCP-4'-NPE")
#- 8.2.2: o-Toluidine
pvcd("CP3017", "BL_12082022_093", "BP3-1_2", rtr = c(8.05, 8.25),ppm_tolerance = 1000, png_name = "quant_o-Toluidine")
#- 8.2.3: Cyfluthrin
pvcd("CP3153", "BL_12082022_021", "BP3-1_1", rtr = c(8.05, 8.25), ppm_tolerance = 1000, png_name = "quant_Cyfluthrin")
#- 8.2.4: Guthion
pvcd("CP2002", "BL_12082022_001", "BP2-1_2", rtr = c(15.05, 15.2), ppm_tolerance = 1000, png_name = "quant_Guthion")
#- 8.2.4: DNOP
pvcd("CP2187", "BL_12082022_011", "BP2-1_1", rtr = c(17.75, 18.15), ppm_tolerance = 100, png_name = "quant_DNOP")
#+ 8.2: Manual spectral extraction (IARC)
#- 8.2.1: 2-Naphthylamine
#_Cadaver
pvcd("CP3014", "BL_08222024_ThyroidTssue_003", "BP3-1_1", 
     study = "cadaver", 
     source = "IARC",
     rtr = c(7.1, 7.7), 
     ppm_tolerance = 1000, 
     png_name = "iarc_cadaver_2-Naphthylamine")
#_Tumor
pvcd("CP3014", "BL_12082022_115", "BP3-1_2", 
     study = "tumor", 
     source = "IARC",
     rtr = c(6.9, 7.6), 
     ppm_tolerance = 1000, 
     png_name = "iarc_tumor_2-Naphthylamine")
#- 8.2.2: 4-aminobiphenyl
#_Cadaver
pvcd("CP3002", "BL_08222024_ThyroidTssue_007", "BP3-1_2", 
     study = "cadaver", 
     source = "IARC",
     rtr = c(5, 5.6), 
     ppm_tolerance = 1000, 
     png_name = "iarc_cadaver_4-ABP")
#_Tumor
pvcd("CP3002", "BL_12082022_004", "BP3-1_2", 
     study = "tumor", 
     source = "IARC",
     rtr = c(4.9, 5.7), 
     ppm_tolerance = 1000, 
     png_name = "iarc_tumor_4-ABP")
#- 8.2.3: Benz[a]pyrene
#_Cadaver
pvcd("CP3028", "BL_08222024_ThyroidTssue_003", "BP3-1_2", 
     study = "cadaver", 
     source = "IARC",
     rtr = c(16.35, 16.8), 
     ppm_tolerance = 1000, 
     png_name = "iarc_cadaver_Benz[a]pyrene")
#_Tumor
pvcd("CP3028", "BL_12082022_023", "BP3-1_1", 
     study = "tumor", 
     source = "IARC",
     rtr = c(16.2, 17.1), 
     ppm_tolerance = 1000, 
     png_name = "iarc_tumor_Benz[a]pyrene")
#- 8.2.4: Benzidine
#_Cadaver
pvcd("CP2215", "BL_08222024_ThyroidTssue_003", "BP2-1_2", 
     study = "cadaver", 
     source = "IARC",
     rtr = c(9.7, 10.3), 
     ppm_tolerance = 1000, 
     png_name = "iarc_cadaver_benzidine")
#_Tumor
pvcd("CP2215", "BL_12082022_003", "BP2-1_2", 
     study = "tumor", 
     source = "IARC",
     rtr = c(9.7, 10.3), 
     ppm_tolerance = 1000, 
     png_name = "iarc_tumor_benzidine")
#- 8.2.4: o-Toluidine
#_Cadaver
pvcd("CP3017", "BL_08222024_ThyroidTssue_003", "BP3-1_2", 
     study = "cadaver", 
     source = "IARC",
     rtr = c(7.7, 8), 
     ppm_tolerance = 1000, 
     png_name = "iarc_cadaver_toluidine")
#_Tumor
pvcd("CP3017", "BL_12082022_093", "BP3-1_2", 
     study = "tumor", 
     source = "IARC",
     rtr = c(8.05, 8.25),
     ppm_tolerance = 1000, 
     png_name = "iarc_tumor_toluidine")
