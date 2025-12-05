#* Test
source("R/Utilities/Validation/pvc_rtx_multi.R")
pvc_rtx_multi("CP2002", "BL_12082022_037", 
              mz_manual = c(132.0445, 104.0495),
              rtr = c(14.97, 15.02),
              ppm_tolerance = 10, 
              standard = FALSE,
              stick = FALSE,
              png_name = "test_manual",
              spectrum = TRUE,
              rt = 14.99)
#* Other
source("R/Utilities/Validation/pvc_mzx.R")
s2.3.1 <- pvc_mzx("CP2002", "BL_12082022_037", "BP2-1_1", rt = 15.0, rt_window = 0.4, mzr = c(100, 150), ppm_filter = 10, png_name = "test_guthion")
source("R/Utilities/Validation/pvc_rtx.R")
s2.3.1 <- pvc_rtx("CP2002", "BL_12082022_037", rtr = c(14.97, 15.02), ppm_tolerance = 10, png_name = "test_guthion", standard = FALSE, stick = TRUE, max_i = TRUE)