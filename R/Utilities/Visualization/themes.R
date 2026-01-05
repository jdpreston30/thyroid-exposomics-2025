#* Visualization Themes and Color Palettes
#' Defines consistent color schemes and themes for thyroid cancer variant analysis
#' Colors are set for Follicular, FV-PTC, and Papillary variants

#+ Variant Color Palette
variant_colors <- c(
  "Follicular" = "#294B88",
  "FV-PTC" = "#23744E",
  "Papillary" = "#DF8D09"
)

#+ IARC Classification Color Palette
IARC_colors <- c(
  "Group 1" = "#BF2D39",
  "Group 2A" = "#CF5349",
  "Group 2B" = "#EAA269",
  "Group 3" = "#FAD97E",
  "Not Classified" = "#FFFFFF"
)

#+ EDC Classification Color Palette
EDC_colors <- c(
  "Non-EDC" = "#7CA0D4",
  "Potential EDC" = "#BADD86"
)

#+ Carcinogen Classification Color Palette
carcinogen_colors <- c(
  "Known Carcinogen" = "#BF2D39",
  "Likely Carcinogen" = "#CF5349",
  "Possible Carcinogen" = "#EAA269",
  "Uncertain Risk" = "#FAD97E"
)

#+ Tumor vs Non-Cancer Color Palette
tumor_noncancer_colors <- c(
  "Thyroid Tumor Tissue" = "#BE4E4D",
  "Non-Cancer Cadaver Thyroids" = "#0C5EA5"
)

#+ Color definitions are kept above - used by other functions
