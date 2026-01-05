#* UNUSED COLOR/SCALE FUNCTIONS FROM themes.R
#* Moved to utilities_hold for testing pipeline without them

#+ Helper Functions (UNUSED - Moved from themes.R)
#' Apply variant colors to ggplot
#' @param discrete Logical, if TRUE uses scale_color_manual, if FALSE uses scale_fill_manual
#' @return ggplot2 scale object
scale_variant_color <- function(discrete = TRUE) {
  if (discrete) {
    scale_color_manual(values = variant_colors)
  } else {
    scale_fill_manual(values = variant_colors)
  }
}

#' Apply variant fill colors to ggplot
#' @return ggplot2 scale object
scale_variant_fill <- function() {
  scale_fill_manual(values = variant_colors)
}

#' Get color for specific variant
#' @param variant_name Character string of variant name
#' @return Hex color code
get_variant_color <- function(variant_name) {
  variant_colors[variant_name]
}

#' Apply IARC classification colors to ggplot
#' @return ggplot2 scale object
scale_IARC_color <- function() {
  scale_color_manual(values = IARC_colors)
}

#' Apply IARC classification fill colors to ggplot
#' @return ggplot2 scale object
scale_IARC_fill <- function() {
  scale_fill_manual(values = IARC_colors)
}

#' Get color for specific IARC group
#' @param group_name Character string of IARC group name
#' @return Hex color code
get_IARC_color <- function(group_name) {
  IARC_colors[group_name]
}

#' Apply EDC classification colors to ggplot
#' @return ggplot2 scale object
scale_EDC_color <- function() {
  scale_color_manual(values = EDC_colors)
}

#' Apply EDC classification fill colors to ggplot
#' @return ggplot2 scale object
scale_EDC_fill <- function() {
  scale_fill_manual(values = EDC_colors)
}

#' Get color for specific EDC classification
#' @param edc_class Character string of EDC classification
#' @return Hex color code
get_EDC_color <- function(edc_class) {
  EDC_colors[edc_class]
}

#' Apply carcinogen classification colors to ggplot
#' @return ggplot2 scale object
scale_carcinogen_color <- function() {
  scale_color_manual(values = carcinogen_colors)
}

#' Apply carcinogen classification fill colors to ggplot
#' @return ggplot2 scale object
scale_carcinogen_fill <- function() {
  scale_fill_manual(values = carcinogen_colors)
}

#' Get color for specific carcinogen classification
#' @param carc_class Character string of carcinogen classification
#' @return Hex color code
get_carcinogen_color <- function(carc_class) {
  carcinogen_colors[carc_class]
}

#' Apply tumor vs non-cancer colors to ggplot
#' @return ggplot2 scale object
scale_tumor_noncancer_color <- function() {
  scale_color_manual(values = tumor_noncancer_colors)
}

#' Apply tumor vs non-cancer fill colors to ggplot
#' @return ggplot2 scale object
scale_tumor_noncancer_fill <- function() {
  scale_fill_manual(values = tumor_noncancer_colors)
}

#' Get color for tumor or non-cancer tissue
#' @param tissue_type Character string of tissue type
#' @return Hex color code
get_tumor_noncancer_color <- function(tissue_type) {
  tumor_noncancer_colors[tissue_type]
}

#+ Base Theme (UNUSED)
#' Custom ggplot2 theme for publication-ready figures
#' Standardized formatting: Arial font, size 10 bold for axis text, size 12 bold for titles
theme_thyroid <- function() {
  theme_minimal(base_family = "Arial") +
    theme(
      # Axis text: size 10, bold
      axis.text.x = element_text(size = 10, face = "bold", family = "Arial"),
      axis.text.y = element_text(size = 10, face = "bold", family = "Arial"),
      # Axis titles: size 12, bold
      axis.title.x = element_text(size = 12, face = "bold", family = "Arial"),
      axis.title.y = element_text(size = 12, face = "bold", family = "Arial"),
      # Legend text: size 10, plain (not bold)
      legend.text = element_text(size = 10, face = "plain", family = "Arial"),
      legend.title = element_text(size = 10, face = "bold", family = "Arial"),
      legend.position = "right",
      # Other elements
      plot.title = element_text(size = 16, face = "bold", family = "Arial"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 12, face = "bold", family = "Arial")
    )
}
