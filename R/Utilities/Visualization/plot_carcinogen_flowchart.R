#' Render the Carcinogenicity Classification Decision-Tree Flowchart
#'
#' Builds the decision-tree schematic that documents the logic of
#' \code{\link{classify_carcinogenicity}} and writes it to PDF (and PNG) for
#' inclusion in the supplementary material as Figure S1. The diagram is a
#' Graphviz (DOT) flowchart rendered via DiagrammeR; node colors are locked to
#' the project \code{carcinogen_colors} / \code{IARC_colors} palettes so the
#' figure matches the other supplementary figures.
#'
#' @param output_dir_pdf Character path to the PDF output directory
#'   (default: "Supplementary/Components/Figures/PDF"). The figure is written as
#'   \code{S1.pdf}.
#' @param output_dir_png Character path to the PNG output directory
#'   (default: "Supplementary/Components/Figures/PNG"). The figure is also
#'   written as \code{S1.png} for quick preview / parity with print_to_png_pdf().
#' @param font_family Character font family for all labels (default: "Arial",
#'   matching the ggplot figures).
#' @param dpi Numeric resolution for the PNG copy (default: 600).
#'
#' @details
#' Branch order mirrors \code{classify_carcinogenicity()} exactly: decisions
#' cascade down the spine (the "No" path) while each "Yes" branches right to its
#' outcome box. Outcome boxes repeat where two conditions map to the same class
#' (e.g. IARC 2A and H350 >= 50\% both -> Likely). \code{splines=ortho} gives
#' clean right-angle connectors that dock at the diamond tips. Uniform diamond
#' and box dimensions keep the "Yes" arrows equal length and the outcome boxes
#' aligned in a column.
#'
#' The final diamond is phrased as the logical complement of
#' \code{classify_carcinogenicity()}'s actual test ("IARC group and/or GHS
#' statement present?" rather than "neither present?") so that Uncertain Risk
#' -- one of the four classes shown in the stacked bar figure -- reaches the
#' right-hand column via "Yes" like every other bar-chart class. This is the
#' same boolean partition, just restated as its negation: the code's only path
#' to Unclassified requires BOTH no IARC group and no GHS statement, and every
#' other remaining case (including its catch-all) resolves to Uncertain Risk,
#' so "IARC and/or GHS present" exactly (not approximately) separates the two
#' classes with no gap. Unclassified (excluded from the stacked bar) is the
#' sole class left at the bottom of the spine via "No". "and/or" (rather than
#' bare "or") is used here and in \code{q_partial} because both diamonds test
#' multiple sub-conditions that can co-occur in the same chemical's GHS
#' string, not mutually exclusive alternatives.
#'
#' Requires the \code{carcinogen_colors} and \code{IARC_colors} palettes to be
#' in scope (sourced from themes.R by the pipeline, as with the other plot
#' functions).
#'
#' @return Invisibly returns a list with the PDF and PNG file paths.
#'
#' @examples
#' \dontrun{
#'   render_carcinogen_flowchart()
#' }
#'
#' @importFrom DiagrammeR grViz
#' @importFrom DiagrammeRsvg export_svg
#' @importFrom rsvg rsvg_pdf rsvg_png
#' @export
render_carcinogen_flowchart <- function(
    output_dir_pdf = "Supplementary/Components/Figures/PDF",
    output_dir_png = "Supplementary/Components/Figures/PNG",
    font_family = "Arial",
    dpi = 600) {

  # Palette (locked to the bar graphs) ----------------------------------------
  known_fill     <- unname(carcinogen_colors["Known Carcinogen"])     # #BF2D39
  likely_fill    <- unname(carcinogen_colors["Likely Carcinogen"])    # #CF5349
  possible_fill  <- unname(carcinogen_colors["Possible Carcinogen"])  # #EAA269
  uncertain_fill <- unname(carcinogen_colors["Uncertain Risk"])       # #FAD97E
  unclass_fill   <- unname(IARC_colors["Not Classified"])             # #FFFFFF

  # Uniform shape dimensions --------------------------------------------------
  dia_w <- 3.2        # diamond width  (uniform)
  dia_h <- 1.25       # diamond height (uniform; fits up to 3 label lines)
  box_w <- 2.1        # outcome box width  (uniform)
  box_h <- 0.55       # outcome box height (uniform)
  png_w <- 7.5        # PNG width in inches (PDF uses the SVG's natural size)

  # DOT graph -----------------------------------------------------------------
  dot <- '
digraph carcinogenicity {
  graph [rankdir = TB, splines = ortho, nodesep = 0.35, ranksep = 0.45,
         fontname = "Helvetica", bgcolor = "white"]
  node  [fontname = "Helvetica", fontsize = 12,
         labelloc = "c", labeljust = "c"]
  edge  [fontname = "Helvetica", fontsize = 11, color = "#555555",
         arrowsize = 0.8]

  // --- decision nodes (uniform diamonds) ---
  node [shape = diamond, style = "filled", fillcolor = "#EEF2F7",
        color = "#8DA0B8", penwidth = 1, fixedsize = true,
        width = __DIAW__, height = __DIAH__]
  start    [label = "Chemical differing\\nbetween tumor types",
            shape = box, style = "filled,rounded", fillcolor = "#DCE6F1",
            color = "#8DA0B8", fixedsize = false, width = 2.8, height = 0.7]
  q_iarc1  [label = "IARC =\\nGroup 1?"]
  q_iarc2a [label = "IARC =\\nGroup 2A?"]
  q_h350   [label = "H350 and/or H350i\\n≥ 50%?"]
  q_iarc2b [label = "IARC =\\nGroup 2B?"]
  q_partial[label = "0% < H350/H350i\\n< 50%, and/or\\nH351 > 0%?"]
  q_none   [label = "IARC group and/or\\nGHS statement\\npresent?"]

  // --- terminal nodes (uniform boxes), colored by carcinogen_colors ---
  node [shape = box, style = "filled,rounded", penwidth = 0, fixedsize = true,
        width = __BOXW__, height = __BOXH__]
  known     [label = "Known Carcinogen",    fillcolor = "__KNOWN__",    fontcolor = "white"]
  likely1   [label = "Likely Carcinogen",   fillcolor = "__LIKELY__",   fontcolor = "white"]
  likely2   [label = "Likely Carcinogen",   fillcolor = "__LIKELY__",   fontcolor = "white"]
  possible1 [label = "Possible Carcinogen", fillcolor = "__POSSIBLE__", fontcolor = "#333333"]
  possible2 [label = "Possible Carcinogen", fillcolor = "__POSSIBLE__", fontcolor = "#333333"]
  uncertn2  [label = "Uncertain Risk",      fillcolor = "__UNCERTAIN__", fontcolor = "#333333"]
  unclass   [label = "Unclassified",        fillcolor = "__UNCLASS__",  fontcolor = "#333333",
             penwidth = 1, color = "#8DA0B8"]

  // --- "Yes" label holders: labelloc=t floats the text above the arrow; width
  //     ~0 so the two split-arrow segments meet with no visible white gap. ---
  node [shape = plaintext, style = "", fixedsize = true, width = 0.001,
        height = 0.5, labelloc = "t", fontsize = 11, fontcolor = "#000000"]
  yl1 [label = "Yes"]
  yl2 [label = "Yes"]
  yl3 [label = "Yes"]
  yl4 [label = "Yes"]
  yl5 [label = "Yes"]
  yl7 [label = "Yes"]

  // --- spine: the "No" path straight down ---
  start    -> q_iarc1
  q_iarc1  -> q_iarc2a  [label = "    No"]
  q_iarc2a -> q_h350    [label = "    No"]
  q_h350   -> q_iarc2b  [label = "    No"]
  q_iarc2b -> q_partial [label = "    No"]
  q_partial-> q_none    [label = "    No"]
  q_none   -> unclass   [label = "    No"]

  // --- "Yes" branches out to the right ---
  // headclip/tailclip=false -> both segments run to the label node centre and
  // meet there, so the split arrow is truly gapless (no white speck).
  q_iarc1  -> yl1 [arrowhead = none, headclip = false]   yl1 -> known     [tailclip = false]
  q_iarc2a -> yl2 [arrowhead = none, headclip = false]   yl2 -> likely1   [tailclip = false]
  q_h350   -> yl3 [arrowhead = none, headclip = false]   yl3 -> likely2   [tailclip = false]
  q_iarc2b -> yl4 [arrowhead = none, headclip = false]   yl4 -> possible1 [tailclip = false]
  q_partial-> yl5 [arrowhead = none, headclip = false]   yl5 -> possible2 [tailclip = false]
  q_none   -> yl7 [arrowhead = none, headclip = false]   yl7 -> uncertn2  [tailclip = false]

  { rank = same; q_iarc1;  yl1; known }
  { rank = same; q_iarc2a; yl2; likely1 }
  { rank = same; q_h350;   yl3; likely2 }
  { rank = same; q_iarc2b; yl4; possible1 }
  { rank = same; q_partial; yl5; possible2 }
  { rank = same; q_none;   yl7; uncertn2 }
}
'

  # Inject font, sizes, palette (sub/gsub avoids sprintf clashing with "%") ----
  dot <- gsub("Helvetica",     font_family,    dot, fixed = TRUE)
  dot <- gsub("__DIAW__",      dia_w,          dot, fixed = TRUE)
  dot <- gsub("__DIAH__",      dia_h,          dot, fixed = TRUE)
  dot <- gsub("__BOXW__",      box_w,          dot, fixed = TRUE)
  dot <- gsub("__BOXH__",      box_h,          dot, fixed = TRUE)
  dot <- sub("__KNOWN__",      known_fill,     dot, fixed = TRUE)
  dot <- gsub("__LIKELY__",    likely_fill,    dot, fixed = TRUE)
  dot <- gsub("__POSSIBLE__",  possible_fill,  dot, fixed = TRUE)
  dot <- gsub("__UNCERTAIN__", uncertain_fill, dot, fixed = TRUE)
  dot <- sub("__UNCLASS__",    unclass_fill,   dot, fixed = TRUE)

  # Render --------------------------------------------------------------------
  if (!dir.exists(output_dir_pdf)) dir.create(output_dir_pdf, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir_png)) dir.create(output_dir_png, recursive = TRUE, showWarnings = FALSE)

  svg <- DiagrammeRsvg::export_svg(DiagrammeR::grViz(dot))
  raw <- charToRaw(svg)

  pdf_path <- file.path(output_dir_pdf, "S1.pdf")
  png_path <- file.path(output_dir_png, "S1.png")

  # PDF uses the SVG's natural (vector) size; LaTeX scales it to the page.
  rsvg::rsvg_pdf(raw, file = pdf_path)
  # PNG width only -> height follows the aspect ratio (no stretch).
  rsvg::rsvg_png(raw, file = png_path, width = round(png_w * dpi))

  message("Created: ", pdf_path)
  invisible(list(pdf = pdf_path, png = png_path))
}
