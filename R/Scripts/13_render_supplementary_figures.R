#* 13: Render Supplementary Figures
#+ 13.1: Supplementary Figure 1
#- 13.1.1: Align plots
s1_aligned <- plot_grid(s1A, s1B, 
                           ncol = 2, 
                           align = "hv",
                           axis = "tblr")
#- 13.1.2: Create Supplementary Figure 1
sup_fig1 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  draw_plot(s1_aligned,
    x = -0.6, y = 4.85,
    width = 9.2, height = 5.98
  ) +
  # Labels
  figure_labels(list(
    A = c(0.075, 10.88),
    B = c(4.25, 10.88)
  ))
#+ 13.2: Supplementary Figure 2
#- 13.2.1: Page 1
sup_fig2.1 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # Row 1: 2,6-DCP-4'-NPE
  draw_plot(s2.1.1, x = -0.245, y = 7.229, width = 4.485, height = 3.7375) +
  draw_plot(s2.1.2, x = 4.23, y = 7.229, width = 4.485, height = 3.7375) +
  # Row 2: Cyfluthrin (shifted down 0.75 total)
  draw_plot(s2.2.1, x = -0.245, y = 3.4915, width = 4.485, height = 3.7375) +
  draw_plot(s2.2.2, x = 4.23, y = 3.4915, width = 4.485, height = 3.7375) +
  # Row 3: Guthion (shifted down 1.5 total)
  draw_plot(s2.3.1, x = -0.245, y = -0.246, width = 4.485, height = 3.7375) +
  draw_plot(s2.3.2, x = 4.23, y = -0.246, width = 4.485, height = 3.7375) +
  # Panel labels
  draw_label("A", x = 0.075, y = 10.88, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0.5) +
  draw_label("B", x = 0.075, y = 7.1425, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0.5) +
  draw_label("C", x = 0.075, y = 3.405, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0.5) +
  # Chemical name labels
  draw_label("2,6-DCP-4'-NPE", x = 4.33, y = 10.88, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0.5) +
  geom_segment(aes(x = 3.6, xend = 5.06, y = 10.78, yend = 10.78), size = 0.5) +
  draw_label("Cyfluthrin", x = 4.33, y = 7.1425, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0.5) +
  geom_segment(aes(x = 3.85, xend = 4.81, y = 7.0425, yend = 7.0425), size = 0.5) +
  draw_label("Guthion", x = 4.33, y = 3.405, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0.5) +
  geom_segment(aes(x = 3.95, xend = 4.71, y = 3.305, yend = 3.305), size = 0.5)
#- 13.2.2: Page 2
sup_fig2.2 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # Row 1: o-Toluidine
  draw_plot(s2.4.1, x = -0.245, y = 7.229, width = 4.485, height = 3.7375) +
  draw_plot(s2.4.2, x = 4.23, y = 7.229, width = 4.485, height = 3.7375) +
  # Row 2: DNOP (shifted down 0.75 total)
  draw_plot(s2.5.1, x = -0.245, y = 3.4915, width = 4.485, height = 3.7375) +
  draw_plot(s2.5.2, x = 4.23, y = 3.4915, width = 4.485, height = 3.7375) +
  # Panel labels
  draw_label("D", x = 0.075, y = 10.88, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0.5) +
  draw_label("E", x = 0.075, y = 7.1425, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0.5) +
  # Chemical name labels
  draw_label("o-Toluidine", x = 4.33, y = 10.88, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0.5) +
  geom_segment(aes(x = 3.8, xend = 4.86, y = 10.78, yend = 10.78), size = 0.5) +
  draw_label("DNOP", x = 4.33, y = 7.1425, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0.5) +
  geom_segment(aes(x = 4.08, xend = 4.58, y = 7.0425, yend = 7.0425), size = 0.5)
#+ 13.3: Supplementary Figure 3
sup_fig3 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # Row 1
  draw_plot(s3.1.1, x = -0.245, y = 7.229, width = 4.485, height = 3.7375) +
  draw_plot(s3.1.2, x = 4.23, y = 7.229, width = 4.485, height = 3.7375) +
  # Row 2
  draw_plot(s3.2.1, x = -0.245, y = 3.4915, width = 4.485, height = 3.7375) +
  draw_plot(s3.2.2, x = 4.23, y = 3.4915, width = 4.485, height = 3.7375) +
  # Panel labels
  draw_label("A", x = 0.075, y = 10.88, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0.5) +
  draw_label("B", x = 0.075, y = 7.1425, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0.5) +
  # Chemical name labels
  draw_label("o-Toluidine", x = 4.33, y = 10.88, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0.5) +
  geom_segment(aes(x = 3.8, xend = 4.86, y = 10.78, yend = 10.78), size = 0.5) +
  draw_label("o-Toluidine", x = 4.33, y = 7.1425, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0.5) +
  geom_segment(aes(x = 3.8, xend = 4.86, y = 7.0425, yend = 7.0425), size = 0.5)
#+ 13.4: Print Supplementary Figures (PNG and PDF)
#- 13.4.1: Supplementary Figure 1
print_to_png_pdf(sup_fig1, "S1")
#- 13.4.2: Supplementary Figure 2a
print_to_png_pdf(sup_fig2.1, "S2.1")
#- 13.4.3: Supplementary Figure 2b
print_to_png_pdf(sup_fig2.2, "S2.2")
#- 13.4.4: Supplementary Figure 3
print_to_png_pdf(sup_fig3, "S3")
