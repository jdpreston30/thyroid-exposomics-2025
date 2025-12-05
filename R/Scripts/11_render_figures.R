#* 11 Render Figures
#+ 11.1: Figure 1
fig1 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # 1A
  draw_plot(p1A, x = 1.18, y = 6.774, 
                 width = 2.98, height = 3.055) +
  # 1B
  draw_plot(p1B, x = 4.23, y = 7.285,
                 width = 3.16, height = 2.6) +
  # Labels
  figure_labels(list(
    A = c(1.26, 10.00),
    B = c(4.35, 10.00),
    "Figure 1" = c(0.49, 10.43)
  ))
#+ 11.2: Figure 2
fig2 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # 2A
  draw_plot(p2A,
    x = 0.53, y = 3.699,
    width = 4.6, height = 6.28
  ) +
  # 2B
  draw_plot(p2B,
    x = 4.55, y = 7.016,
    width = 2.7945, height = 2.03
  ) +
  # 2C
  draw_plot(p2C,
    x = 4.55, y = 5.016,
    width = 2.8, height = 2.03
  ) +
  # Labels
  figure_labels(list(
    A = c(0.75, 10.00),
    B = c(4.6, 8.97),
    C = c(4.6, 6.97),
    "Figure 2" = c(0.49, 10.43)
  ))
#+ 11.3: Figure 3
# fig3 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
#   # 1A
#   draw_plot(p3A,
#     x = 1.18, y = 6.774,
#     width = 2.98, height = 3.055
#   ) +
#   # 1B
#   draw_plot(p3B,
#     x = 4.23, y = 7.285,
#     width = 3.16, height = 2.6
#   ) +
#   # 1C 
#   # ggsave("balloon3C.svg", plot = p3C, width = 4.4, height = 3, dpi = 1200)
#   draw_plot(p3C,
#     x = 4.55, y = 7.016,
#     width = 2.7945, height = 2.03
#   ) +
#   # 1D
#   draw_plot(p3D,
#     x = 4.55, y = 5.016,
#     width = 2.8, height = 2.03
#   ) +
#   # 1E
#   draw_plot(p3E,
#     x = 0.53, y = 3.699,
#     width = 4.4, height = 3
#   ) +
#   # Labels
#   figure_labels(list(
#     A = c(1.26, 10.00),
#     B = c(4.35, 10.00),
#     C = c(4.6, 8.97),
#     D = c(4.6, 6.97),
#     E = c(0.75, 10.00),
#     "Figure 3" = c(0.49, 10.43)
#   ))
#+ 11.4: Supplementary Figure 1
#- 11.4.1: Align plots
s1_aligned <- plot_grid(s1A, s1B, 
                           ncol = 2, 
                           align = "hv",
                           axis = "tblr")
#- 11.4.2: Create Supplementary Figure 1
sup_fig1 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  draw_plot(s1_aligned,
    x = -0.1, y = 4,
    width = 7.9, height = 5.98
  ) +
  # Labels
  figure_labels(list(
    A = c(0.75, 10.00),
    B = c(4.25, 10.00)
  ))
#+ 11.5: Supplementary Figure 2 (Page 1)
sup_fig2.1 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # Row 1: 2,6-DCP-4'-NPE
  draw_plot(s2.1.1, x = 0.555, y = 6.774, width = 3.9, height = 3.25) +
  draw_plot(s2.1.2, x = 4.105, y = 6.774, width = 3.9, height = 3.25) +
  # Row 2: Cyfluthrin (shifted down 0.75 total)
  draw_plot(s2.2.1, x = 0.555, y = 3.524, width = 3.9, height = 3.25) +
  draw_plot(s2.2.2, x = 4.105, y = 3.524, width = 3.9, height = 3.25) +
  # Row 3: Guthion (shifted down 1.5 total)
  draw_plot(s2.3.1, x = 0.555, y = 0.274, width = 3.9, height = 3.25) +
  draw_plot(s2.3.2, x = 4.105, y = 0.274, width = 3.9, height = 3.25) +
  # Labels
  figure_labels(list(
    "Supplementary Figure 2.1" = c(0.49, 10.43)
  )) +
  # Chemical name labels
  draw_label("2,6-DCP-4'-NPE", x = 4.33, y = 10.05, size = 14, fontface = "bold.italic", fontfamily = "Arial", hjust = 0.5) +
  geom_segment(aes(x = 3.6, xend = 5.06, y = 9.95, yend = 9.95), size = 0.5) +
  draw_label("Cyfluthrin", x = 4.33, y = 6.8, size = 14, fontface = "bold.italic", fontfamily = "Arial", hjust = 0.5) +
  geom_segment(aes(x = 3.85, xend = 4.81, y = 6.7, yend = 6.7), size = 0.5) +
  draw_label("Guthion", x = 4.33, y = 3.55, size = 14, fontface = "bold.italic", fontfamily = "Arial", hjust = 0.5) +
  geom_segment(aes(x = 3.95, xend = 4.71, y = 3.45, yend = 3.45), size = 0.5)
#+ 11.5b: Supplementary Figure 2 (Page 2)
sup_fig2.2 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # Row 1: o-Toluidine
  draw_plot(s2.4.1, x = 0.555, y = 6.774, width = 3.9, height = 3.25) +
  draw_plot(s2.4.2, x = 4.105, y = 6.774, width = 3.9, height = 3.25) +
  # Row 2: DNOP (shifted down 0.75 total)
  draw_plot(s2.5.1, x = 0.555, y = 3.524, width = 3.9, height = 3.25) +
  draw_plot(s2.5.2, x = 4.105, y = 3.524, width = 3.9, height = 3.25) +
  # Labels
  figure_labels(list(
    "Supplementary Figure 2.2" = c(0.49, 10.43)
  )) +
  # Chemical name labels
  draw_label("o-Toluidine", x = 4.33, y = 10.05, size = 14, fontface = "bold.italic", fontfamily = "Arial", hjust = 0.5) +
  geom_segment(aes(x = 3.8, xend = 4.86, y = 9.95, yend = 9.95), size = 0.5) +
  draw_label("DNOP", x = 4.33, y = 6.8, size = 14, fontface = "bold.italic", fontfamily = "Arial", hjust = 0.5) +
  geom_segment(aes(x = 4.08, xend = 4.58, y = 6.7, yend = 6.7), size = 0.5)
#+ 11.6: Supplementary Figure 3
#- 11.6.2: Create Supplementary Figure 3
sup_fig3 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # Row 1
  draw_plot(s3.1.1, x = 0.555, y = 6.774, width = 3.9, height = 3.25) +
  draw_plot(s3.1.2, x = 4.105, y = 6.774, width = 3.9, height = 3.25) +
  # Row 2
  draw_plot(s3.2.1, x = 0.555, y = 3.524, width = 3.9, height = 3.25) +
  draw_plot(s3.2.2, x = 4.105, y = 3.524, width = 3.9, height = 3.25) +
  # Labels
  figure_labels(list(
    "Supplementary Figure 3" = c(0.49, 10.43)
  )) +
  # Chemical name labels
  draw_label("o-Toluidine", x = 4.33, y = 10.05, size = 14, fontface = "bold.italic", fontfamily = "Arial", hjust = 0.5) +
  geom_segment(aes(x = 3.8, xend = 4.86, y = 9.95, yend = 9.95), size = 0.5) +
  draw_label("o-Toluidine", x = 4.33, y = 6.8, size = 14, fontface = "bold.italic", fontfamily = "Arial", hjust = 0.5) +
  geom_segment(aes(x = 3.8, xend = 4.86, y = 6.7, yend = 6.7), size = 0.5)
#+ 11.7: Print All Main Figures
print_to_png(fig1, "Fig1.png", output_dir = "Outputs/Figures")
print_to_png(fig2, "Fig2.png", output_dir = "Outputs/Figures")
#+ 11.8: Print Supplementary Figures (PNG and PDF)
#- 11.8.1: Supplementary Figure 1
print_to_png_pdf(sup_fig1, "S1")
#- 11.8.2: Supplementary Figure 2a
print_to_png_pdf(sup_fig2.1, "S2.1")
#- 11.8.3: Supplementary Figure 2b
print_to_png_pdf(sup_fig2.2, "S2.2")
#- 11.8.4: Supplementary Figure 3
print_to_png_pdf(sup_fig3, "S3")
#+ 11.9: Combine S1, S2.1, S2.2 into single PDF

renv::snapshot()
renv::clean() 