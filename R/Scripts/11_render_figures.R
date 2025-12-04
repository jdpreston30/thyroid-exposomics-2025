#* 5 Render Figures
#+ 5.1: Figure 1
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
#+ 5.2: Figure 2
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
#+ 5.3: Figure 3
fig3 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # 1A
  draw_plot(p3A,
    x = 1.18, y = 6.774,
    width = 2.98, height = 3.055
  ) +
  # 1B
  draw_plot(p3B,
    x = 4.23, y = 7.285,
    width = 3.16, height = 2.6
  ) +
  # 1C 
  # ggsave("balloon3C.svg", plot = p3C, width = 4.4, height = 3, dpi = 1200)
  draw_plot(p3C,
    x = 4.55, y = 7.016,
    width = 2.7945, height = 2.03
  ) +
  # 1D
  draw_plot(p3D,
    x = 4.55, y = 5.016,
    width = 2.8, height = 2.03
  ) +
  # 1E
  draw_plot(p3E,
    x = 0.53, y = 3.699,
    width = 4.4, height = 3
  ) +
  # Labels
  figure_labels(list(
    A = c(1.26, 10.00),
    B = c(4.35, 10.00),
    C = c(4.6, 8.97),
    D = c(4.6, 6.97),
    E = c(0.75, 10.00),
    "Figure 3" = c(0.49, 10.43)
  ))
#+ 5.4: Supplementary Figure 1
#- 5.4.1: Align plots
aligned_plots <- plot_grid(s1A, s1B, 
                           ncol = 2, 
                           align = "hv",
                           axis = "tblr")
#- 5.4.2: Create Supplementary Figure 1
sup_fig1 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  draw_plot(aligned_plots,
    x = -0.1, y = 4,
    width = 7.9, height = 5.98
  ) +
  # Labels
  figure_labels(list(
    A = c(0.75, 10.00),
    B = c(4.25, 10.00)
  ))
#+ 5.5: Print All Main Figures
print_to_png(fig1, "Fig1.png", output_dir = "Outputs/Figures")
print_to_png(fig2, "Fig2.png", output_dir = "Outputs/Figures")
#+ 5.6: Print Supplementary Figure 1 (PNG and PDF)
#- 5.6.1: Print as PNG
print_to_png(sup_fig1, "S1.png", output_dir = "Supplementary/Components/Figures/PNG")
#- 5.6.2: Print as PDF
image_write(image_read("Supplementary/Components/Figures/PNG/S1.png"), "Supplementary/Components/Figures/PDF/S1.pdf", format = "pdf", density = 600)
