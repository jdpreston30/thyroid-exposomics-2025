#* 5 Render Figures
#+ 5.1: Figure 1
#- 5.1.1: Align all plots as a 2x2 grid to ensure both x and y axis alignment
aligned_plots <- align_plots(p1A, p1B, align = "hv", axis = "tblr")
#- 5.1.2: Render Figure 1
fig1 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # 1A
  draw_plot(aligned_plots[[1]], x = 0.57, y = 7.41, width = 4, height = 2.4) +
  # 1B
  draw_plot(aligned_plots[[2]], x = 4.7, y = 7.41, width = 3, height = 2.4) +
  # Labels
  figure_labels(list(
    A = c(0.74, 10.00),
    B = c(4.87, 10.00),
    "Figure 1" = c(0.49, 10.43)
  ))
#+ 5.6: Print All Main Figures
#- 5.7.1: Figure 1
print_to_png(fig1 + grdgd(),
  "Fig1.png",
  width = 8.5, height = 11, dpi = 300, output_dir = "Outputs/Figures"
)