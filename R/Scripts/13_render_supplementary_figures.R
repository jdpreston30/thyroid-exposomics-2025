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
#+ 13.2: Print Supplementary Figures (PNG and PDF)
#- 13.2.1: Supplementary Figure 1
print_to_png_pdf(sup_fig1, "S1")
#- 13.2.2: Supplementary Figure 2.1
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "2.1")
#- 13.2.2: Supplementary Figure 2.2
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "2.2")
#- 13.2.2: Supplementary Figure 3.1
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "3.1")
#- 13.2.2: Supplementary Figure 3.2
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "3.2")
#- 13.2.2: Supplementary Figure 4.1
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "4.1")
#- 13.2.2: Supplementary Figure 4.2
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "4.2")
