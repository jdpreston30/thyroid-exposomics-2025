#* 13: Render Supplementary Figures
# Remove any existing PDFs to prevent corruption from previous failed runs
if (file.exists("Supplementary/Components/Figures/PDF/S1.pdf")) {
  file.remove(list.files("Supplementary/Components/Figures/PDF", pattern = "^S[0-9]", full.names = TRUE))
}
#+ 13.1: Print Supplementary Figures (PNG and PDF)
#- 13.1.1: Supplementary Figure 1.1
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "1.1")
#- 13.1.2: Supplementary Figure 1.2
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "1.2")
#- 13.1.3: Supplementary Figure 2.1
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "2.1")
#- 13.1.4: Supplementary Figure 2.2
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "2.2")
#- 13.1.5: Supplementary Figure 3.1
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "3.1")
#- 13.1.6: Supplementary Figure 3.2
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "3.2")
