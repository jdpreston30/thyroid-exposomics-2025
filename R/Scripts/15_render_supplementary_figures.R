#* 15: Render Supplementary Figures
# Remove any existing PDFs to prevent corruption from previous failed runs
if (file.exists("Supplementary/Components/Figures/PDF/S1.pdf")) {
  file.remove(list.files("Supplementary/Components/Figures/PDF", pattern = "^S[0-9]", full.names = TRUE))
}
#+ 15.1: Print Supplementary Figures (PNG and PDF)
#- 15.1.1: Supplementary figure 1 — carcinogenicity classification decision tree
#! Standalone Graphviz schematic of classify_carcinogenicity(); writes S1.pdf
render_carcinogen_flowchart()
#- 15.1.2: Supplementary figure 2.1 (validation; formerly 1.1)
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "2.1")
#- 15.1.3: Supplementary figure 2.2 (validation; formerly 1.2)
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "2.2")
#- 15.1.4: Supplementary figure 3.1 (validation; formerly 2.1)
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "3.1")
#- 15.1.5: Supplementary figure 3.2 (validation; formerly 2.2)
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "3.2")
#- 15.1.6: Supplementary figure 4.1 (validation; formerly 3.1)
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "4.1")
#- 15.1.7: Supplementary figure 4.2 (validation; formerly 3.2)
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "4.2")
