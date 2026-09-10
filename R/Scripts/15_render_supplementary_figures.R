#* 15: Render Supplementary Figures
# Remove any existing PDFs to prevent corruption from previous failed runs
if (file.exists("Supplementary/Components/Figures/PDF/S1.pdf")) {
  file.remove(list.files("Supplementary/Components/Figures/PDF", pattern = "^S[0-9]", full.names = TRUE))
}
#+ 15.1: Print Supplementary Figures (PNG and PDF)
#! Re-read the validation grobs HERE, not at 00c. 00c_FTs.R loads validation_plot_metadata_ordered$grob at run start, but script 09 rewrites Outputs/Validation/revised/grobs/ in the same run whenever run_validation_step is true -- so every end-to-end rebuild compiled the supplement from the PREVIOUS run's grobs. Found 2026-09-10: script 09 had written the corrected p.28 Menthone axis title at 16:38 and the supplement built at 17:04 still carried the old one. Harmless only when the grobs did not change between runs.
validation_plot_metadata_ordered$grob <- purrr::map(validation_plot_metadata_ordered$full_path, readRDS)
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
