validation_plot_metadata <- read_xlsx(config$paths$validation, sheet = "figure_order")
#* 10: Validation Plots Organization
#+ 10.0: Notes
# The workflow was first to create the combined PDFs from 08. These were inspected visually and any obvious failures were assigned state = "failed" in the validation.xlsx file. Alternatively, any that clearly passed on first inspection were assigned "final".
# A quality score was given on the validation sheet- either F for fail, 1 for level 1, 2 for level 2. An asterisk on those indicates that there was an alternative ID for the same chemical which looked better (or state = "alternative").
# Level 1 quality was considered co-elution of 2 or more fragments of standard or sample within a reasonable time of each other. Far far spread in time was not considered a pass but close co-elution that could be reasonably explained by matrix effects was. Note that there were select cases where the quantiative fragment (asterisk, see below) could not be seen in the standard but could be seen in the sample, and the standard and sample had other co-eluting fragments. These were still considered level 1.
# Single fragments, even if co-eluting in standard and sample, were considered fails.
# Level 2 quality was when you could detect co-elution of 2 or more fragments in the sample but not the standard, or the standard co-eluting fragments were shifted extremely far away.
# The peak (asterisk) used for statistics/quant was considered closely in this process, such that if you zoomed in, even if it was a smaller peak in all of the co-eluting fragments, it still had to follow the same general trajectory as its other fragments. If it could really only be considered noise, then it was considered a fail.
#+ 10.1: Load and prepare validation plot metadata
validation_plot_metadata_ordered <- validation_plot_metadata %>%
  arrange(order) %>%
  mutate(
    sf_sub = paste(figure, subfigure, sep = "."),
    full_path = here::here(plot),
    grob = map(full_path, readRDS)
  ) %>%
  select(order, id, short_name, figure, subfigure, sf_sub, panel, plot, full_path, grob)
#+ 10.2: Compile validation plots into PDFs
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "2.1")
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "2.2")
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "3.1")
compile_sf_sub_pdf(validation_plot_metadata_ordered, sf_sub_value = "3.2")
