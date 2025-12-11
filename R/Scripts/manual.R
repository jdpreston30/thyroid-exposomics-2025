source("R/Utilities/Validation/finalize_validation_plots.R")
finalized_validation <- finalize_validation_plots(
  validation_curated = validation_check_curate,
  variant_rtx_folder = "variant_rtx",
  iarc_tumor_rtx_folder = "iarc_tumor_rtx",
  config = config,
  output_dir = "Outputs/Validation",
  pdf_name = "finalized_validation.pdf"
)

#- 8.3.1: Set vector of problematic plots

#- 8.3.0: Read in manual validation results
validation_check_curate <- validation_check |>
  filter(state == "final") |>
  mutate(rt_range = (rtu-rtl)/2) |>
  select(-c(to_do, note, rtl, rtu)) |>
  arrange(order)
#- 8.3.2: Set vector of problematic plots
validation_check_problems <- validation_check |>
  filter(id %in% problematic_ids)
copy_raw_validation_plots(
  validation_curated = validation_check_curate,
  config = config,
  output_dir = config$paths$validation_plots_raw
)











fragment_ids <- c("CP2487", "CP2382")
source("R/Utilities/Visualization/individual_fragments.R")
dir.create("Outputs/Validation/Intermediate", recursive = TRUE, showWarnings = FALSE)

for (id in fragment_ids) {
  # Pull all plots for this ID
  id_plots <- validation_plots_adjusted[grepl(paste0("_", id, "$"), names(validation_plots_adjusted))]
  
  if (length(id_plots) == 0) {
    cat(sprintf("No plots found for %s - skipping\n", id))
    next
  }
  
  cat(sprintf("\nProcessing %d plots for %s\n", length(id_plots), id))
  
  # Create separated fragment plots for each
  for (plot_name in names(id_plots)) {
    id_separated <- individual_fragments(id_plots[[plot_name]], layout = "vertical")
    ggsave(
      filename = file.path("Outputs/Validation/Intermediate", paste0(plot_name, "_separated.pdf")),
      plot = id_separated,
      width = 8.5,
      height = 11,
      units = "in"
    )
    cat(sprintf("  Saved: %s_separated.pdf\n", plot_name))
  }
}
