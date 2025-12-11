source("R/Utilities/Validation/finalize_validation_plots.R")
finalized_validation <- finalize_validation_plots(
  validation_curated = validation_check_curate,
  variant_rtx_folder = "variant_rtx",
  iarc_tumor_rtx_folder = "iarc_tumor_rtx",
  config = config,
  output_dir = "Outputs/Validation",
  pdf_name = "finalized_validation.pdf"
)