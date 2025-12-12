plot_path <- file.path(
  config$paths$validation_plot_directory_onedrive,
  "variant_rtx",
  "F6_S1_CP3182.rds"
)

# Check if file exists
file.exists(plot_path)

# Read the plot (will take ~18-20 seconds due to 80 MB size)
F6_S1_CP3182 <- readRDS(plot_path)
write_small(F6_S1_CP3182)

