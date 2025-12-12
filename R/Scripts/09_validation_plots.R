#* 9: Validation Plots Adjustment and Manual Review
#+ 9.4: Manually copy over files, read in, adjust x ranges
#!!!!!
validation_check <- read_xlsx(config$paths$variant_validation, sheet = "validation")
#- 9.4.0: Read in manual validation results metadata
validation_check_files <- validation_check |>
  filter(!state %in% c("failed", "not used")) |>
  mutate(rt_range = (rtu-rtl)/2) |>
  select(-c(modification, note, rtl, rtu)) |>
  arrange(order)
#- 9.4.1: Derive a list of all unique plots to read in
variant_plot_list <- validation_check_files %>%
  filter(source != "IARC") %>%
  pull(plot) %>%
  str_split(",\\s*") %>%
  unlist() %>%
  unique()
#- 9.3.2: Set up information for copy from OneDrive
{
  ggplot_raw_dir <- "Outputs/Validation/ggplot_objects_raw"
  dir.create(ggplot_raw_dir, showWarnings = FALSE, recursive = TRUE)
  onedrive_base <- config$paths$validation_plot_directory_onedrive
  remaining_plots <- variant_plot_list
  copied_count <- 0
}
#- 9.3.3: If else to pull all files from variant_rtx and iarc_tumor_rtx
{
  # Search in variant_rtx subfolder
  variant_rtx_dir <- file.path(onedrive_base, "variant_rtx")
  for (plot_name in variant_plot_list) {
    source_file <- file.path(variant_rtx_dir, paste0(plot_name, ".rds"))
    if (file.exists(source_file)) {
      dest_file <- file.path(ggplot_raw_dir, paste0(plot_name, ".rds"))
      file.copy(source_file, dest_file, overwrite = TRUE)
      remaining_plots <- setdiff(remaining_plots, plot_name)
      copied_count <- copied_count + 1
      cat(sprintf("  ✓ Copied from variant_rtx: %s\n", plot_name))
    }
  }
  cat(sprintf("\nCopied %d files from variant_rtx\n", copied_count))
  # Search in iarc_tumor_rtx subfolder for remaining plots
  if (length(remaining_plots) > 0) {
    iarc_count <- 0
    iarc_tumor_dir <- file.path(onedrive_base, "iarc_tumor_rtx")
    for (plot_name in remaining_plots) {
      source_file <- file.path(iarc_tumor_dir, paste0(plot_name, ".rds"))
      if (file.exists(source_file)) {
        dest_file <- file.path(output_dir, paste0(plot_name, ".rds"))
        file.copy(source_file, dest_file, overwrite = TRUE)
        remaining_plots <- setdiff(remaining_plots, plot_name)
        iarc_count <- iarc_count + 1
        cat(sprintf("  ✓ Copied from iarc_tumor_rtx: %s\n", plot_name))
      }
    }
    cat(sprintf("\nCopied %d files from iarc_tumor_rtx\n", iarc_count))
  }
  # Report results
  if (length(remaining_plots) == 0) {
    cat("\n✅ Matched all plots!\n")
  } else {
    cat(sprintf("\n⚠️  Could not find %d plots:\n", length(remaining_plots)))
    print(remaining_plots)
  }
}
#- 9.3.4: Read all copied RDS files into a single object (parallel)
{
  rds_files <- list.files(ggplot_raw_dir, pattern = "\\.rds$", full.names = TRUE)
  cat(sprintf("Reading %d RDS files in parallel...\n", length(rds_files)))
  cl <- makeCluster(8)
  validation_plots <- parLapply(cl, rds_files, function(rds_file) {
    tryCatch({
      plot_tag <- tools::file_path_sans_ext(basename(rds_file))
      plot_data <- readRDS(rds_file)
      list(tag = plot_tag, data = plot_data, success = TRUE)
    }, error = function(e) {
      list(tag = tools::file_path_sans_ext(basename(rds_file)), 
           data = NULL, success = FALSE, error = e$message)
    })
  })
  stopCluster(cl)
  names(validation_plots) <- sapply(validation_plots, function(x) x$tag)
  validation_plots <- lapply(validation_plots, function(x) x$data)
  failed <- sapply(validation_plots, is.null)
  if (any(failed)) {
    warning(sprintf("⚠️  Failed to load %d files", sum(failed)))
  }
  validation_plots <- validation_plots[!failed]
  cat(sprintf("✓ Loaded %d plots\n", length(validation_plots)))
}
#- 9.3.4: Adjust x-axis RT ranges for each plot
validation_plots_adjusted <- adjust_validation_plot_ranges(
  validation_plots = validation_plots,
  validation_curated = validation_check_files
)
#+ 9.4: Manual Adjustment of Specific Plots
source("R/Utilities/Helpers/remove_standard.R")
source("R/Utilities/Helpers/write_small.R")
#- 8.4.1: Set vector of plots that need adjustment
adjust_ids <- c("CP2382", "CP3007", "CP2486", "CP2212", "CP1090", "CP3113", "CP3193", "CP3182", "CP3148", "CP2365", "CP3174", "CP3017", "CP3021", "CP2487", "CP1016", "CP2107", "CP3066")
#- 8.3.4: Pull all adjust plots for each ID into separate lists
aps <- list()
for (id in adjust_ids) {
  id_plots <- validation_plots_adjusted[grepl(paste0("_", id, "$"), names(validation_plots_adjusted))]
  if (length(id_plots) > 0) {
    aps[[id]] <- id_plots
    cat(sprintf("Pulled %d plots for %s\n", length(id_plots), id))
  }
}
#+ 8.4: Plots that failed after review
#- 8.4.1: Molinate (CP2486)
F3_S1_CP2486_revised <- remove_standard(aps$CP2486$F1_S1_CP2486, xl = 5.1, xu = 5.175)
#- 8.4.2: Atrazine (CP3113)
write_small(aps$CP3113)
#- 8.4.16: Bupirimate (CP2107)
write_small(aps$CP2107)
#- 8.4.10: Resmethrin (CP3174)
write_small(aps$CP3174)
#+ 8.6: Plots that required revision
#- 8.4.1: MDA (CP3007)
#! removed standard, zoomed Y axis, calling level 2
F3_S1_CP3007_revised <- remove_standard(aps$CP3007$F3_S1_CP3007, xl = 3.8, xu = 4)
#- 8.4.5: N-MeFOSAA (CP3193)
#- 8.4.3: MEHP
#- 8.4.7: TEEP (CP3182)
#- 8.4.8: Menthone (CP3148)
#- 8.4.8: Menthone (CP3148)
#- 8.4.9 : Prosulfuron (CP2365)
#- 8.4.10: Resmethrin (CP3174)
#- 8.4.11: o-Toluidine (CP3017)
#- 8.4.12: o-Anisidine (CP3021)
#- 8.4.13: Vernolate (CP2487)
#- 8.4.16: Bupirimate (CP2107)
#- 8.4.17: o-Cresol (CP3066)
#+ 8.5 Save 
