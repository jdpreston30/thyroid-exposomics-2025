#* 0b: Configuration Setup
#+ 0b.1: Set Up R Options and Repositories
options(expressions = 10000)
#+ 0b.2: Load Utility Functions First (needed for Dynamic config)
utils_path <- "R/Utilities/"
if (dir.exists(utils_path)) {
  purrr::walk(
    list.files(utils_path, pattern = "\\.[rR]$", full.names = TRUE, recursive = TRUE),
    source
  )
  cat("🔧 Loaded utility functions\n")
}
#+ 0b.3: Load Dynamic Project Configuration
.GlobalEnv$CONFIG <- config
#+ 0b.4: Set Up Global Paths from Config
output_path <- config$paths$output  
scripts_path <- config$paths$scripts
#+ 0b.5: Create Output Directory If It doesn't Exist
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
  cat("📁 Created output directory:", output_path, "\n")
}
#+ 0b.6: Set Up R Environment Preferences
#- 0b.6.1: Tibble preferences
options(
  tibble.print_max = config$analysis$tibble_options$print_max,
  tibble.print_min = config$analysis$tibble_options$print_min,
  pillar.sigfig = config$analysis$tibble_options$sigfig
)
#- 0b.6.2: Data.table preferences
if (!is.null(config$analysis$datatable_options)) {
  options(
    datatable.print.class = config$analysis$datatable_options$print_class,
    datatable.print.keys = config$analysis$datatable_options$print_keys
  )
  .datatable.aware = config$analysis$datatable_options$aware
}
cat("✅ Configuration and environment setup complete!\n")
