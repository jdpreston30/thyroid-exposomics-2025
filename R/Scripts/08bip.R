#* Get all RDS files in validation plots directory
validation_path <- "/Users/JoshsMacbook2015/Library/CloudStorage/OneDrive-EmoryUniversity/Research/MS_raw_data/GCMS/Thyroid_GCMS/DTC_2022/validation_plots"
rds_files <- list.files(path = validation_path, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
rds_names <- tibble(file_name = basename(rds_files) %>% str_remove("\\.rds$"))
#* Check file counts per CP ID
file_counts <- rds_names %>%
  mutate(cp_id = str_extract(file_name, "CP\\d{4}")) %>%
  count(cp_id, name = "n_files")
#* Join with validation_check_files to see completeness
validation_file_check <- validation_check_files %>%
  left_join(file_counts, by = c("id" = "cp_id")) %>%
  mutate(n_files = replace_na(n_files, 0)) %>%
  select(id, short_name, n_files, state) %>%
  arrange(desc(n_files), id)
