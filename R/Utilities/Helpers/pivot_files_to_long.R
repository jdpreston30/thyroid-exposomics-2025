pivot_files_to_long <- function(data) {
  #' Pivot file columns to long format for consistent ordering
  #' 
  #' Converts wide format with file1, file2, ..., fileN columns into long format
  #' where each row represents one file for one compound, with its specific RT values
  #' 
  #' @param data Tibble with file1-fileN, f1_rt-fN_rt, f1_rt_range-fN_rt_range columns
  #' @return Tibble in long format with columns: id, short_name, monoisotopic, all_files, 
  #'         file, file_rt, file_rt_range, source, asterisk, mz0-mz9, s1_rt, s2_rt, s1_rt_range, s2_rt_range
  
  # Get all file column names (file1, file2, ..., fileN)
  file_cols <- names(data)[grepl("^file[0-9]+$", names(data))]
  
  # Extract file numbers
  file_nums <- as.numeric(gsub("file", "", file_cols))
  
  # Pivot each set of file/rt/rt_range columns
  file_data_list <- list()
  
  for (file_num in file_nums) {
    file_col <- paste0("file", file_num)
    rt_col <- paste0("f", file_num, "_rt")
    rt_range_col <- paste0("f", file_num, "_rt_range")
    
    # Check if columns exist
    if (all(c(file_col, rt_col, rt_range_col) %in% names(data))) {
      temp <- data |>
        select(
          id, short_name, monoisotopic, all_files, 
          file = !!file_col, 
          file_rt = !!rt_col, 
          file_rt_range = !!rt_range_col,
          starts_with("mz"), s1_rt, s2_rt, s1_rt_range, s2_rt_range, 
          source, asterisk
        ) |>
        filter(!is.na(file))  # Remove rows where file is NA
      
      file_data_list[[length(file_data_list) + 1]] <- temp
    }
  }
  
  # Combine all file data
  result <- bind_rows(file_data_list) |>
    # Group by id to ensure consistent file ordering within each compound
    group_by(id) |>
    arrange(id, file) |>
    ungroup()
  
  return(result)
}
