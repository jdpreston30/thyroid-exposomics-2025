#' Read Validation Plots Directly from OneDrive
#'
#' Searches OneDrive validation plot directories and reads RDS files directly
#' without local copying. Searches variant_rtx first, then iarc_tumor_rtx.
#'
#' @param plot_names Character vector of plot names (without .rds extension)
#' @param onedrive_base_path Path to OneDrive validation plot directory
#' @param parallel Logical, whether to use parallel processing (default TRUE)
#'
#' @return Named list of ggplot objects
#' @export
read_validation_plots <- function(plot_names, onedrive_base_path, parallel = TRUE) {
  
  #- Setup directories
  variant_rtx_dir <- file.path(onedrive_base_path, "variant_rtx")
  iarc_tumor_dir <- file.path(onedrive_base_path, "iarc_tumor_rtx")
  
  #- Build file path list
  file_paths <- sapply(plot_names, function(plot_name) {
    variant_file <- file.path(variant_rtx_dir, paste0(plot_name, ".rds"))
    if (file.exists(variant_file)) return(variant_file)
    iarc_file <- file.path(iarc_tumor_dir, paste0(plot_name, ".rds"))
    if (file.exists(iarc_file)) return(iarc_file)
    return(NA_character_)
  })
  
  #- Check for missing files
  missing <- is.na(file_paths)
  if (any(missing)) {
    warning(sprintf("⚠️  Could not find %d plots:\n", sum(missing)))
    print(plot_names[missing])
  }
  
  #- Filter to existing files
  valid_paths <- file_paths[!missing]
  valid_names <- plot_names[!missing]
  
  cat(sprintf("Reading %d validation plots from OneDrive (%s)...\n", 
              length(valid_paths), ifelse(parallel, "parallel", "sequential")))
  
  #- Read files
  if (parallel) {
    cl <- makeCluster(8)
    validation_plots <- parLapply(cl, valid_paths, readRDS)
    stopCluster(cl)
  } else {
    validation_plots <- lapply(valid_paths, readRDS)
  }
  
  #- Name list
  names(validation_plots) <- valid_names
  
  #- Report results
  variant_count <- sum(grepl("variant_rtx", valid_paths))
  iarc_count <- sum(grepl("iarc_tumor_rtx", valid_paths))
  cat(sprintf("✅ Read %d plots from variant_rtx\n", variant_count))
  cat(sprintf("✅ Read %d plots from iarc_tumor_rtx\n", iarc_count))
  
  if (length(not_found) > 0) {
    warning(sprintf("⚠️  Could not find %d plots:\n", length(not_found)))
    print(not_found)
  } else {
    cat("✅ Found all requested plots!\n")
  }
  
  return(validation_plots)
}
