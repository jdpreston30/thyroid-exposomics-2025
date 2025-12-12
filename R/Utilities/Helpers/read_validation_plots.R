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
  
  total_files <- length(valid_paths)
  cat(sprintf("\n▶ Reading %d validation plots from OneDrive (%s)...\n\n", 
              total_files, ifelse(parallel, "parallel", "sequential")))
  
  start_time <- Sys.time()
  
  #- Read files
  if (parallel) {
    #+ Parallel processing with progress
    cl <- makeCluster(8)
    clusterExport(cl, c("valid_paths"), envir = environment())
    
    completed <- 0
    validation_plots <- parLapply(cl, valid_paths, function(path) {
      plot <- readRDS(path)
      return(plot)
    })
    stopCluster(cl)
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("\r[Parallel] ✓ Completed: %d plots | Total time: %02d:%02d:%02d\n\n",
                total_files, floor(elapsed / 3600), floor((elapsed %% 3600) / 60), floor(elapsed %% 60)))
    
  } else {
    #+ Sequential processing with progress bar
    validation_plots <- list()
    
    for (i in seq_along(valid_paths)) {
      #- Read file
      validation_plots[[i]] <- readRDS(valid_paths[i])
      
      #- Update progress
      current_time <- Sys.time()
      elapsed <- as.numeric(difftime(current_time, start_time, units = "secs"))
      
      if (i > 0) {
        avg_time_per_file <- elapsed / i
        remaining_files <- total_files - i
        eta_secs <- avg_time_per_file * remaining_files
        eta_str <- sprintf("%02d:%02d:%02d", 
                          floor(eta_secs / 3600), 
                          floor((eta_secs %% 3600) / 60), 
                          floor(eta_secs %% 60))
      } else {
        eta_str <- "calculating..."
      }
      
      elapsed_str <- sprintf("%02d:%02d:%02d", 
                            floor(elapsed / 3600), 
                            floor((elapsed %% 3600) / 60), 
                            floor(elapsed %% 60))
      
      #- Progress bar
      pct <- i / total_files
      bar_width <- 50
      filled <- floor(pct * bar_width)
      bar <- paste0(paste0(rep("█", filled), collapse = ""), 
                    paste0(rep("░", bar_width - filled), collapse = ""))
      
      #- Current file info
      current_file <- basename(valid_paths[i])
      
      #- Single line progress update
      cat(sprintf("\r[%s] %3.0f%% | File %d/%d | %s | Time: %s | ETA: %s          ", 
                  bar, pct * 100, i, total_files,
                  current_file, elapsed_str, eta_str))
      
      flush.console()
    }
    
    cat("\n\n")
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
