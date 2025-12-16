#' Validation Plot Adjustment Function for Top Fragments
#'
#' Simplified version of vp() that automatically:
#' - Filters to only the top fragment for the compound
#' - Removes standards
#' - Centers x-axis on RT ± 0.16666
#'
#' @param plot_obj Plot object with 'plot' field
#' @param frag_lookup Tibble with columns: id, top_frag (fragment number to keep). 
#'   Default: filters validation_check_files to non-NA top_frag values
#' @param yl Lower y-axis limit (default: 0, requires yu)
#' @param yu Upper y-axis limit (requires yl)
#' @param title_add Optional text to append to existing title in parentheses
#' @param subfolder Subfolder for write_small output (default: "revised")
#' @param save_grob Logical, whether to save plot as grob RDS (default: TRUE)
#'
#' @return Modified plot object with all specified adjustments applied
#' @export
vp_top <- function(plot_obj, 
                   frag_lookup = validation_check_files |> 
                     filter(!is.na(top_frag)) |> 
                     select(id, top_frag),
                   yl = 0, 
                   yu = NULL,
                   title_add = NULL,
                   subfolder = "revised",
                   save_grob = TRUE) {
  
  # If plot_obj is a symbol name (unquoted), load it first
  plot_arg <- substitute(plot_obj)
  if (is.symbol(plot_arg)) {
    plot_name <- deparse(plot_arg)
    
    # Check if already exists in global environment
    if (exists(plot_name, envir = .GlobalEnv)) {
      cat(sprintf("✓ Using existing plot: %s\n", plot_name))
      plot_obj <- get(plot_name, envir = .GlobalEnv)
      
      # Fix plot_tag if filename has C_ prefix but internal plot_tag doesn't
      if (startsWith(plot_name, "C_") && !startsWith(plot_obj$plot_tag, "C_")) {
        plot_obj$plot_tag <- paste0("C_", plot_obj$plot_tag)
        assign(plot_name, plot_obj, envir = .GlobalEnv)
        cat(sprintf("  → Updated plot_tag to include C_ prefix\n"))
      }
    } else {
      cat(sprintf("📂 Loading plot: %s\n", plot_name))
      
      # Search through directories in order
      base_path <- config$paths$validation_plot_directory_onedrive
      search_dirs <- c(
        "variant_rtx",
        "iarc_tumor_rtx",
        "iarc_cadaver_rtx"
      )
      
      rds_path <- NULL
      for (dir in search_dirs) {
        test_path <- file.path(base_path, dir, paste0(plot_name, ".rds"))
        if (file.exists(test_path)) {
          rds_path <- test_path
          cat(sprintf("  Found in: %s\n", dir))
          break
        }
      }
      
      if (is.null(rds_path)) {
        stop(sprintf("Plot file not found in any directory: %s", plot_name))
      }
      
      plot_obj <- readRDS(rds_path)
      
      # Fix plot_tag if loaded filename has C_ prefix but internal plot_tag doesn't
      if (startsWith(plot_name, "C_") && !startsWith(plot_obj$plot_tag, "C_")) {
        plot_obj$plot_tag <- paste0("C_", plot_obj$plot_tag)
        cat(sprintf("  → Updated plot_tag to include C_ prefix\n"))
      }
      
      assign(plot_name, plot_obj, envir = .GlobalEnv)
      cat(sprintf("✅ Loaded %s from RDS\n", plot_name))
    }
  }
  
  # Validate input
  if (!is.list(plot_obj) || !"plot" %in% names(plot_obj)) {
    stop("Invalid plot object. Expected a plot object with 'plot' field.")
  }
  
  # Extract compound ID and look up top fragment
  compound_id <- plot_obj$id
  top_frag_row <- frag_lookup[frag_lookup$id == compound_id, ]
  
  if (nrow(top_frag_row) == 0) {
    stop(sprintf("Compound ID %s not found in fragment lookup table", compound_id))
  }
  
  top_frag <- top_frag_row$top_frag[1]
  cat(sprintf("→ Filtering to fragment %d for %s\n", top_frag, compound_id))
  
  # Extract RT range and calculate center
  rt_range <- plot_obj$rt_range
  rt_center <- mean(rt_range)
  xl <- rt_center - 0.16666
  xu <- rt_center + 0.16666
  cat(sprintf("→ Setting x-axis: RT center = %.4f, range = [%.4f, %.4f]\n", rt_center, xl, xu))
  
  modified_plot <- plot_obj
  
  #- Step 0: Apply global formatting
  cat("→ Applying global formatting...\n")
  
  # Define color palette
  cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  names(cbp2) <- paste0("mz", 0:8)
  
  # Colors
  if ("mz_label" %in% names(modified_plot$plot$data)) {
    present_mz <- unique(modified_plot$plot$data$mz_label)
    present_mz <- present_mz[!is.na(present_mz)]
    mz_numbers <- gsub("^(mz[0-9]+):.*", "\\1", present_mz)
    colors_to_use <- cbp2[mz_numbers]
    names(colors_to_use) <- present_mz
    modified_plot$plot <- modified_plot$plot + 
      ggplot2::scale_colour_manual(values = colors_to_use)
  }
  
  # Y-axis title, grid removal, scientific notation
  modified_plot$plot <- modified_plot$plot +
    ggplot2::labs(y = "Tumor") +
    ggplot2::theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    ggplot2::scale_y_continuous(labels = scales::scientific)
  
  #- Step 1: Remove standard (always true for vp_top)
  cat("→ Removing standard peak...\n")
  modified_plot <- remove_standard(modified_plot, xl = xl, xu = xu, subfolder = subfolder, write_output = FALSE)
  
  #- Step 2: Filter to top fragment only
  cat(sprintf("→ Filtering to fragment %d...\n", top_frag))
  modified_plot <- zoom_fragment(modified_plot, mz_fragment = top_frag, subfolder = subfolder, write_output = FALSE)
  
  #- Step 3: Zoom y-axis if specified
  if (!is.null(yu)) {
    cat("→ Zooming y-axis...\n")
    modified_plot <- zoom_y(modified_plot, yl = yl, yu = yu, subfolder = subfolder, write_output = FALSE)
  }
  
  #- Step 4: Add to title if specified
  if (!is.null(title_add)) {
    cat("→ Updating title...\n")
    current_title <- modified_plot$plot$labels$title
    new_title <- paste0(current_title, " (", title_add, ")")
    modified_plot$plot <- modified_plot$plot + ggplot2::ggtitle(new_title)
  }
  
  #- Step 5: Update legend formatting
  cat("→ Updating legend formatting...\n")
  modified_plot$plot <- modified_plot$plot +
    ggplot2::guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
    ggplot2::theme(
      legend.text = ggtext::element_markdown(size = 5),
      legend.spacing.x = unit(0.02, "cm"),
      legend.box.margin = margin(0, 0, -20, 0)
    )
  
  #- Step 6: Replace γ with "gamma" in title for PDF compatibility
  if (!is.null(modified_plot$plot$labels$title)) {
    modified_plot$plot$labels$title <- gsub("γ", "gamma", modified_plot$plot$labels$title)
  }
  
  #- Write final output with _TOP suffix
  write_small(modified_plot, subfolder = subfolder, suffix = "_TOP")
  
  #- Save as grob if requested
  if (save_grob) {
    grob_output_dir <- file.path("Outputs/Validation", subfolder, "grobs")
    dir.create(grob_output_dir, recursive = TRUE, showWarnings = FALSE)
    plot_tag <- modified_plot$plot_tag
    grob_filename <- paste0(plot_tag, "_TOP.rds")
    grob_path <- file.path(grob_output_dir, grob_filename)
    saveRDS(ggplotGrob(modified_plot$plot), grob_path)
    cat(sprintf("→ Saved grob: %s\n", grob_path))
  }
  
  cat("✓ All adjustments complete\n")
  return(modified_plot)
}
