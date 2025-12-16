#' Validation Plot Adjustment Function
#'
#' Applies comprehensive adjustments to a validation plot:
#' global formatting, remove standard, zoom x-axis, zoom y-axis, zoom fragments.
#' All parameters are optional - only specified adjustments will be applied.
#'
#' @param plot_obj Plot object with 'plot' field
#' @param remove_std Logical, remove standard peak (default: FALSE)
#' @param xl Lower x-axis limit (requires xu)
#' @param xu Upper x-axis limit (requires xl)
#' @param yl Lower y-axis limit (default: 0, requires yu)
#' @param yu Upper y-axis limit (requires yl)
#' @param mz_fragment Single fragment or vector of fragments to zoom to/remove.
#'   Positive values = keep only those fragments, negative = remove those fragments.
#'   Can mix positive and negative (e.g., c(1, 2, -3) = keep 1,2 and remove 3)
#' @param title_add Optional text to append to existing title in parentheses (e.g., "Fragment")
#' @param subfolder Subfolder for write_small output (default: "revised")
#' @param save_grob Logical, whether to save plot as grob RDS (default: FALSE)
#' @param grob_dir Directory to save grobs (default: "Outputs/Validation/revised_grobs")
#'
#' @return Modified plot object with all specified adjustments applied
#' @export
vp <- function(plot_obj, 
               remove_std = FALSE,
               xl = NULL, xu = NULL,
               yl = 0, yu = NULL,
               mz_fragment = NULL,
               title_add = NULL,
               subfolder = "revised",
               save_grob = TRUE,
               grob_dir = "Outputs/Validation/revised_grobs") {
  
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
  
  modified_plot <- plot_obj
  
  #- Step 0: Apply global formatting (colors, y-axis title, grid removal, scientific notation)
  cat("→ Applying global formatting...\n")
  
  # Define color palette (9 colors for fragments mz0-mz8)
  cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  names(cbp2) <- paste0("mz", 0:8)
  
  # Colors (from adjust_VP_colors logic)
  if ("mz_label" %in% names(modified_plot$plot$data)) {
    present_mz <- unique(modified_plot$plot$data$mz_label)
    present_mz <- present_mz[!is.na(present_mz)]
    
    # Extract mz numbers from full labels (e.g., "mz2: 98.0964" -> "mz2")
    mz_numbers <- gsub("^(mz[0-9]+):.*", "\\1", present_mz)
    
    # Get colors for each fragment based on mz number
    colors_to_use <- cbp2[mz_numbers]
    names(colors_to_use) <- present_mz  # Name with full labels for matching
    
    modified_plot$plot <- modified_plot$plot + 
      ggplot2::scale_colour_manual(values = colors_to_use)
  }
  
  # Y-axis title, grid removal, scientific notation
  modified_plot$plot <- modified_plot$plot +
    ggplot2::labs(y = "← Standard | Tumor →") +
    ggplot2::theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    ggplot2::scale_y_continuous(labels = scales::scientific)
  
  #- Step 1: Remove standard if requested
  if (remove_std) {
    if (is.null(xl) || is.null(xu)) {
      stop("remove_std requires xl and xu to be specified for x-axis range.")
    }
    cat("→ Removing standard peak...\n")
    modified_plot <- remove_standard(modified_plot, xl = xl, xu = xu, subfolder = subfolder, write_output = FALSE)
    # xl, xu consumed by remove_standard, don't apply zoom_x again
    xl <- NULL
    xu <- NULL
  }
  
  #- Step 2: Zoom x-axis if specified
  if (!is.null(xl) && !is.null(xu)) {
    cat("→ Zooming x-axis...\n")
    modified_plot <- zoom_x(modified_plot, xl = xl, xu = xu, subfolder = subfolder, write_output = FALSE)
    
    # Apply dynamic x-axis tick spacing based on range width
    x_range <- xu - xl
    if (x_range <= 0.1) {
      # Narrow range: use finer ticks (every 0.025)
      tick_breaks <- seq(xl, xu, by = 0.025)
      modified_plot$plot <- modified_plot$plot + 
        ggplot2::scale_x_continuous(breaks = tick_breaks, limits = c(xl, xu))
    } else {
      # Wider range: use standard ticks (every 0.04)
      tick_breaks <- seq(xl, xu, by = 0.04)
      modified_plot$plot <- modified_plot$plot + 
        ggplot2::scale_x_continuous(breaks = tick_breaks, limits = c(xl, xu))
    }
  }
  
  #- Step 3: Zoom fragments if specified (before y-axis so y-axis has final say)
  if (!is.null(mz_fragment)) {
    cat("→ Adjusting fragments...\n")
    for (frag in mz_fragment) {
      modified_plot <- zoom_fragment(modified_plot, mz_fragment = frag, subfolder = subfolder, write_output = FALSE)
    }
  }
  
  #- Step 4: Zoom y-axis if specified (applied last to override fragment auto-scaling)
  if (!is.null(yu)) {
    cat("→ Zooming y-axis...\n")
    modified_plot <- zoom_y(modified_plot, yl = yl, yu = yu, subfolder = subfolder, write_output = FALSE)
  }
  
  #- Step 5: Add to title if specified
  if (!is.null(title_add)) {
    cat("→ Updating title...\n")
    current_title <- modified_plot$plot$labels$title
    new_title <- paste0(current_title, " (", title_add, ")")
    modified_plot$plot <- modified_plot$plot + ggplot2::ggtitle(new_title)
  }
  
  #- Step 6: Add analyzed fragment to subtitle (the one with asterisk)
  if ("mz_label" %in% names(modified_plot$plot$data)) {
    asterisk_fragment <- grep("\\*\\*\\\\\\*\\*\\*", modified_plot$plot$data$mz_label, value = TRUE)
    if (length(asterisk_fragment) > 0) {
      # Extract just the mz number (e.g., "mz4: 0.1234 **\***" -> "mz4")
      mz_num <- gsub("^(mz\\d+):.*", "\\1", asterisk_fragment[1])
      current_subtitle <- modified_plot$plot$labels$subtitle
      
      # Only add if "Analyzed Fragment" is not already in subtitle
      if (!is.null(current_subtitle) && !is.na(current_subtitle) && !grepl("Analyzed Fragment:", current_subtitle)) {
        if (current_subtitle != "") {
          new_subtitle <- paste0(current_subtitle, " | Analyzed Fragment: ", mz_num)
        } else {
          new_subtitle <- paste0("Analyzed Fragment: ", mz_num)
        }
        modified_plot$plot <- modified_plot$plot + ggplot2::labs(subtitle = new_subtitle)
        cat(sprintf("→ Added analyzed fragment (%s) to subtitle\n", mz_num))
      } else if (is.null(current_subtitle) || is.na(current_subtitle)) {
        new_subtitle <- paste0("Analyzed Fragment: ", mz_num)
        modified_plot$plot <- modified_plot$plot + ggplot2::labs(subtitle = new_subtitle)
        cat(sprintf("→ Added analyzed fragment (%s) to subtitle\n", mz_num))
      }
    }
  }
  
  #- Step 7: Update legend formatting (always 2 rows for consistent vertical spacing and scientific y-axis)
  cat("→ Updating legend formatting...\n")
  
  # Always use 2 rows to ensure consistent vertical spacing across all plots
  modified_plot$plot <- modified_plot$plot +
    ggplot2::guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    ggplot2::theme(
      legend.text = ggtext::element_markdown(size = 5),
      legend.spacing.x = unit(0.02, "cm"),
      legend.box.margin = margin(0, 0, -20, 0)  # Keep top margin 0, negative bottom to pull plot closer
    )
  
  #- Step 7.5: Replace γ with "gamma" in title for PDF compatibility
  if (!is.null(modified_plot$plot$labels$title)) {
    modified_plot$plot$labels$title <- gsub("γ", "gamma", modified_plot$plot$labels$title)
  }
  
  #- Write final output (with _F suffix if fragment adjustment was used)
  suffix <- if (!is.null(mz_fragment)) "_F" else ""
  write_small(modified_plot, subfolder = subfolder, suffix = suffix)
  
  #- Save as grob if requested (to nested grobs/ directory)
  if (save_grob) {
    grob_output_dir <- file.path("Outputs/Validation", subfolder, "grobs")
    dir.create(grob_output_dir, recursive = TRUE, showWarnings = FALSE)
    plot_tag <- modified_plot$plot_tag
    grob_filename <- paste0(plot_tag, suffix, ".rds")
    grob_path <- file.path(grob_output_dir, grob_filename)
    saveRDS(ggplotGrob(modified_plot$plot), grob_path)
    cat(sprintf("→ Saved grob: %s\n", grob_path))
  }
  
  cat("✓ All adjustments complete\n")
  return(modified_plot)
}
