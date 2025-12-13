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
  
  # Validate input
  if (!is.list(plot_obj) || !"plot" %in% names(plot_obj)) {
    stop("Invalid plot object. Expected a plot object with 'plot' field.")
  }
  
  modified_plot <- plot_obj
  
  #- Step 0: Apply global formatting (colors, y-axis title, grid removal, scientific notation)
  cat("→ Applying global formatting...\n")
  
  # Colors (from adjust_VP_colors logic)
  if ("mz_label" %in% names(modified_plot$plot$data)) {
    mz_labels <- unique(modified_plot$plot$data$mz_label)
    color_mapping <- setNames(color_palette[1:length(mz_labels)], mz_labels)
    modified_plot$plot <- modified_plot$plot + 
      ggplot2::scale_colour_manual(values = color_mapping)
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
      if (!is.null(current_subtitle) && current_subtitle != "") {
        new_subtitle <- paste0(current_subtitle, " | Analyzed Fragment: ", mz_num)
      } else {
        new_subtitle <- paste0("Analyzed Fragment: ", mz_num)
      }
      modified_plot$plot <- modified_plot$plot + ggplot2::labs(subtitle = new_subtitle)
      cat(sprintf("→ Added analyzed fragment (%s) to subtitle\n", mz_num))
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
