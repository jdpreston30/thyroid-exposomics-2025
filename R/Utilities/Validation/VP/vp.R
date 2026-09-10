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
  #- Step 0a: Assert the asterisk marker is well formed
#! 5b82667 (2026-01-05) made process_single_compound.R emit `**LABEL \*****` instead of `LABEL **\***`, which broke the markdown, the `^(mz[0-9]+):` colour lookup (the marked fragment vanished with an NA colour) and the Step 6 "Analyzed Fragment" detector. Fixed at source on 2026-09-04, and every grob the supplement reads was regenerated on 2026-09-10, so the load-time repair that lived here is gone; a malformed label now stops the run instead of being patched over.
  if ("mz_label" %in% names(modified_plot$plot$data) && any(grepl("^\\*\\*", modified_plot$plot$data$mz_label))) {
    stop(sprintf("vp(): %s carries the pre-2026-09-04 asterisk marker (label starts with **). Regenerate the grob; do not patch it here.", modified_plot$plot_tag))
  }
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
#! Cohort is DERIVED, not hardcoded. This line read "Tumor" unconditionally and so relabelled every cadaver plot as tumor -- the grob arrives from process_single_compound.R:622 correctly reading "Cadaver", and this overwrote it. Visible on supplement pp. 34-36 (six labels, GC097 samples, no tumor samples on those pages). Masked from grep because the rendered text is the expanded "Standard", while the source form is the abbreviated "Std". Cadaver grobs carry a C_ prefix on plot_tag, enforced at the load step above.
  .cohort <- if (is.character(modified_plot$plot_tag) &&
                 length(modified_plot$plot_tag) == 1L &&
                 startsWith(modified_plot$plot_tag, "C_")) "Cadaver" else "Tumor"
#! The mirrored title only when the mirror is still there. remove_standard() drops every Standard row and relabels the axis "Intensity (Sample)"; script 09 then calls vp() a second time on that object for the fragment-isolated panel, and this line used to put the mirrored title back on a plot with no standard in it (supplement p.28, Menthone, F5_S1_CP3148_F). The data decides: Standard rows present -> mirrored title; absent -> sample-only title.
  .has_std <- "type" %in% names(modified_plot$plot$data) && any(modified_plot$plot$data$type == "Standard", na.rm = TRUE)
  modified_plot$plot <- modified_plot$plot +
    ggplot2::labs(y = if (.has_std) sprintf("← Standard | %s →", .cohort) else "Intensity (Sample)") +
    ggplot2::theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
#! abs() because the mirror is drawn by negating the standard (process_single_compound.R:523), so the
#! lower half is standard intensity, not a negative quantity. Without it this line overrides the correct
#! labeller set at creation, and it is only repaired if zoom_y, remove_standard or zoom_fragment runs.
#! Identical to the labeller those three already use, so plots that do reach them are unaffected.
    ggplot2::scale_y_continuous(labels = function(x) scales::scientific(abs(x)))
  
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
  
  #- Step 6.5: Assert the subtitle retention time equals the measured value
#! Until 2026-09-10 this step silently REPAIRED the label: the generator printed mean(sample_rt_range) from %.2f-rounded endpoints, and remove_standard()/zoom_fragment() rebuilt it from window midpoints, so up to a third of labels were off by 0.01 (o-Toluidine sat on the 8.085 tie). All three writers now use the measured RT (subtitle_rt / subtitle_rt_of()), and the 2026-09-10 full rebuild produced 0 corrections across all 83 loaded plots with the supplement text identical to the verified build. The lookup is kept as an assertion so a regression fails the run instead of being patched over in silence. Same lookup as before: compound id + sample id against the PeakWalk RT table built in 00d, arrange() before [1] for determinism.
  {
    .st <- modified_plot$plot$labels$subtitle
    .tag <- modified_plot$plot_tag
    if (!is.null(.st) && !is.na(.st) && grepl("RT = [0-9.]+ min", .st) &&
        is.character(.tag) && length(.tag) == 1L) {
      .id  <- sub(".*?(CP\\d+).*", "\\1", .tag)
      .smp <- sub(".*Sample:\\s*([^ |]+).*", "\\1", .st)
      .src <- if (startsWith(.tag, "C_")) "cadaver_rt_long" else "tumor_rt_long"
      if (grepl("^CP\\d+$", .id) && nzchar(.smp) && exists(.src, envir = .GlobalEnv)) {
        .rt <- get(.src, envir = .GlobalEnv) |>
          dplyr::filter(grepl(paste0("^", .id, "_"), id_subid), file == .smp) |>
          dplyr::arrange(id_subid) |> dplyr::pull(rt)
        if (length(.rt) && !is.na(.rt[1])) {
          .expected <- sprintf("RT = %.2f min", .rt[1])
          if (!grepl(.expected, .st, fixed = TRUE)) {
            stop(sprintf("vp(): subtitle RT for %s does not match the measured RT (%s). Subtitle: %s. The generator or a vp helper is printing a window midpoint again -- fix it at the source (process_single_compound.R / subtitle_rt_of()), do not patch here.", .tag, .expected, .st))
          }
        }
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
  #- Step 7.6: Italicise chemical locants in the title
#! rtx bakes title = short_name and renders it with element_text, which parses no markup, so locants
#! came out roman while the caption above the same plot italicised them. Applied here, not in rtx,
#! because the title is still a live ggplot label until ggplotGrob() runs below. Alphabetic locants
#! only -- numeric (4,4') and Greek (lambda-) locants stay roman, per the manuscript convention.
#! IDEMPOTENT. Script 09 calls vp() twice on the same plot whenever title_add is used (e.g. F3_S2_CP3017_RF <- vp(F3_S2_CP3017_R, title_add = ...)). Without the guard the second pass ran the escape below over markup the first pass added, turning *o*-Toluidine into \*o\*-Toluidine, which ggtext renders as literal asterisks. Four titles in the supplement carried that defect. The escape must still run on a virgin title so that any literal asterisk in short_name is not read as markup, so the guard tests for markup rather than skipping the block outright. element_markdown is applied on every pass because the theme is what makes the first pass's markup render.
  if (!is.null(modified_plot$plot$labels$title)) {
    t <- modified_plot$plot$labels$title
    if (!grepl("\\*[A-Za-z]+\\*", t)) {
      t <- gsub("*", "\\*", t, fixed = TRUE)
      t <- gsub("(?<![A-Za-z0-9])(sec|tert|cis|trans)-", "*\\1*-", t, perl = TRUE)
      t <- gsub("(?<![A-Za-z0-9*])([NOSomnp])-(?=[A-Za-z])", "*\\1*-", t, perl = TRUE)
      t <- gsub("\\[([a-z]),([a-z])\\]", "[*\\1*,*\\2*]", t, perl = TRUE)
      t <- gsub("\\[([a-z])\\]", "[*\\1*]", t, perl = TRUE)
      modified_plot$plot$labels$title <- t
    }
    modified_plot$plot <- modified_plot$plot +
      ggplot2::theme(plot.title = ggtext::element_markdown(hjust = 0.5, face = "bold", size = 9,
                                                           margin = margin(0, 0, 2, 0)))
  }

  #- Write final output (with _F suffix if fragment adjustment was used), and the grob, in a child R process
#! Both the PNG (ggsave inside write_small) and the gtable (ggplotGrob) BUILD the plot, and under ggplot2 4.0 every build of a plot that came in from RDS leaves ~165 MB reachable in the building process's heap for good (measured 2026-09-10: linear over 60 builds; rm/gc/dev.off do not release it). Script 09 builds 83 of them, which put ~13 GB into swap and turned scripts 09-20 into a two-hour crawl. So the build happens in a short-lived child that writes exactly the files this block always wrote and exits; the parent keeps the unbuilt ggplot, which is what the caller uses. A callr subprocess, NOT a fork: ragg's text rendering goes through CoreText, and macOS kills a forked child the moment it initialises a Cocoa class (OBJC_DISABLE_INITIALIZE_FORK_SAFETY does not rescue it -- tested). A fresh R process renders the PNG byte-identically to the parent. A child error is re-raised here so the run still fails loudly.
  suffix <- if (!is.null(mz_fragment)) "_F" else ""
  .built <- tryCatch(callr::r(function(mp, subfolder, suffix, save_grob, wd) {
    setwd(wd)
    suppressPackageStartupMessages({ library(ggplot2); library(ggtext) })
    source("R/Utilities/Validation/VP/write_small.R")
    write_small(mp, subfolder = subfolder, suffix = suffix)
    if (save_grob) {
      grob_output_dir <- file.path("Outputs/Validation", subfolder, "grobs")
      dir.create(grob_output_dir, recursive = TRUE, showWarnings = FALSE)
      grob_path <- file.path(grob_output_dir, paste0(mp$plot_tag, suffix, ".rds"))
      saveRDS(ggplotGrob(mp$plot), grob_path)
      cat(sprintf("→ Saved grob: %s\n", grob_path))
    }
    TRUE
  }, args = list(mp = modified_plot, subfolder = subfolder, suffix = suffix, save_grob = save_grob, wd = getwd()),
  libpath = .libPaths(), show = TRUE), error = function(e) conditionMessage(e))
  if (!isTRUE(.built)) stop(sprintf("vp(): rendering %s failed in the child process: %s", modified_plot$plot_tag,
                                    if (is.character(.built)) .built else "no result returned"))
  
  cat("✓ All adjustments complete\n")
  return(modified_plot)
}
