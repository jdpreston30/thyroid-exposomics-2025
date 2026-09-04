#' Plot Top 5 Quantitative Compounds by Variant
#'
#' Creates a faceted bar plot showing peak area z-scores for the top 5 
#' quantitative compounds across thyroid cancer variants
#'
#' @param data Tibble with variant column and name_sub_lib_id columns (z-scores)
#' @param compound_names Named vector mapping name_sub_lib_id to short_name (optional)
#' @param return_legend Logical, if TRUE returns list with plot and legend grobs (default FALSE)
#' @param add_cld Logical, if TRUE adds compact letter display from posthoc_table_pvalues (default FALSE)
#' @return ggplot object or list with plot and legend if return_legend = TRUE
#' @export
#' 
plot_top5_quant <- function(data, compound_names = NULL, return_legend = FALSE, add_cld = FALSE,
                            axis_text_size = 10, axis_title_size = 12) {
  
  # Get short names from MT_final if not provided
  if (is.null(compound_names)) {
    if (exists("MT_final", envir = .GlobalEnv)) {
      mt <- get("MT_final", envir = .GlobalEnv)
      compound_names <- setNames(mt$short_name, mt$name_sub_lib_id)
    } else {
      stop("compound_names must be provided or MT_final must exist in global environment")
    }
  }
  
  # Remap raw "FV-PTC" to display label before any factoring
  data <- data |> mutate(variant = if_else(variant == "FV-PTC", "IEFVPTC", variant))
  # Get column order from original data (preserves order)
  col_order <- colnames(data)[colnames(data) != "variant"]
  
  # Map to short names for ordering
  chemical_order <- sapply(col_order, function(x) {
    chem <- compound_names[x]
    # Apply same overrides
    if (!is.na(chem) && chem == "DNOP") return("Di-n-octyl phthalate")
    if (!is.na(chem) && chem == "MEHP") return("Mono(2-ethylhexyl) phthalate")
    return(chem)
  }, USE.NAMES = FALSE)
  
  # Pivot to long format
  data_long <- data %>%
    pivot_longer(
      cols = -variant,
      names_to = "name_sub_lib_id",
      values_to = "z_score"
    ) %>%
    mutate(
      # Map to short names
      chemical = compound_names[name_sub_lib_id],
      # Override specific chemical names for display
      chemical = case_when(
        chemical == "DNOP" ~ "Di-n-octyl phthalate",
        chemical == "MEHP" ~ "Mono(2-ethylhexyl) phthalate",
        TRUE ~ chemical
      ),
      # Factor chemical with original column order
      chemical = factor(chemical, levels = chemical_order),
      # Factor variant with correct order (reversed for y-axis top-to-bottom display)
      variant = factor(variant,
                      levels = c("Papillary", "IEFVPTC", "Follicular"),
                      labels = c("Papillary", "IEFVPTC", "Follicular"))
    )
  
  # Get variant colors from theme
  variant_colors <- c(
    "Follicular" = "#294B88",
    "IEFVPTC" = "#23744E",
    "Papillary" = "#DF8D09"
  )
  
  # Create alpha version for fills
  variant_colors_fill <- alpha(variant_colors, 0.5)
  
  # Get CLD letters if requested
  cld_data <- NULL
  if (add_cld && exists("posthoc_table_pvalues", envir = .GlobalEnv)) {
    posthoc <- get("posthoc_table_pvalues", envir = .GlobalEnv)
    
    # Get unique compound names and extract CLD letters
    cld_data <- posthoc %>%
      select(Name, FTC_let, FV_PTC_let, PTC_let) %>%
      distinct() %>%
      pivot_longer(cols = c(FTC_let, FV_PTC_let, PTC_let), 
                   names_to = "variant_col", values_to = "cld") %>%
      mutate(
        variant = case_when(
          variant_col == "FTC_let" ~ "Follicular",
          variant_col == "FV_PTC_let" ~ "IEFVPTC",
          variant_col == "PTC_let" ~ "Papillary"
        ),
        # Remove superscripts and leading values from CLD
        cld = gsub("^-?[0-9.]+", "", cld),
        # Map chemical names
        chemical = case_when(
          Name == "DNOP" ~ "Di-n-octyl phthalate",
          Name == "MEHP" ~ "Mono(2-ethylhexyl) phthalate",
          TRUE ~ Name
        ),
        # Factor chemical with same order as main plot
        chemical = factor(chemical, levels = chemical_order),
        # Factor variant with same order as main plot
        variant = factor(variant,
                        levels = c("Papillary", "IEFVPTC", "Follicular"),
                        labels = c("Papillary", "IEFVPTC", "Follicular"))
      ) %>%
      filter(chemical %in% chemical_order)
  }
  
  # Create plot (flipped coordinates)
  p <- ggplot(data_long, aes(x = z_score, y = variant, fill = variant, color = variant)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, linewidth = 0.5, show.legend = TRUE) +
#! Seeded, so 3A's dots land identically on every draw -- unseeded jitter re-randomises per DRAW, not per build, so PNG/TIFF/PDF of one figure object disagreed. Matches seed = 42 in plot_iarc and plot_detection_scatter. geom_point, since geom_jitter() errors if given both position and width/height.
    geom_point(position = position_jitter(height = 0.15, width = 0, seed = 42),
               size = 0.5, alpha = 1, show.legend = FALSE) +
#! Italics applied at the labeller, NOT by re-levelling `chemical`: chemical_order and the compact-letter
#! lookup both key off those levels, so they must stay plain. label_parsed renders the plotmath strings.
    facet_wrap(~chemical, ncol = 1, scales = "fixed", strip.position = "left", drop = FALSE,
               labeller = ggplot2::as_labeller(locant_plotmath, default = ggplot2::label_parsed)) +
    scale_fill_manual(
      values = variant_colors_fill,
      breaks = c("Follicular", "IEFVPTC", "Papillary"),
      labels = c("Follicular ", "IEFVPTC   ", "Papillary"),
      guide = guide_legend(nrow = 1, byrow = TRUE, label.hjust = 0)
    ) +
    scale_color_manual(
      values = variant_colors,
      breaks = c("Follicular", "IEFVPTC", "Papillary"),
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = seq(-2, 2, by = 1),
      limits = c(-2.3, 3)
    ) +
    labs(
      #! z-score is a variable symbol, so lowercase italic z per AMA/APA and Methods 2.6; expression() rather than ggtext to avoid adding a theme dependency for one label
      x = expression(paste(bold("Peak Area "), bolditalic(z), bold("-score"))),
      y = NULL,
      fill = NULL
    ) +
    # Add CLD letters if requested
    {if (!is.null(cld_data)) {
      geom_text(data = cld_data, aes(x = 2.95, y = variant, label = cld),
                inherit.aes = FALSE, size = 3, fontface = "italic", family = "Arial")
    }} +
    guides(
      fill = guide_legend(
        keywidth = unit(0.4, "cm"),
        keyheight = unit(0.25, "cm"),
        byrow = TRUE,
        label.position = "right",
        label.hjust = 0,
        nrow = 1,
        override.aes = list(color = variant_colors, alpha = 0.5)
      ),
      color = "none"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 10, base_family = "Arial") +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 9),
      strip.placement = "outside",
#! 10 to match 3B/3D/3E/3F's tick labels. Was 8 -- the smallest text in the figure, on the widest
#! panel, which inverted the relationship between available space and type size.
      axis.text.x = element_text(face = "bold", size = axis_text_size),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
#! .x specifically, because axis.ticks.y is blanked above -- a bare axis.ticks would depend on order.
#! 0.8 matches plot_iarc (3E/3F) and plot_carcinogen_stacked (3D). Previously unset, so it inherited
#! theme_classic's base_line_size = base_size/22 = 0.4545, which rendered 4 px vs their 7 px at 300 dpi.
      axis.ticks.x = element_line(color = "black", linewidth = 0.8),
#! Absolute 0.15 cm to match 3D/3E/3F. theme_classic stores this as rel(0.5) against base_size, so at
#! base_size = 10 it resolved to 2.5 pt (10.4 px at 300 dpi) versus their 17.7 px -- 59% as long.
      axis.ticks.length = unit(0.15, "cm"),
#! 12 to match 3D/3E/3F's axis titles. Sits directly above 3B/3C, so check clearance at y ~ 6.56.
      axis.title.x = element_text(face = "bold", size = axis_title_size),
      axis.line = element_blank(),
      legend.position = "top",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.box.margin = margin(l = 5/300, unit = "in"),
      legend.text = element_text(size = 8, face = "bold", family = "Arial"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "white", linewidth = 0.5),
      legend.key.size = unit(0.25, "cm"),
      legend.spacing.x = unit(0.15, "cm"),
      legend.spacing.y = unit(0.02, "cm"),
      legend.box.spacing = unit(0, "cm"),
#! Back to 0 so adjacent panel.border edges coincide and each internal divider renders as ONE line at
#! the border linewidth, not a separated pair.
      panel.spacing = unit(0, "lines"),
#! 0.379 = 9 px at 800 dpi, measured, matching tile_linewidth in plot_qualitative_heatmap() so 3A's facet
#! dividers and 3B's cell dividers read as the same weight. Was 0.337 (8 px). thickness_in = linewidth *
#! 2.845276 / 96 holds here because adjacent panel.border edges coincide, so the full stroke is drawn.
#! panel.border is per-panel and cannot differ between the outer edge and the facet dividers, so it is
#! set to the DIVIDER weight here; the heavier outer frame comes from add_panel_frame() in 14.
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.379)
    )
  
  # If return_legend is TRUE, extract legend and return both
  if (return_legend) {
    # Create a temporary plot with legend visible to extract it
    p_temp <- p + theme(legend.position = "right")
    legend <- cowplot::get_legend(p_temp)
    
    # Remove legend from main plot
    p_no_legend <- p + theme(legend.position = "none")
    
    return(list(plot = p_no_legend, legend = legend))
  }
  
  return(p)
}
