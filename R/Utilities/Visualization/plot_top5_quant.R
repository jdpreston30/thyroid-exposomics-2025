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
plot_top5_quant <- function(data, compound_names = NULL, return_legend = FALSE, add_cld = FALSE) {
  
  # Get short names from MT_final if not provided
  if (is.null(compound_names)) {
    if (exists("MT_final", envir = .GlobalEnv)) {
      mt <- get("MT_final", envir = .GlobalEnv)
      compound_names <- setNames(mt$short_name, mt$name_sub_lib_id)
    } else {
      stop("compound_names must be provided or MT_final must exist in global environment")
    }
  }
  
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
      # Factor variant with correct order (reversed for y-axis top-to-bottom display)
      variant = factor(variant,
                      levels = c("Papillary", "FV-PTC", "Follicular"),
                      labels = c("Papillary", "FV-PTC", "Follicular"))
    )
  
  # Get variant colors from theme
  variant_colors <- c(
    "Follicular" = "#294B88",
    "FV-PTC" = "#23744E",
    "Papillary" = "#DF8D09"
  )
  
  # Create alpha version for fills
  variant_colors_fill <- alpha(variant_colors, 0.5)
  
  # Order chemicals by their appearance in the original data columns (bottom-to-top)
  chemical_order <- unique(data_long$chemical)
  data_long <- data_long %>%
    mutate(chemical = factor(chemical, levels = chemical_order))
  
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
          variant_col == "FV_PTC_let" ~ "FV-PTC",
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
        # Factor variant with same order as main plot
        variant = factor(variant,
                        levels = c("Papillary", "FV-PTC", "Follicular"),
                        labels = c("Papillary", "FV-PTC", "Follicular"))
      ) %>%
      filter(chemical %in% chemical_order)
  }
  
  # Create plot (flipped coordinates)
  p <- ggplot(data_long, aes(x = z_score, y = variant, fill = variant, color = variant)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, linewidth = 0.5, show.legend = TRUE) +
    geom_jitter(height = 0.15, width = 0, size = 0.5, alpha = 1, show.legend = FALSE) +
    facet_wrap(~chemical, ncol = 1, scales = "fixed", strip.position = "left") +
    scale_fill_manual(
      values = variant_colors_fill,
      breaks = c("Follicular", "FV-PTC", "Papillary"),
      labels = c("Follicular ", "FV-PTC   ", "Papillary"),
      guide = guide_legend(nrow = 1, byrow = TRUE, label.hjust = 0)
    ) +
    scale_color_manual(
      values = variant_colors,
      breaks = c("Follicular", "FV-PTC", "Papillary"),
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = seq(-2, 2, by = 1),
      limits = c(-2.3, 3)
    ) +
    labs(
      x = "Peak Area Z-Score",
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
      axis.text.x = element_text(face = "bold", size = 8),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(face = "bold", size = 10),
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
      panel.spacing = unit(0, "lines"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
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
