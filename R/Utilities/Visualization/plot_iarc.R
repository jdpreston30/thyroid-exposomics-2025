#' Plot IARC Validation Concentrations
#'
#' Creates a violin plot with overlaid points showing chemical concentrations
#' in tumor tissue vs. control (cadaver) samples on a logarithmic scale.
#'
#' @param data A data frame with the following columns:
#'   \describe{
#'     \item{tumor_vs_ctrl}{Character indicating sample type ("Thyroid Tumor Tissue" or "Non-Cancer Cadaver Thyroids")}
#'     \item{concentration}{Numeric concentration value in PPM}
#'   }
#' @param chemical_name Character string for the chemical name to display on x-axis
#' @param p_value Optional numeric p-value from statistical test (e.g., Wilcoxon).
#'   If provided, will be displayed in the top-left corner of the plot.
#'
#' @return A ggplot object showing violin plots with overlaid points, colored by
#'   tumor vs. control tissue type on a logarithmic y-axis scale.
#'
#' @details
#' The plot includes:
#' \itemize{
#'   \item Y-axis: Logarithmic scale from 10^-4 to 10^3 PPM
#'   \item X-axis: Chemical name (no axis label)
#'   \item Violin plots with 50% alpha fill
#'   \item Overlaid jittered points
#'   \item Tumor vs. control colors from tumor_noncancer_colors palette
#'   \item Optional p-value annotation in top-left corner
#'   \item Arial font, size 12 bold for axis titles, size 10 bold for axis text
#' }
#'
#' @examples
#' \dontrun{
#' # Prepare data
#' plot_data <- full_joiner %>%
#'   select(sample_ID, variant, tumor_vs_ctrl, `o-Toluidine_0_BP3.GC2_CP3017`) %>%
#'   rename(concentration = `o-Toluidine_0_BP3.GC2_CP3017`)
#' 
#' # Create plot with statistical test
#' wilcox_p <- wilcox.test(concentration ~ tumor_vs_ctrl, data = plot_data)$p.value
#' p <- plot_iarc(plot_data, chemical_name = "o-Toluidine", p_value = wilcox_p)
#' ggsave("iarc_validation.pdf", p, width = 4, height = 5)
#' }
#'
#' @seealso \code{\link{tumor_noncancer_colors}} for color palette
#'
#' @export
plot_iarc <- function(data, chemical_name, p_value = NULL) {
  library(ggplot2)
  library(dplyr)
  
  # Recode tumor_vs_ctrl to match expected values and create labels with line breaks
  plot_data <- data |>
    mutate(
      # Map "Control"/"Tumor" to full names for color matching
      tumor_vs_ctrl_full = case_when(
        tumor_vs_ctrl == "Control" ~ "Non-Cancer Cadaver Thyroids",
        tumor_vs_ctrl == "Tumor" ~ "Thyroid Tumor Tissue",
        TRUE ~ tumor_vs_ctrl
      ),
      tumor_vs_ctrl_full = factor(tumor_vs_ctrl_full, 
                                  levels = c("Non-Cancer Cadaver Thyroids", "Thyroid Tumor Tissue")),
      # Create x-axis labels with line breaks
      x_label = factor(tumor_vs_ctrl_full,
                      levels = c("Non-Cancer Cadaver Thyroids", "Thyroid Tumor Tissue"),
                      labels = c("Cadaver\nThyroid", "Thyroid\nTumor"))
    )
  
  # Calculate dynamic y-axis limits based on data range
  min_val <- min(plot_data$concentration, na.rm = TRUE)
  max_val <- max(plot_data$concentration, na.rm = TRUE)
  
  # Find appropriate log10 scale bounds (round down min, round up max)
  y_min <- 10^floor(log10(min_val))
  y_max <- 10^ceiling(log10(max_val))
  
  # Generate breaks from y_min to y_max
  y_breaks <- 10^seq(log10(y_min), log10(y_max), by = 1)
  
  # Calculate summary statistics (mean and SEM for error bars)
  summary_data <- plot_data |>
    group_by(tumor_vs_ctrl_full, x_label) |>
    summarise(
      mean_conc = mean(concentration, na.rm = TRUE),
      sem = sd(concentration, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Create bar plot with error bars and points
  p <- ggplot() +
    # Bar plot (mean) - rectangles from bottom of scale to mean
    geom_rect(
      data = summary_data,
      aes(xmin = as.numeric(x_label) - 0.35, xmax = as.numeric(x_label) + 0.35,
          ymin = y_min, ymax = mean_conc, fill = tumor_vs_ctrl_full),
      alpha = 0.5,
      color = NA
    ) +
    # Bar outlines
    geom_rect(
      data = summary_data,
      aes(xmin = as.numeric(x_label) - 0.35, xmax = as.numeric(x_label) + 0.35,
          ymin = y_min, ymax = mean_conc, color = tumor_vs_ctrl_full),
      fill = NA,
      linewidth = 0.6
    ) +
    # Individual points
    geom_point(
      data = plot_data,
      aes(x = x_label, y = concentration, color = tumor_vs_ctrl_full),
      size = 1,
      alpha = 1,
      shape = 16,
      position = position_jitter(width = 0.15, seed = 42)
    ) +
    scale_fill_manual(
      values = tumor_noncancer_colors,
      breaks = c("Thyroid Tumor Tissue", "Non-Cancer Cadaver Thyroids")
    ) +
    scale_color_manual(
      values = tumor_noncancer_colors,
      breaks = c("Thyroid Tumor Tissue", "Non-Cancer Cadaver Thyroids")
    ) +
    scale_y_log10(
      limits = c(y_min, y_max),
      breaks = y_breaks,
      labels = function(x) sprintf("%.0e", x),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = NULL,
      y = "log(PPM)",
      title = chemical_name,
      fill = NULL,
      color = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.position = "none",
      plot.title = element_text(face = "bold", color = "black", size = 12, hjust = 0.5),
      axis.text.x = element_text(face = "bold", color = "black", size = 10),
      axis.text.y = element_text(face = "bold", color = "black", size = 10),
      axis.title.y = element_text(face = "bold", color = "black", size = 12, margin = margin(r = 10)),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  # Add p-value annotation if provided
  if (!is.null(p_value)) {
    p_text <- if (p_value < 0.001) {
      "p < 0.001"
    } else {
      sprintf("p = %.3f", p_value)
    }
    
    p <- p + annotate(
      "text",
      x = -Inf,
      y = y_max,
      label = p_text,
      hjust = -0.25,
      vjust = 1.5,
      size = 8 / .pt,
      family = "Arial",
      fontface = "plain"
    )
  }
  
  return(p)
}
