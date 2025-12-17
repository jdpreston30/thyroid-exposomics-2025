#' Plot Stacked Bar Chart of Carcinogen Classifications by Variant
#'
#' Creates a stacked bar plot showing the count of chemicals in each carcinogen
#' classification category for each thyroid cancer variant.
#'
#' @param data A data frame with columns:
#'   \describe{
#'     \item{Variant}{Character indicating variant type (Follicular, FV-PTC, Papillary)}
#'     \item{Known Carcinogen}{Numeric count of known carcinogens}
#'     \item{Likely Carcinogen}{Numeric count of likely carcinogens}
#'     \item{Possible Carcinogen}{Numeric count of possible carcinogens}
#'     \item{Uncertain Risk}{Numeric count of uncertain risk chemicals (optional)}
#'   }
#'
#' @return A ggplot object showing a stacked bar plot with carcinogen categories
#'   colored according to carcinogen_colors palette.
#'
#' @details
#' The plot includes:
#' \itemize{
#'   \item X-axis: Variants in order (Follicular, FV-PTC, Papillary), rotated 45 degrees
#'   \item Y-axis: Count of chemicals highest in variant
#'   \item Stacked bars with carcinogen classification colors
#'   \item Black borders on bars (0.4 linewidth)
#'   \item Horizontal legend at top left
#'   \item Arial font, size 12 bold for axis titles, size 10 bold for axis text
#' }
#'
#' @examples
#' \dontrun{
#' # Prepare data
#' carc_by_variant <- MT_final %>%
#'   filter(!is.na(highest)) %>%
#'   group_by(Variant, carcinogen_status) %>%
#'   summarise(n = n()) %>%
#'   pivot_wider(names_from = carcinogen_status, values_from = n, values_fill = 0)
#' 
#' # Create plot
#' p <- plot_carcinogen_stacked(carc_by_variant)
#' ggsave("carcinogen_by_variant.pdf", p, width = 4, height = 5)
#' }
#'
#' @seealso \code{\link{carcinogen_colors}} for color palette
#'
#' @export
plot_carcinogen_stacked <- function(data) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # Convert data to long format for stacking
  plot_data <- data |>
    pivot_longer(
      cols = -Variant,
      names_to = "carcinogen_status",
      values_to = "count"
    ) |>
    mutate(
      Variant = factor(Variant, levels = c("Follicular", "FV-PTC", "Papillary")),
      carcinogen_status = factor(carcinogen_status, 
                                 levels = c("Known Carcinogen", "Likely Carcinogen", 
                                           "Possible Carcinogen", "Uncertain Risk"))
    )
  
  # Create stacked bar plot
  p <- ggplot(plot_data, aes(x = Variant, y = count, fill = carcinogen_status)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.6, width = 0.6) +
    scale_fill_manual(
      values = carcinogen_colors,
      breaks = c("Known Carcinogen", "Likely Carcinogen", 
                "Possible Carcinogen", "Uncertain Risk")
    ) +
    guides(fill = guide_legend(
      keywidth = unit(0.25, "cm"),
      keyheight = unit(0.25, "cm"),
      byrow = FALSE,
      label.position = "right",
      ncol = 1
    )) +
    scale_x_discrete(
      labels = c("Follicular" = "Follicular\n", "FV-PTC" = "FV-PTC\n", "Papillary" = "Papillary\n")
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05), add = 0)
    ) +
    labs(
      x = NULL,
      y = "# Highest in Variant",
      title = " ",
      fill = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.title = element_text(face = "bold", color = "transparent", size = 12, hjust = 0.5),
      legend.position = c(0.0, 1.05),
      legend.justification = c(0, 1),
      legend.direction = "horizontal",
      legend.text = element_text(size = 8, face = "plain", family = "Arial"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "white", color = "black", linewidth = 0.05),
      legend.key.size = unit(0.25, "cm"),
      legend.spacing.x = unit(0.02, "cm"),
      legend.key.spacing.y = unit(0.07, "cm"),
      legend.box.spacing = unit(0, "cm"),
      axis.text.x = element_text(face = "bold", color = "black", size = 10),
      axis.text.y = element_text(face = "bold", color = "black", size = 10),
      axis.title.y = element_text(face = "bold", color = "black", size = 12, margin = margin(r = 10)),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8)
    )
  
  return(p)
}
