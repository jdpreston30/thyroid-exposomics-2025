#' Create Balloon Plot for Variant Analysis
#'
#' Generates a balloon plot displaying count data across thyroid cancer variants
#' and categories. Bubble size represents count magnitude, and colors correspond
#' to variant-specific palette.
#'
#' @param data A data frame containing the balloon plot data with columns:
#'   \code{Variant} (character), \code{usage_class} or category column (character),
#'   and \code{n} (numeric count).
#' @param y_var Character string specifying the y-axis variable (e.g., "usage_class").
#'   Default is "usage_class".
#' @param size_range Numeric vector of length 2 specifying min and max bubble sizes.
#'   Default is c(0.9, 4.5).
#' @param size_breaks Numeric vector specifying breaks for size legend.
#'   Default is c(1, 3, 5, 7).
#' @param show_x_labels Logical indicating whether to display x-axis variant labels.
#'   Default is FALSE.
#' @param y_text_size Numeric specifying y-axis text size. Default is 8.
#'
#' @return A ggplot2 object representing the balloon plot.
#'
#' @details
#' The balloon plot uses the following visual encoding:
#' \itemize{
#'   \item \strong{Bubble size:} Magnitude of count (n)
#'   \item \strong{Bubble fill:} Variant-specific colors from \code{variant_colors} palette
#'   \item \strong{Bubble outline:} Black border for visual clarity
#'   \item \strong{Background:} Minimal theme with subtle grid lines
#' }
#'
#' Visual specifications:
#' \itemize{
#'   \item Font: Arial bold at various sizes
#'   \item Y-axis labels: 8pt bold (default)
#'   \item Legend text: 10pt bold
#'   \item Grid: Major (grey80, 0.3 linewidth) and minor (grey90, 0.2 linewidth)
#'   \item Size legend: White fill with black outline (stroke = 0.5)
#' }
#'
#' @examples
#' \dontrun{
#' # Basic balloon plot with default settings
#' balloon_data <- data.frame(
#'   Variant = rep(c("Follicular", "FV-PTC", "Papillary"), each = 5),
#'   usage_class = rep(c("Class A", "Class B", "Class C", "Class D", "Class E"), 3),
#'   n = sample(1:10, 15, replace = TRUE)
#' )
#' plot_balloon(balloon_data)
#'
#' # Custom size range and breaks
#' plot_balloon(balloon_data, size_range = c(1, 6), size_breaks = c(2, 4, 6, 8))
#'
#' # Show x-axis labels
#' plot_balloon(balloon_data, show_x_labels = TRUE)
#'
#' # Custom y-axis category variable
#' balloon_data2 <- data.frame(
#'   Variant = rep(c("Follicular", "FV-PTC", "Papillary"), each = 3),
#'   chemical_class = rep(c("Pesticide", "PAH", "Phthalate"), 3),
#'   n = sample(1:8, 9, replace = TRUE)
#' )
#' plot_balloon(balloon_data2, y_var = "chemical_class")
#' }
#'
#' @seealso \code{\link{variant_colors}} for color palette,
#'   \code{\link{plot_detection_scatter}} for related visualization functions
#'
#' @export
plot_balloon <- function(data,
                         y_var = "usage_class",
                         size_range = c(0.9, 4.5),
                         size_breaks = c(1, 2, 3),
                         show_x_labels = FALSE,
                         y_text_size = 9) {
  # Validate inputs
  required_cols <- c("Variant", y_var, "n")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Create base plot
  p <- ggplot(data, aes(x = Variant, y = .data[[y_var]], size = n, fill = Variant)) +
    geom_point(shape = 21, color = "black") +
    scale_size_continuous(
      range = size_range,
      limits = c(0, max(data$n, na.rm = TRUE)),
      breaks = size_breaks,
      labels = as.character(size_breaks)
    ) +
    scale_fill_manual(values = variant_colors, guide = "none") +
    labs(size = "Count") +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial", face = "bold", size = 14, color = "black"),
      axis.text.y = element_text(family = "Arial", face = "bold", size = y_text_size, color = "black"),
      legend.text = element_text(family = "Arial", face = "bold", size = 8, color = "black"),
      legend.title = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.3, linetype = "solid"),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2, linetype = "solid")
    ) +
    guides(size = guide_legend(override.aes = list(fill = "white", color = "black", stroke = 0.5)))

  # Handle x-axis labels
  if (show_x_labels) {
    p <- p + theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        family = "Arial",
        face = "bold",
        size = 10,
        color = "black"
      )
    )
  } else {
    p <- p + theme(axis.text.x = element_blank())
  }

  return(p)
}
