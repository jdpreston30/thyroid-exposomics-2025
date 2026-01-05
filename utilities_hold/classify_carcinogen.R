#' Classify Chemical Carcinogenicity Based on GHS and IARC Data
#'
#' Systematically categorizes chemicals into carcinogenicity risk levels using
#' a hierarchical decision tree that integrates IARC (International Agency for
#' Research on Cancer) classifications and GHS (Globally Harmonized System)
#' hazard statements.
#'
#' @param ghs Character string containing GHS hazard statement(s) with percentage
#'   confidence values in parentheses. Expected format: "H350 (75.0)", 
#'   "H350i (50.0)", or "H351 (25.0)". Can be NA, "NA", or "no ghs carcinogen statement".
#' @param iarc Character string indicating IARC carcinogen group classification:
#'   "1" (carcinogenic to humans), "2A" (probably carcinogenic), "2B" (possibly
#'   carcinogenic), "3" (not classifiable), or NA/Not classified.
#'
#' @return A character string indicating the carcinogenicity classification:
#'   \describe{
#'     \item{Known Carcinogen}{IARC Group 1}
#'     \item{Likely Carcinogen}{IARC Group 2A OR GHS H350/H350i ≥ 50\%}
#'     \item{Possible Carcinogen}{IARC Group 2B OR GHS H350/H350i 0-50\% OR H351 > 0\%}
#'     \item{Uncertain Risk}{IARC Group 3 with no supporting GHS data}
#'     \item{Unclassified}{No IARC or GHS data available}
#'   }
#'
#' @details
#' The classification follows a hierarchical priority system:
#' \enumerate{
#'   \item \strong{Highest Priority:} IARC Group 1 → "Known Carcinogen"
#'   \item \strong{High Priority:} IARC 2A OR H350/H350i ≥50\% → "Likely Carcinogen"
#'   \item \strong{Moderate Priority:} IARC 2B OR H350/H350i <50\% OR H351 >0\% → "Possible Carcinogen"
#'   \item \strong{Low Priority:} IARC 3 (no GHS support) → "Uncertain Risk"
#'   \item \strong{No Data:} No classifications available → "Unclassified"
#' }
#'
#' GHS Hazard Statements interpreted:
#' \itemize{
#'   \item H350: May cause cancer
#'   \item H350i: May cause cancer by inhalation
#'   \item H351: Suspected of causing cancer
#' }
#'
#' @examples
#' \dontrun{
#' # Known carcinogen (IARC Group 1)
#' classify_carcinogenicity(ghs = NA, iarc = "1")  # "Known Carcinogen"
#'
#' # Likely carcinogen (high GHS confidence)
#' classify_carcinogenicity(ghs = "H350 (75.0)", iarc = "3")  # "Likely Carcinogen"
#'
#' # Possible carcinogen (IARC 2B)
#' classify_carcinogenicity(ghs = NA, iarc = "2B")  # "Possible Carcinogen"
#'
#' # Unclassified
#' classify_carcinogenicity(ghs = NA, iarc = NA)  # "Unclassified"
#' }
#'
#' @export
classify_carcinogenicity <- function(ghs, iarc) {
  # Treat explicit "NA" as true NA for consistency
  if (is.na(ghs) || ghs == "NA" || ghs == "no ghs carcinogen statement") {
    ghs <- NA
  }
  if (is.na(iarc) || iarc == "NA" || iarc == "Not classified") {
    iarc <- NA
  }
  # Extract numeric values from H350, H350i, H351
  h350_perc <- as.numeric(str_extract(tolower(ghs), "(?<=h350 \\()\\d+\\.?\\d*"))
  h350i_perc <- as.numeric(str_extract(tolower(ghs), "(?<=h350i \\()\\d+\\.?\\d*"))
  h351_perc <- as.numeric(str_extract(tolower(ghs), "(?<=h351 \\()\\d+\\.?\\d*"))
  # Known Carcinogen (Highest Priority)
  if (!is.na(iarc) && iarc == "1") {
    return("Known Carcinogen")
  }
  # Likely Carcinogen (High Priority)
  if (!is.na(iarc) && iarc == "2A") {
    return("Likely Carcinogen")
  }
  if (!is.na(h350_perc) && h350_perc >= 50) {
    return("Likely Carcinogen")
  }
  if (!is.na(h350i_perc) && h350i_perc >= 50) {
    return("Likely Carcinogen")
  }
  # Possible Carcinogen (Moderate Priority)
  if (!is.na(iarc) && iarc == "2B") {
    return("Possible Carcinogen")
  }
  if (!is.na(h350_perc) && h350_perc > 0 && h350_perc < 50) {
    return("Possible Carcinogen")
  }
  if (!is.na(h350i_perc) && h350i_perc > 0 && h350i_perc < 50) {
    return("Possible Carcinogen")
  }
  if (!is.na(h351_perc) && h351_perc > 0) {
    return("Possible Carcinogen")
  }
  # Uncertain Risk (Low Priority)
  if (!is.na(iarc) && iarc == "3" && is.na(h350_perc) && is.na(h350i_perc) && is.na(h351_perc)) {
    return("Uncertain Risk")
  }
  # Unclassified (No IARC and No GHS)
  if (is.na(iarc) && is.na(ghs)) {
    return("Unclassified")
  }
  # Catch-all for unexpected cases
  return("Uncertain Risk")
}