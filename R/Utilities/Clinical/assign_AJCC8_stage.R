#' Assign AJCC 8th Edition Overall Stage for Differentiated Thyroid Cancer
#'
#' Calculates Stage I–IV for papillary and follicular thyroid carcinoma
#' based on AJCC/UICC TNM 8th edition criteria, with an option to ignore age.
#'
#' @param df Data frame with columns for T, N, M, and Age (in years).
#' @param t_col Character, column name containing T stage (e.g. "T_stage_comp").
#' @param n_col Character, column name containing N stage (e.g. "N").
#' @param m_col Character, column name containing M stage (e.g. "M").
#' @param age_col Character, column name with patient age (default = "age").
#' @param out_col Character, name of output stage column (default = "AJCC8_Stage").
#' @param consider_age Logical. If TRUE (default), apply AJCC age-specific rules (<55 vs ≥55).
#'   If FALSE, all patients are staged using ≥55-year criteria.
#'
#' @return Input data frame with a new ordered factor column for AJCC Stage.
#' @export
assign_AJCC8_stage <- function(
    df,
    t_col = "T_stage_comp",
    n_col = "N",
    m_col = "M",
    age_col = "age",
    out_col = "AJCC8_Stage",
    consider_age = TRUE) {
  T_stage <- df[[t_col]]
  N_stage <- df[[n_col]]
  M_stage <- df[[m_col]]
  Age <- df[[age_col]]

  if (!consider_age) {
    # Treat everyone as ≥55
    Age <- rep(60, length(Age))
  }

  df[[out_col]] <- dplyr::case_when(
    # ---- Age <55 ----
    Age < 55 & M_stage == "M0" ~ "I",
    Age < 55 & M_stage == "M1" ~ "II",

    # ---- Age ≥55 ----
    Age >= 55 & M_stage == "M1" ~ "IV",
    Age >= 55 & T_stage %in% c("T4", "T4a", "T4b") & M_stage == "M0" ~ "III",
    Age >= 55 & T_stage %in% c("T1", "T2") & N_stage %in% c("N1", "N1a", "N1b") & M_stage == "M0" ~ "II",
    Age >= 55 & T_stage %in% c("T3") & M_stage == "M0" ~ "II",
    Age >= 55 & T_stage %in% c("T1", "T2") &
      (N_stage %in% c("N0", "NX") | is.na(N_stage)) & M_stage == "M0" ~ "I",
    TRUE ~ NA_character_
  )

  df[[out_col]] <- factor(df[[out_col]],
    levels = c("I", "II", "III", "IV"),
    ordered = TRUE
  )

  return(df)
}
