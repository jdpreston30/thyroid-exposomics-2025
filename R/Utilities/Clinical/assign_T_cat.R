#' Assign AJCC 8th-edition T stage (T1–T4) for differentiated thyroid cancer
#'
#' Determines the primary tumor T category using largest tumor dimension (LD)
#' and extrathyroidal extension (ETE) per AJCC/UICC TNM 8th edition.
#'
#' @section Assumption:
#' ETE is coded as 0 = none, 1 = microscopic/minimal, 2 = gross/extensive.
#' This function assumes ETE = 2 always represents gross ETE beyond strap muscles
#' and therefore upgrades to at least T4 if size-based T-stage is lower.
#'
#' @param df Data frame containing tumor size and ETE columns.
#' @param ld_col Character. Column name with largest tumor dimension in cm (default "LD").
#' @param ete_col Character. Column name with ETE code 0/1/2 (default "ETE").
#' @param units Either "cm" (default) or "mm" for the size column.
#' @param out_col Character. Name of output column to write (default "T_stage").
#'
#' @details
#' Size-based rules (microscopic ETE ignored in AJCC8):
#' \itemize{
#'   \item T1: \eqn{\le} 2 cm, confined to thyroid
#'   \item T2: > 2 cm and \eqn{\le} 4 cm, confined to thyroid
#'   \item T3: > 4 cm, confined to thyroid
#' }
#' ETE handling:
#' \itemize{
#'   \item ETE == 0 or 1: keep size-based T (microscopic ETE does not upstage)
#'   \item ETE == 2: upgrade to T4 if size-based T-stage is lower
#' }
#' Lymphovascular invasion (LVI) and multifocality are not used for T1–T4 assignment.
#'
#' @return A copy of \code{df} with an ordered factor column \code{out_col} in \{T1,T2,T3,T4\}.
#' @export
assign_T_cat <- function(
    df,
    ld_col = "LD",
    ete_col = "ETE",
    units = c("cm", "mm"),
    out_col = "T_stage") {
  units <- match.arg(units)

  LD <- df[[ld_col]]
  ETE <- df[[ete_col]]

  # Convert to cm if needed
  LD_cm <- if (units == "mm") LD / 10 else LD

  # Size-only T (AJCC8; microscopic ETE ignored)
  T_size <- dplyr::case_when(
    is.na(LD_cm) ~ NA_character_,
    LD_cm <= 2 ~ "T1",
    LD_cm > 2 & LD_cm <= 4 ~ "T2",
    LD_cm > 4 ~ "T3",
    TRUE ~ NA_character_
  )

  # Upgrade logic for gross ETE (ETE==2 → always T4 if lower)
  ord <- c(T1 = 1L, T2 = 2L, T3 = 3L, T4 = 4L)
  T_final <- dplyr::case_when(
    is.na(T_size) ~ NA_character_,
    is.na(ETE) ~ T_size,
    ETE %in% c(0, 1) ~ T_size,
    ETE == 2 ~ {
      cur_ord <- ord[T_size]
      ifelse(cur_ord < ord["T4"], "T4", T_size)
    },
    TRUE ~ T_size
  )

  df[[out_col]] <- factor(T_final, levels = c("T1", "T2", "T3", "T4"), ordered = TRUE)
  df
}
