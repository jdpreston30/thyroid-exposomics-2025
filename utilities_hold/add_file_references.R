#' Add File References to Validation Table
#'
#' @description
#' Joins file references from file_list to check columns (check1-check5) and
#' relocates file columns next to their corresponding check columns.
#'
#' @param data Tibble with check1-check5 columns containing sample IDs
#' @param file_list Tibble with ID and files columns
#'
#' @return Tibble with file1-file5 columns added and relocated
#'
#' @examples
#' tumor_top5_iarc |> add_file_references(file_list)
add_file_references <- function(data, file_list) {
  data |>
    left_join(file_list |> select(ID, file1 = files), by = c("check1" = "ID")) |>
    left_join(file_list |> select(ID, file2 = files), by = c("check2" = "ID")) |>
    left_join(file_list |> select(ID, file3 = files), by = c("check3" = "ID")) |>
    left_join(file_list |> select(ID, file4 = files), by = c("check4" = "ID")) |>
    left_join(file_list |> select(ID, file5 = files), by = c("check5" = "ID")) |>
    relocate(file1, .after = check1) |>
    relocate(file2, .after = check2) |>
    relocate(file3, .after = check3) |>
    relocate(file4, .after = check4) |>
    relocate(file5, .after = check5)
}
