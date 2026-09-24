#' Decide which IDs to query
#' @description
#' This is a helper function, which other functions use to work out which
#' database identifier to query.
#'
#' @param data a dataframe
#' @param col chr. The name of a column in `data` which contains database
#' identifiers
#'
#' @returns a logical vector
.get_logic <- function(data, col) {
  result <- !is.na(data[col])
}
