#' Ensure HMDB identifiers have 7 digits
#'
#' @description
#' PubChem search requires 7 digits in HMDB identifiers. This function repairs them.
#'
#' @param string A single chr, beginning with "HMDB"
#'
#' @returns A single chr, in the format "HMDBxxxxxxx"
#'
#' @import stringr
#'
#' @export
#'
#' @examples
#' add_zeroes_to_HMDB("HMDB1")
#' # [1] "HMDB0000001"
add_zeroes_to_HMDB <- function(string) {
  # Check that `string` begins with "HMDB"
  if(!str_detect(string, "^HMDB")){
    stop("That doesn't look to be an HMDB identifier")
  }
  # Add zeroes to ensure there are seven digits
  digits <- string |>
    stringr::str_extract("(\\d+)$") |>
    stringr::str_pad(width = 7, side = "left", pad = "0")

  result <- paste0("HMDB", digits)

  result
}
