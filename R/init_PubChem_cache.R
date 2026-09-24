#' Initialise a PubChem cache
#'
#' @param path chr. A file path to store the database.
#' Default NULL will set the database location to "./Documents/PubChem_DB"
#'
#' @returns An RDS environment
#'
#' @import storr
#' @export
#'
#' @examples
#' cache <- init_PubChem_cache()    # Initialises cache in default location
#'
#' filepath <- file.path(Sys.getenv("HOME"), "cache")
#' DB <- init_PubChem_cache(filepath)     # Initialises database in "./cache"
init_PubChem_cache <- function(path = NULL) {
  # Check `storr` is installed
  if(!requireNamespace("storr")) {
    stop("Package `storr` is not loaded!")
  }
  # Set the cache path if no path is provided
  if(is.null(path)) {
    user_home <- if (.Platform$OS.type == "windows") {
      Sys.getenv("USERPROFILE")
    } else {
      Sys.getenv("HOME")
    }
    path <- file.path(user_home, "Documents", "PubChem_DB")
  }
  # Initialise the RDS cache environment
  cache <- storr::storr_rds(path)
  # Return the RDS environment
  cache
}
