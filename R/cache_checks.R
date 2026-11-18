#' Check your cache exists and is the right type
#'
#' @param cache a putative `storr` cache object
#'
#' @returns a correctly set up PubChem cache
#' @export
#'
.cache_checks <- function(cache) {

  if(is.null(cache)) {

  # Initialise the cache if one is not specified
    cache <- init_PubChem_cache()
    message("No cache specified. Initialising a cache in the default location.")
  } else if(! "storr" %in% class(cache)) {

  # Initialise a cache if the cache object is the wrong type
    cache <- init_PubChem_cache()
    message("Incorrect cache object type. Initialising a cache in the default location.")
  } else {cache <- cache}

  cache

}
