#' Standardise CAS, LIPID-MAPS and HMDB IDs
#' @description
#' Standardises identifiers (CAS, HMDB, LIPID-MAPS) to the PubChem CID, InChIKey
#' and title. Returns a single row dataframe with these as outputs.
#' Requires internet to query the PubChem database.
#'
#' HMDB identifiers must be in the format "HMDBxxxxxxx" with seven digits.
#' A helper function, `add_zeroes_to_HMDB` is available to coerce HMDB IDs to this format.
#'
#'
#' @param identifier a database identifier to search for (e.g. "58-18-4", "LMFA01100004" or "HMDB0002100")
#' @param type the database you are searching. Either "CAS", "LMP", or "HMDB"
#' @param cache an RDS cache environment. Generated with init_PubChem_cache.
#'
#' @returns a one-row dataframe with the PubChem compound ID, InChIKey, and PubChem title
#'
#' @import curl
#' @import storr
#'
#' @export
#'
#' @examples
#' # Set up a cache in the default location
#' init_PubChem_cache()
#'
#' # Search for a CAS identifier
#' cas_id <- identifier_standardiser(
#'   identifier = "58-18-4",
#'   type = "CAS",
#'   cache = cache)
#'
#' # Search for a LIPID-MAPS identifier
#' lipidmap_id <- identifier_standardiser(
#'   identifier = "LMFA01100004",
#'   type = "LMP",
#'   cache = cache)
#'
#' # Search for an HMDB identifier
#' hmdb_id <- identifier_standardiser(
#'   identifier = "HMDB0002100",
#'   type = "HMDB",
#'   cache = cache)
identifier_standardiser <- function(
    identifier,
    type = c("CAS", "LMP", "HMDB"),
    cache) {
  # Check that there is an internet connection. Stop if there isn't
  if(!curl::has_internet()) {
    stop("No internet connection available. Stopping PubChem query")
  }

  # Check that `type` is one of the values I'm looking for
  type <- match.arg(type)

  # Make sure there is a correctly connected cache
  if(is.null(cache) || (! "storr" %in% class(cache))) {
    stop(paste0('\n No cache found. Set up a cache using `init_PubChem_cache()'))
  }

  # This tryCatch checks the local cache first, then if it errors (i.e. no match)
  # it queries PubChem
  result <- tryCatch({
    # Check the local cache for the result
    cache$get(identifier)

  }, error = function(e) {
    # Rephrase `type` into what PubChem wants to see, putting the result in a new
    # variable: `database`
    if (type == "CAS") {
      database <- "CAS"
    } else if (type == "LMP") {
      database <- "Lipid+Maps+ID+(LM_ID)"
    } else if (type == "HMDB") {
      database <- "HMDB+ID"
    } else {
      stop("I don't know how `match.arg` missed this but `type` must be one of 'CAS',
        'LMP' or 'HMDB'")
    }

    # tryCatch runs HTTP request and saves "NA" if it fails
    query <- tryCatch({

      # Send a programmatic request to PubChem using the PUG REST API. This query
      # reruns the compound ID (CID), InChIKey and Title.
      resp <-
        httr2::request(
          paste0(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/identifier/",
            identifier,
            "/property/InChIKey,Title/JSON?identifier_type=",
            database)) |>
        # PubChem asks applications to make no more than 5 requests per second
        httr2::req_throttle(
          capacity = 5,
          fill_time_s = 1) |>
        httr2::req_perform()

      # Save the data to a more usable format
      query_df <- as.data.frame(httr2::resp_body_json(resp)[[1]][[1]][[1]])
      query_df[[type]] <- identifier
      query_df$PubChem_Message <- NA

      query_df

    }, error = function(e) {

      # If there was an error, print the error message ...
      message("PubChem lookup failed: ", e$message)

      # ... save an identical dataframe albeit with NA instead of results ...
      query_df <- data.frame(
        CID = c(NA_integer_),
        InChIKey = c(NA_character_),
        DB_Name = c(NA_character_))
      query_df[[type]] <- identifier
      query_df$PubChem_Message <- e$message

      query_df
    })
    # Save the result to the local cache
    cache$set(identifier, query)

    # Return query
    query
  })
}
