#' Standardise metabolite names
#'
#' @description
#' Standardises a metabolite name to the PubChem CID, InChIKey
#' and title using PubChem synonym search.
#' Returns a single row dataframe with these as outputs.
#' Requires internet to query the PubChem database.
#'
#' @param name chr. A metabolite name to search for.
#' @param cache an RDS cache environment. Generated with `init_PubChem_cache`.
#'
#' @returns a one-row dataframe with the PubChem compound ID, InChIKey, and PubChem title
#'
#' @import curl
#' @import httr2
#'
#' @export
#'
#' @examples
#' # Set up a cache in the default location
#' cache <- init_PubChem_cache()
#'
#' # Search using PubChem synonyms
#' syn_id <- synonym_search(
#'   name = "Alanine",
#'   cache = cache)
synonym_search <- function(
    name,
    cache) {

  # Check that there is an internet connection. Stop if there isn't
  if(!curl::has_internet()) {
    stop("No internet connection available. Stopping PubChem query")
  }

  # Make sure there is a correctly connected cache
  if(is.null(cache) || (! "storr" %in% class(cache))) {
    stop(paste0('\n No cache found. Set up a cache using `init_PubChem_cache()'))
  }

  # This tryCatch checks the local cache first, then if it errors (i.e. no match)
  # it queries PubChem
  result <- tryCatch({
    # Check the local cache for the result
    # The gsub call replaces any characters forbidden from Windows filenames
    # with an underscore
    cache$get(gsub('[<>:"/\\\\|?*]', "_", name))

  }, error = function(e) {
    # tryCatch runs HTTP request and saves "NA" if it fails
    query <- tryCatch({
      # Send a programmatic request to PubChem using the PUG REST API. This query
      # reruns the compound ID (CID), InChIKey and Title.
      resp <- httr2::request(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/property/InChIKey,Title/JSON") |>
        httr2::req_url_query(name = name) |>
        httr2::req_perform()

      # Save the data to a more usable format
      query_df <- as.data.frame(httr2::resp_body_json(resp)[[1]][[1]][[1]])
      query_df[["Name"]] <- name
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
      query_df[["Name"]] <- name
      query_df$PubChem_Message <- e$message

      query_df
    })
    # Save the result to the local cache
    # The gsub call replaces any characters forbidden from Windows filenames
    # with an underscore
    cache$set(gsub('[<>:"/\\\\|?*]', "_", name), query)

    # Return query
    query
  })
}
