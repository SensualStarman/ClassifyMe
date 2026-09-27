#' Send multiple queries to pubchem and return a dataframe
#'
#' @description
#' This is designed to be a helper function which is repeatedly called via
#' `query_DB` to query the PubChem database.
#'
#' @param queries chr vector. A set of identifiers (HMDB, LIPID-MAPS or CAS) or
#' metabolite names to query PubChem with
#' @param query_fun a PubChem querying function. Either `identifier_standardiser`
#' or `synonym_search`
#' @param id_colname chr. A name of a column where the queries are stored
#' @param label chr. The type of identifier ("CAS", "HMDB", "LMP" or "Name")
#'
#' @returns a dataframe containing the PubChem CID, InChIKey, PubChem Title and
#' the identifier used for searching. Also includes PubChem error messages
#' @export
.query_and_bind <- function(queries, query_fun, id_colname, label) {
  # Query the database for each element of `queries` using `query_fun`
  results <- lapply(queries, function(query) {
    print(paste0("Querying local cache and PubChem for ", label, "# ", query, " . . ."))
    result <- query_fun(query)
    result[[4]] <- query
    print("Done!")
    # Put the results into a one-row dataframe for binding later
    data.frame(
      CID = result[[1]],
      InChIKey = result[[2]],
      Title = result[[3]],
      ID = result[[4]],
      PubChem_Message = result[[5]],
      stringsAsFactors = FALSE
    )
  })
  # Bind and return the results
  out <- do.call(rbind, results)
  names(out)[4] <- id_colname
  rownames(out) <- NULL
  out
}
