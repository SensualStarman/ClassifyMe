#' Batch query the PubChem database for Title, CID and InChIKey
#'
#' @description
#' `query_pubchem` is designed to work with MassHunter Explorer outputs, but can
#' receive any dataframe containing columns with metabolite names, CAS identifiers,
#' HMDB identifiers, or LIPID-MAPS identifiers. The function systematically
#' queries the PubChem database using the PUG-REST API and returns harmonised
#' names (the PubChem Title), PubChem compound IDs (CIDs) and an InChIKey.
#'
#' @param data dataframe. A metabolomics dataset from MassHunter explorer. To
#' use default settings, your dataframe should contain the columns:
#' \describe{
#'   \item{Name}{chr. The putative metaobolite identification as output by MassHunter explorer}
#'   \item{CAS}{chr. The CAS identifier}
#'   \item{HMDB}{chr. The HMDB identifier}
#'   \item{LMP}{chr. The LIPID-MAPS identifier}
#' }
#' @param cache an RDS cache environment. Generated with `init_PubChem_cache`.
#' @param name chr. The name of a column in `data` to override looking for
#' metabolite names in the default column "Name". Set to `NULL` to skip querying
#' PubChem with names.
#' @param cas chr. Change the "CAS" column. Behaves the same as `name`.
#' @param lipidmaps chr. Change the "LMP" column. Behaves the same as `name`.
#' @param hmdb chr. Change the "HMDB" column. Behaves the same as `name`.
#'
#' @returns the original dataframe, appended with the InChIKey, PubChem title,
#' PubChem compound ID (CID), and any error messages arising from PubChem queries
#'
#' @import purrr
#' @import stringr
#' @import dplyr
#'
#' @export
#'
#' @examples
#' # Set up a cache in the default location
#' cache <- init_PubChem_cache()
#'
#' # Query the PubChem database
#' batch_query <- query_pubchem(
#'    data = explorer_metadata,
#'    cache = cache
#' )
#' # Query the PubChem database with metabolite names only
#' name_query <- query_pubchem(
#'    data = explorer_metadata,
#'    cache = cache,
#'    cas = NULL, lipidmaps = NULL, hmdb = NULL
#' )
query_pubchem <- function(
    data,
    cache,
    name = "Name",
    cas = "CAS",
    lipidmaps = "LMP",
    hmdb = "HMDB") {
  # Calculate the number of rows
  data_rows <- nrow(data)

  # Initialise a blank dataframe to sue for null results
  null_result <- data.frame(
    CID = integer(),
    InChIKey = character(),
    Title = character(),
    ID = character(),
    PubChem_Message = character(),
    stringsAsFactors = FALSE
  )

  # CAS lookup ---
  if(is.null(cas)) {
    # Give `cas` a value, for the `by` argument of `left_join`
    cas <- "CAS"
    # Make CAS_logic all false - var required later
    CAS_logic <- vector(mode = "logical", length = data_rows)
    # Produce a blank results dataframe
    CAS_results <- null_result
    names(CAS_results)[4] <- cas

    message("No column specified for `cas`")

  } else {
    CAS_logic <- .get_logic(data, cas)

    # Query PubChem with CAS IDs
    print("Querying local cache and PubChem with CAS IDs . . .")
    CAS_results <- .query_and_bind(
      queries = data[[cas]][CAS_logic],
      query_fun = function(q) identifier_standardiser(
        identifier = q,
        type = cas,
        cache = cache),
      id_colname = cas,
      label = cas
    )
  }

  # LIPID-MAPS lookup ----
  if(is.null(lipidmaps)) {
    # Give `lipidmaps` a value, for the `by` argument of `left_join`
    lipidmaps <- "LMP"
    # Make LMP_logic all false
    LMP_logic <- vector(mode = "logical", length = data_rows)
    # Produce a blank results dataframe
    LMP_results <- null_result
    names(LMP_results)[4] <- lipidmaps

    message("No column specified for `lipidmaps`")

  } else {
    LMP_logic <- .get_logic(data, lipidmaps)
    LMP_logic[CAS_logic] <- F

    # Query PubChem with LIPIDMAPS IDs
    print("Querying local cache and PubChem with LIPIDMAPS IDs . . .")
    LMP_results <- .query_and_bind(
      queries = data[[lipidmaps]][LMP_logic],
      query_fun = function(q) identifier_standardiser(
        identifier = q,
        type = lipidmaps,
        cache = cache),
      id_colname = lipidmaps,
      label = lipidmaps
    )
    # Retry failed LIPIDMAPS lookups via synonym search
    print("Using synonym search for LIPIDMAPS queries that failed . . .")
    failed_LMP_indices <- which(is.na(LMP_results$CID))

    if (length(failed_LMP_indices) > 0) {
      retried <- lapply(LMP_results[[lipidmaps]][failed_LMP_indices], function(query) {
        print(paste0("Querying local cache and PubChem for LMP# ", query, " . . ."))
        result <- synonym_search(
          query,
          cache = cache)
        result[[4]] <- query

        data.frame(
          CID = result[[1]],
          InChIKey = result[[2]],
          Title = result[[3]],
          LMP = result[[4]],
          PubChem_Message = result[[5]],
          stringsAsFactors = FALSE
        )
      })
      LMP_results[failed_LMP_indices, ] <- do.call(rbind, retried)
    }
  }

  # HMDB lookup ----
  if(is.null(hmdb)) {
    # Give `hmdb` a value, for the `by` argument of `left_join`
    hmdb <- "HMDB"
    # Make HMDB_logic all false - var required later
    HMDB_logic <- vector(mode = "logical", length = data_rows)
    # Produce a blank results dataframe
    HMDB_results <- null_result
    names(HMDB_results)[4] <- hmdb

    message("No column specified for `hmdb`")

  } else {
    HMDB_logic <- .get_logic(data, hmdb)
    HMDB_logic[CAS_logic | LMP_logic] <- F

    # Query PubChem with HMDB IDs
    print("Querying local cache and PubChem with HMDB IDs . . .")
    HMDB_results <- .query_and_bind(
      queries = data[[hmdb]][HMDB_logic],
      query_fun  = function(q) {
        q_padded <- add_zeroes_to_HMDB(q)
        identifier_standardiser(
          identifier = q_padded,
          type = hmdb,
          cache = cache)
      },
      id_colname = hmdb,
      label = hmdb
    )
  }

  # Lookup remaining names ----
  if(is.null(name)) {
    # Give `name` a value, for the `by` argument of `left_join`
    name <- "Name"
    # Make NO_ID_logic all false - var required later
    NO_ID_logic <- vector(mode = "logical", length = data_rows)
    # Produce a blank results dataframe
    Name_results <- null_result
    names(Name_results)[4] <- name

    message("No column specified for `name`")
  } else {
    NO_ID_logic <- !( CAS_logic | LMP_logic | HMDB_logic )

    # Query PubChem with leftover names
    print("Querying local cache and PubChem with names . . .")
    Name_results <- .query_and_bind(
      queries = data[[name]][NO_ID_logic],
      query_fun = function(q) synonym_search(
        name = q,
        cache = cache),
      id_colname = name,
      label = name
    )
  }

  # Merge each set of results with the original dataframe ----
  joined <- data %>%
    dplyr::left_join(CAS_results, by = cas) %>%
    dplyr::left_join(LMP_results, by = lipidmaps) %>%
    dplyr::left_join(HMDB_results, by = hmdb) %>%
    dplyr::left_join(Name_results, by = name)

  # Then merge all the extraneous columns
  to_merge <- c("CID", "InChIKey", "Title", "PubChem_Message")

  name_sets <-
    lapply(
      to_merge,
      function(nm) names(joined)[stringr::str_detect(names(joined), nm)]) %>%
    purrr::set_names(to_merge)

  coalesced <- joined %>%
    .coalesce_complimentary_cols(cols = name_sets$InChIKey, output_col_name = "InChIKey") %>%
    .coalesce_complimentary_cols(cols = name_sets$Title, output_col_name = "Title") %>%
    .coalesce_complimentary_cols(cols = name_sets$CID, output_col_name = "CID") %>%
    .coalesce_complimentary_cols(cols = name_sets$PubChem_Message, output_col_name = "PubChem_Message")

  coalesced
}
