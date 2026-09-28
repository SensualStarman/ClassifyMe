#' Batch query the PubChem database for Title, CID and InChIKey
#'
#' @description
#' #
#'
#'
#' @param data dataframe. A metabolomics dataset from MassHunter explorer. Should
#' contain the columns:
#' \describe{
#'   \item{Name}{chr. The putative metaobolite identification as output by MassHunter explorer}
#'   \item{Formula}{chr. The chemcial formula of the metabolite}
#'   \item{CAS}{chr. The CAS identifier}
#'   \item{HMDB}{chr. The HMDB identifier}
#'   \item{LMP}{chr. The LIPID-MAPS identifier}
#' }
#' @param cache an RDS cache environment. Generated with `init_PubChem_cache`.
#'
#' @returns the original dataframe, appended with the InChIKey, PubChem title,
#' PubChem compound ID (CID), and any error messages arising from PubChem queries
#'
#' @import purrr
#'
#' @export
#'
#' @examples
# Set up a cache in the default location
#' cache <- init_PubChem_cache()
#'
#' # Query the PubChem database
#' batch_query <- query_DB(
#'    data = explorer_metadata,
#'    cache = cache
#' )
query_DB <- function(data, cache) {
  data_rows <- nrow(data)

  CAS_logic <- .get_logic(data, "CAS")

  LMP_logic <- .get_logic(data, "LMP")
  LMP_logic[CAS_logic] <- F

  HMDB_logic <- .get_logic(data, "HMDB")
  HMDB_logic[CAS_logic | LMP_logic] <- F

  NO_ID_logic <- !( CAS_logic | LMP_logic | HMDB_logic )

  # Query PubChem with CAS IDs ----
  print("Querying local cache and PubChem with CAS IDs . . .")
  CAS_results <- .query_and_bind(
    queries    = data[["CAS"]][CAS_logic],
    query_fun  = function(q) identifier_standardiser(
      identifier = q,
      type = "CAS",
      cache = cache),
    id_colname = "CAS",
    label      = "CAS"
  )

  # Query PubChem with LIPIDMAPS IDs ----
  print("Querying local cache and PubChem with LIPIDMAPS IDs . . .")
  LMP_results <- .query_and_bind(
    queries    = data[["LMP"]][LMP_logic],
    query_fun  = function(q) identifier_standardiser(
      identifier = q,
      type = "LMP",
      cache = cache),
    id_colname = "LMP",
    label      = "LMP"
  )

  # Retry failed LIPIDMAPS lookups via synonym search ----
  print("Using synonym search for LIPIDMAPS queries that failed . . .")
  failed_LMP_indices <- which(is.na(LMP_results$CID))

  if (length(failed_LMP_indices) > 0) {
    retried <- lapply(LMP_results[["LMP"]][failed_LMP_indices], function(query) {
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

  # Query PubChem with HMDB IDs ----
  print("Querying local cache and PubChem with HMDB IDs . . .")
  HMDB_results <- .query_and_bind(
    queries    = data[["HMDB"]][HMDB_logic],
    query_fun  = function(q) {
      q_padded <- add_zeroes_to_HMDB(q)
      identifier_standardiser(
        identifier = q_padded,
        type = "HMDB",
        cache = cache)
    },
    id_colname = "HMDB",
    label      = "HMDB"
  )

  # Query PubChem with leftover names ----
  print("Querying local cache and PubChem with names . . .")
  Name_results <- .query_and_bind(
    queries    = data[["Name"]][NO_ID_logic],
    query_fun  = function(q) synonym_search(
      name = q,
      cache = cache),
    id_colname = "Name",
    label      = "Name"
  )

  # Merge each set of results with the original dataframe ----
  joined <- data %>%
    left_join(CAS_results, by = "CAS") %>%
    left_join(LMP_results, by = "LMP") %>%
    left_join(HMDB_results, by = "HMDB") %>%
    left_join(Name_results, by = "Name")

  # Then merge all the extraneous columns
  to_merge <- c("CID", "InChIKey", "Title", "PubChem_Message")

  name_sets <- lapply(to_merge, function(nm) names(joined)[str_detect(names(joined), nm)]) %>%
    purrr::set_names(to_merge)

  coalesced <- joined %>%
    .coalesce_complimentary_cols(cols = name_sets$InChIKey, output_col_name = "InChIKey") %>%
    .coalesce_complimentary_cols(cols = name_sets$Title, output_col_name = "Title") %>%
    .coalesce_complimentary_cols(cols = name_sets$CID, output_col_name = "CID") %>%
    .coalesce_complimentary_cols(cols = name_sets$PubChem_Message, output_col_name = "PubChem_Message")

  coalesced
}
