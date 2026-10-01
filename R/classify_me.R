#' Batch query the ClassyFire metabolite classification database
#'
#' @description
#' Receives a vector in InChI keys and queries the ClassyFire metabolite
#' classification database. Returns a simple dataframe with results, ready for
#' merging or further analyses. This function wraps functions from the
#' `classyfireR` package.
#'
#' @param inchikeys chr vector containing only InChI keys and `NA`
#' @param cache the file path of an SQLite cache. Default `NULL` places the
#' cache in the user's home folder.
#'
#' @returns a dataframe with columns containing the InChI key and the first four
#' classification levels from the ClassyFire database
#'
#' @import DBI
#' @import janitor
#' @import magrittr
#' @import purrr
#' @import RSQLite
#' @import classyfireR
#'
#' @export
#'
#' @examples
#' # Some example InChI keys to test
#' keys <- c(
#'   "BACDZNLMIXNCOG-JTCWOHKRSA-N",
#'   "TVHAWOPAFXXIQM-PLNGDYQASA-N",
#'   "ROFVXGGUISEHAM-UHFFFAOYSA-N",
#'   "KQTROLJZWNACJT-VKAVYKQESA-N",
#'   "GMYNCKRSFMEXPG-UHFFFAOYSA-N")
#'
#' # Initialise a cache in the default location
#' classes <- classify_me(keys)
#'
#' # Future calls are faster as they access the local cache instead of the
#' # online classyfire cache
#' classes <- classify_me(keys)
#'
#' # You can setup your cache wherever you like
#' classes_alt <- classify_me(
#'   keys,
#'  cache = file.path(Sys.getenv("HOME"), "Metabolomics", "classyfire_DB/ClassyFireCache.db")
#' )

classify_me <- function(
    inchikeys,
    cache = NULL) {
  # Work out where the cache cache should live
  if(is.null(cache)) {
    # No path given, so use the default location in the user's home folder
    user_home <- if (.Platform$OS.type == "windows") {
      Sys.getenv("USERPROFILE")
    } else {
      Sys.getenv("HOME")
    }

    DB_location <- file.path(user_home, "Documents/classyfire_DB/ClassyFireCache.db")

  } else if(!file.exists(cache)) {
    # A path was given but nothing exists there yet, so warn the user and
    # create a fresh cache file at that location
    message(paste0("'ClassyFireCache.db' not found! Automatically creating a copy.\n\n",
                   "If you have a cache saved elsewhere, press stop and put it in \n'",
                   cache,
                   "' now, or use the `cache` argument.\n",
                   "It will save you a lot of time, I promise."))
    DB_location <- cache

    # Make sure the parent folder exists before creating the cache file
    dir.create(dirname(DB_location), recursive = TRUE, showWarnings = FALSE)

    # Touch the file into existence, then close the connection immediately
    ClassyFireDB <- RSQLite::dbConnect(RSQLite::SQLite(), DB_location)

    RSQLite::dbDisconnect(ClassyFireDB)

  } else {
    # A path was given and the file already exists, so just use it
    DB_location <- cache
  }

  # Open the cache and look up each InChIKey's classification
  ClassyFireCache <- classyfireR::open_cache(dbname = DB_location)

  # Drop any missing InChIKeys before querying
  clean <- inchikeys[!is.na(inchikeys)]

  # Query the cache for every InChIKey; returns one S4 object per key
  clunky_classes <- clean %>%
    purrr::map(classyfireR::get_classification, conn = ClassyFireCache) %>%
    set_names(clean)

  # Done querying, so close the connection
  DBI::dbDisconnect(conn = ClassyFireCache)

  # Some InChIKeys won't have returned a classification; keep only those that did
  not_empty <- names(clunky_classes)[!sapply(clunky_classes, is.null)]

  # Set up an empty list to hold the tidied classification for each key
  classes <- vector(mode = "list", length = length(not_empty)) %>%
    purrr::set_names(not_empty)

  # Pull the classification dataframe out of each S4 result object
  for (index in not_empty) {
    classes[[index]] <- classification(clunky_classes[[index]])
  }

  # Classifications can have different numbers of levels, and that changes
  # how they need to be subsetted below, so split them into two groups
  nlevels <- sapply(classes, function(x) if (is.null(x)) 0 else nrow(x))
  longs <- nlevels >= 4
  shorts <- nlevels < 4

  # Reshape each classification from long (one row per level) to wide
  # (one row per compound), keeping only the top 4 levels
  classes[longs] <- lapply(
    names(classes)[longs],
    .reshape_classification,
    clunky_classes = clunky_classes,
    n = 4) %>%
    set_names(names(classes)[longs])

  # Same reshaping, but keep every level since there are fewer than 4
  classes[shorts] <- lapply(
    names(classes)[shorts],
    .reshape_classification,
    clunky_classes = clunky_classes
  ) %>%
    set_names(names(classes)[shorts])

  # Stack every compound's row into one tidy dataframe and return it
  result <- bind_rows(classes)
  result
}

# Helper: reshape one classification from long to wide, optionally
# keeping only the first `n` levels
.reshape_classification <- function(index, clunky_classes, n = NULL) {
  df <- classyfireR::classification(clunky_classes[[index]])
  if (!is.null(n)) df <- df[1:n, ]

  df %>%
    select(1:2) %>%
    t() %>%
    as.data.frame() %>%
    janitor::row_to_names(1) %>%
    mutate(InChIKey = index, .before = 1) %>%
    set_rownames(NULL)
}
