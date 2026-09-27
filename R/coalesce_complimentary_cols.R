#' Coalesce columns with complementary missingness
#'
#' @description
#' Coalesces dataframe columns where there are NA values in unique locations
#' across columns. Designed to be used after left join, but will work on any
#' set of columns that have one entry per row and are of the same type.
#'
#' @param data a dataframe with columns to be coalesced
#' @param cols a chr vector containing the names of columns to be coalesced. Must be the same type (e.g. `chr`).
#' @param output_col_name a name for the new column containing colasced data. Default is "merged_col"
#'
#' @returns The coalesced column appended to `data`, minus the columns that were coalesced
#' @import dplyr
#' @export
#'
#' @examples
#' NA_3 <- rep(NA_character_, 3)
#'
#' df_1 <- data.frame(
#'   a = c("A", "B", "C", NA_3, NA_3),
#'   b = c(NA_3, "D", "E", "F", NA_3),
#'   c = c(NA_3, NA_3, "G", "H", "I")
#' )
#'
#' # Appends a new column containing coalesced data
#' .coalesce_complimentary_cols(df_1, cols = names(df_1))
#'
#' # If there is mroe than one entry in any row, return NA for that row
#' df_2 <- df_1
#' df_2[1,3] <- "X"
#'
#' .coalesce_complimentary_cols(df_2, cols = names(df_2))
.coalesce_complimentary_cols <-
  function(
    data,
    cols,
    output_col_name = "merged_col") {
    # Subset the data with only the cols required
    subset <- data[cols]

    # Check that all vectors are the same type
    types <- sapply(subset, class)
    if(! length(unique(types)) == 1) {
      stop("Columns to be coalesced are not the same data type\n",
           paste(names(types), collapse = "\t"),
           "\n",
           paste(types), collapse = "\t")
    }

    # Calculate the number of values that aren't NA in each row. This should == 1
    # so that rows can be merged. Where they are > 1, output is NA.
    n_not_NA_logic <- rowSums(!is.na(subset)) > 1

    rownames_not_NA <- rownames(data)[n_not_NA_logic]

    # Print a warning if any rows have more than one entry
    if(sum(n_not_NA_logic > 0)) {
      warning("One of more columns has more than one entry!\n",
              "Returning NA for these row indices:\n\n",
              paste(rownames_not_NA, collapse = ", "))
    }


    # Coalesce the columns
    merged_col <- do.call(dplyr::coalesce, subset)
    merged_col[n_not_NA_logic] <- NA

    # Save and return the result
    result <- data[!names(data) %in% cols]
    result[[output_col_name]] <- merged_col

    result
}
