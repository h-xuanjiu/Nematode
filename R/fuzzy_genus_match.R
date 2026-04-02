#' @importFrom utils adist data
NULL
#' Fuzzy Matching of Nematode Genus Names
#'
#' @description
#' This function performs fuzzy matching of nematode genus names against a reference database
#' using Levenshtein distance (edit distance) with case insensitivity.
#'
#' @param Query.genus A \code{character vector} of genus names to be matched against the reference
#' @param max_dist Maximum allowed Levenshtein distance for matches (default = 2)
#' @param ... Additional parameters (currently unused)
#'
#' @return A data frame containing:
#' \itemize{
#'   \item Query.genus - Original query genus name
#'   \item CorrectName - Matched genus name from reference
#'   \item Distance - Edit distance between query and match
#'   \item Additional columns - All columns from nematode.info for matched records
#' }
#'
#' @examples
#' fuzzy_genus_match(c("Harterta", "Meloidogyne"))
#'
#' @export
fuzzy_genus_match <- function(Query.genus, max_dist = 2, ...) {
  nematode.info <- Nematode::nematode.info
  all_cols <- c("Query.genus", "CorrectName", "Distance", colnames(nematode.info))
  results <- list()

  for (genus in Query.genus) {

    distances <- stringdist::stringdist(
      tolower(genus),
      tolower(nematode.info$Genus),
      method = "lv"
    )
    matches <- nematode.info$Genus[distances <= max_dist]

    if (length(matches) == 0) {
      match_df <- as.data.frame(matrix(ncol = length(all_cols), nrow = 1))
      colnames(match_df) <- all_cols
      match_df$Query.genus <- genus
    } else {
      match_df <- data.frame(
        Query.genus = genus,
        CorrectName = matches,
        Distance = distances[distances <= max_dist],
        stringsAsFactors = FALSE
      )

      match_df <- match_df[order(match_df$Distance), ]

      ref_idx <- match(
        tolower(match_df$CorrectName),
        tolower(nematode.info$Genus)
      )

      ref_data <- nematode.info[ref_idx, , drop = FALSE]
      match_df <- cbind(match_df, ref_data)
    }
    results[[genus]] <- match_df
  }

  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  return(result_df)
}
