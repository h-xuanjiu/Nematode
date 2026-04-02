#' @importFrom utils data
NULL
#' Check validity of nematode genus names against reference database
#'
#' @description
#' This generic function validates nematode genus names by checking
#' their existence in a reference database (nematode.info). It
#' supports multiple input types and provides flexible
#' output formats.
#'
#' @param Query.genus Input to check: can be \code{character vector} or \code{data.frame}
#' @param Query.col When input is \code{data.frame}, specifies column name containing genus names
#'                 (ignored for character input)
#' @param show.details Logical controlling output format:
#' \itemize{
#'   \item \code{TRUE}: returns data.frame with query, existence status, and full reference info
#'   \item \code{FALSE}: returns only invalid/missing genus names
#' }
#' @param ... Additional arguments (currently unused).
#'
#' @return Output varies by input type and show.details:
#'   * For \code{character vector} input:
#'     - show.details = TRUE: data.frame with query, existence, and reference data
#'     - show.details = FALSE: character vector of invalid genera
#'   * For \code{data.frame} input: same as character input for the specified column
#'   * For unsupported types: error message
#'
#' @examples
#' # Check character vector
#' check_nematode_genus(c("Caenorhabditis", "Wrong"))
#'
#' # Check data.frame column
#' df <- data.frame(genus = c("Meloidogyne", "XXX"))
#' check_nematode_genus(Query.genus = df, Query.col = "genus")
#'
#' @export
check_nematode_genus <- function(Query.genus, Query.col = NULL, show.details = TRUE, ...) {
  UseMethod("check_nematode_genus")
}


#' @rdname check_nematode_genus
#' @method check_nematode_genus character
#' @exportS3Method Nematode::check_nematode_genus
check_nematode_genus.character <- function(Query.genus, Query.col = NULL, show.details = TRUE, ...) {
  # Load reference data
  if (length(Query.genus) == 0) {
    return(NULL)
  }
  nematode.info <- Nematode::nematode.info

  if (show.details) {
    # Find matches in reference data (case-insensitive)
    hit_idx <- match(tolower(Query.genus), tolower(nematode.info$Genus))

    # Create results data.frame
    details <- data.frame(
      Query.genus = Query.genus,
      Exist = !is.na(hit_idx),  # Logical indicating if match was found
      stringsAsFactors = FALSE
    )

    # Append reference data for matches
    if (any(details$Exist)) {
      details <- cbind(details, nematode.info[hit_idx, ])
    } else {
      details <- cbind(details, nematode.info[NA_integer_, , drop = FALSE])
    }

    # Reorder columns to put Query.genus and Exist first
    details <- details[, c("Query.genus", "Exist", setdiff(names(details), c("Query.genus", "Exist")))]
    rownames(details) <- NULL
    return(details)
  } else {
    # Return only invalid names (not found in reference)
    return(Query.genus[!tolower(Query.genus) %in% tolower(nematode.info$Genus)])
  }
}


#' @rdname check_nematode_genus
#' @method check_nematode_genus data.frame
#' @exportS3Method Nematode::check_nematode_genus
check_nematode_genus.data.frame <- function(Query.genus, Query.col, show.details = TRUE, ...) {
  # Extract specified column and dispatch to character method
  check_nematode_genus(Query.genus = Query.genus[[Query.col]], show.details = show.details)
}


#' @rdname check_nematode_genus
#' @method check_nematode_genus default
#' @exportS3Method Nematode::check_nematode_genus
check_nematode_genus.default <- function(Query.genus, Query.col = NULL, show.details = TRUE, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(Query.genus), collapse = "/"), call. = FALSE)
}
