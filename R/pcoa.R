.cmdscale <- function(d, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE, list. = eig || add || x.ret, ...) {
  result <- stats::cmdscale(
    d, k = k, eig = eig, add = add, x.ret = x.ret,
    list. = eig || add || x.ret
  )
  if (eig || add || x.ret) {
    return(result)
  } else {
    return(list(points = result))
  }
}


#' Principal Coordinates Analysis (PCoA) Analysis
#'
#' @description
#' This function performs PCoA analysis on a dataset using the specified distance metric,
#' and optionally runs PERMANOVA (adonis2) and ANOSIM tests for group differences.
#' It supports both data.frame and matrix inputs.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names.
#' @param group \code{data.frame}. A data frame with sample names as row names and a single column containing group information for each sample.
#' @param distance Distance metric to use (default: "bray"). See \code{\link[vegan]{vegdist}} for all available options.
#' @param k Number of dimensions for PCoA (default: 2).
#' @param adonis2 Logical; whether to perform PERMANOVA test using \code{\link[vegan]{adonis2}} (default: TRUE).
#' @param anosim Logical; whether to perform ANOSIM test using \code{\link[vegan]{anosim}} (default: TRUE).
#' @param simper Logical; whether to perform SIMPER test using \code{\link[vegan]{simper}} (default: TRUE).
#' @param ... Additional arguments passed to \code{\link[stats]{cmdscale}}, \code{\link[vegan]{decostand}}, \code{\link[vegan]{adonis2}}, \code{\link[vegan]{anosim}}, or \code{\link[vegan]{simper}}.
#'
#' @returns An object of class "PCoA" containing:
#' \itemize{
#'   \item data - List containing the input data and group information
#'   \item call - The function call
#'   \item Points - Sample coordinates in the reduced space.
#'   \item Eigenvalues - Variance explained by each principal coordinate axis.
#'   \item adonis2 - PERMANOVA results (if adonis2 = TRUE)
#'   \item anosim - ANOSIM results (if anosim = TRUE)
#'   \item SIMPER - SIMPER results (if simper = TRUE)
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[stats]{cmdscale}} for details on cmdscale implementation
#'   \item \code{\link[vegan]{vegdist}} for available distance metrics
#'   \item \code{\link[vegan]{adonis2}} for PERMANOVA
#'   \item \code{\link[vegan]{anosim}} for ANOSIM
#'   \item \code{\link[vegan]{simper}} for SIMPER
#' }
#'
#' @examples
#' # Example with default Bray-Curtis distance
#' data <- data.frame(
#'   Cephalobus = c(10, 20, 30, 1, 6, 5),
#'   Eucephalobus = c(5, 10, 12, 30, 1, 6),
#'   Acrobeloides = c(1, 2, 3, 12, 30, 1),
#'   Caenorhabditis = c(5, 8, 15, 2, 3, 12),
#'   Aphelenchus = c(5, 13, 11, 15, 2, 3),
#'   Leptonchus = c(3, 10, 15, 0, 15, 11),
#'   Pratylenchus = c(9, 2, 15, 15, 0, 15),
#'   Tylenchus = c(5, 0, 15, 11, 15, 2),
#'   Mesodorylaimus = c(7, 10, 18, 3, 12, 30),
#'   Discolaimus = c(1, 10, 25, 10, 18, 3),
#'   row.names = c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6")
#' )
#' group_df <- data.frame(
#'   group = c("A", "A", "B", "B", "C", "C"),
#'   row.names = c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6")
#' )
#'
#' pcoa <- runPCoA(data, group = group_df)
#'
#' @export
runPCoA <- function(data, group, k = 2, distance = "bray", adonis2 = TRUE, anosim = TRUE, simper = TRUE, ...) {
  UseMethod("runPCoA")
}

#' @rdname runPCoA
#' @method runPCoA data.frame
#' @exportS3Method Nematode::runPCoA
runPCoA.data.frame <- function(data, group, k = 2, distance = "bray", adonis2 = TRUE, anosim = TRUE, simper = TRUE, ...) {
  if (is.null(rownames(data))) {
    stop("Data must have row names")
  }

  if (!all(row.names(data) %in% row.names(group)) || any(is.na(group[, 1]))) {
    stop("Please provide the group information for all samples!")
  }

  rawdata <- data
  data[is.na(data)] <- 0
  group <- group[row.names(data), , drop = F]

  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  result <- structure(
    list(
      data = list(
        data = rawdata,
        data.std = data,
        group = group
      ),
      call = match.call()
    ),
    class = c("Ordination", "PCoA")
  )
  result$adonis2 <- NULL
  result$anosim <- NULL
  dist <- vegan::vegdist(data, method = distance, ...)
  pcoa <- .cmdscale(dist, k = k, eig = TRUE, ...)
  eig <- summary(vegan::eigenvals(pcoa))
  axis <- paste0("PCo", 1:ncol(eig))
  eig <- data.frame(Axis = axis, t(eig))
  Points <- pcoa$points
  colnames(Points) <- paste0("PCo", 1:ncol(Points))
  result$Points <- Points
  result$Eigenvalues <- eig
  if (adonis2) {
    result$adonis2 <- vegan::adonis2(data ~ group[, 1], method = distance, ...)
  }
  if (anosim) {
    result$anosim <- .anosim(data, group[, 1], distance = distance, ...)
  }
  if (simper) {
    result <- runSimper.Ordination(result, ...)
  }
  return(result)
}

#' @rdname runPCoA
#' @method runPCoA matrix
#' @exportS3Method Nematode::runPCoA
runPCoA.matrix <- function(data, group, k = 2, distance = "bray", adonis2 = TRUE, anosim = TRUE, simper = TRUE, ...) {
  if (is.null(rownames(data))) {
    stop("Data must have row names")
  }

  if (!all(row.names(data) %in% row.names(group)) || any(is.na(group[, 1]))) {
    stop("Please provide the group information for all samples!")
  }

  rawdata <- data
  data[is.na(data)] <- 0
  group <- group[row.names(data), , drop = F]

  if (!is.numeric(data)) {
    stop("Input matrix contains non-numeric values")
  }
  result <- structure(
    list(
      data = list(
        data = rawdata,
        data.std = data,
        group = group
      ),
      call = match.call()
    ),
    class = c("Ordination", "PCoA")
  )
  result$adonis2 <- NULL
  result$anosim <- NULL
  dist <- vegan::vegdist(data, method = distance, ...)
  pcoa <- .cmdscale(dist, k = k, eig = TRUE, ...)
  eig <- summary(vegan::eigenvals(pcoa))
  axis <- paste0("PCo", 1:ncol(eig))
  eig <- data.frame(Axis = axis, t(eig))
  Points <- pcoa$points
  colnames(Points) <- paste0("PCo", 1:ncol(Points))
  result$Points <- Points
  result$Eigenvalues <- eig
  if (adonis2) {
    result$adonis2 <- vegan::adonis2(data ~ group[, 1], method = distance, ...)
  }
  if (anosim) {
    result$anosim <- .anosim(data, group[, 1], distance = distance, ...)
  }
  if (simper) {
    result <- runSimper.Ordination(result, ...)
  }
  return(result)
}

#' @rdname runPCoA
#' @method runPCoA default
#' @exportS3Method Nematode::runPCoA
runPCoA.default <- function(data, group, k = 2, distance = "bray", adonis2 = TRUE, anosim = TRUE, simper = TRUE, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}


#' Summarize PCoA Results
#'
#' Provides a concise summary of Principal Coordinates Analysis (PCoA) analysis results,
#' including PERMANOVA (adonis2) and ANOSIM test statistics.
#'
#' @param object An object of class "PCoA" produced by \code{\link{runPCoA}} function.
#' @param ... Additional arguments (currently not used).
#'
#' @return A list containing:
#' \itemize{
#'   \item points - Sample coordinates in the reduced space.
#'   \item eig - Variance explained by each principal coordinate axis.
#'   \item adonis2 - PERMANOVA results (R2, p-value, significance)
#'   \item anosim - ANOSIM results (R statistic, p-value, significance)
#' }
#'
#' @examples
#' # Example for summary.PCoA
#' pcoa <- runPCoA(data, group = group_df)
#' summary(pcoa)
#'
#' @rdname summary
#' @method summary PCoA
#' @export
summary.PCoA <- function(object, ...) {
  result <- list(
    call = object$call,
    points = object$Points,
    eig = object$Eigenvalues,
    adonis2 = list(
      adonis2.R2 = object$adonis2$R2[1],
      adonis2.p = object$adonis2$`Pr(>F)`[1],
      adonis2.Sig = ifelse(
        object$adonis2$`Pr(>F)`[1] < 0.001, "***", ifelse(
          object$adonis2$`Pr(>F)`[1] < 0.01, "**", ifelse(
            object$adonis2$`Pr(>F)`[1] < 0.05, "*", ""
          )
        )
      )
    ),
    anosim = list(
      anosim.R = object$anosim$statistic,
      anosim.p = object$anosim$signif,
      anosim.Sig = ifelse(
        object$anosim$signif < 0.001, "***", ifelse(
          object$anosim$signif < 0.01, "**", ifelse(
            object$anosim$signif < 0.05, "*", ""
          )
        )
      )
    )
  )
  return(result)
}
