#' @importFrom utils getFromNamespace combn
#' @importFrom vegan metaMDS stressplot adonis2 anosim simper decostand
NULL
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("vegan", quietly = TRUE)) {
    if (exists("plot.metaMDS", envir = asNamespace("vegan"), mode = "function")) {
      plot_method <- get("plot.metaMDS", envir = asNamespace("vegan"))
      registerS3method("plot", "metaMDS", plot_method, envir = asNamespace(pkgname))
    }
    if (exists("summary.simper", envir = asNamespace("vegan"), mode = "function")) {
      summary_method <- get("summary.simper", envir = asNamespace("vegan"))
      registerS3method("summary", "simper", summary_method, envir = asNamespace(pkgname))
    }
    if (exists("summary.eigenvals", envir = asNamespace("vegan"), mode = "function")) {
      summary_method <- get("summary.eigenvals", envir = asNamespace("vegan"))
      registerS3method("summary", "eigenvals", summary_method, envir = asNamespace(pkgname))
    }
  }
}

.onAttach <- function(libname, pkgname) {
  if (interactive() && !getOption("nematode.quiet", FALSE)) {
    if (!"vegan" %in% loadedNamespaces()) {
      packageStartupMessage(
        "Note: For full functionality, please load vegan with:\n",
        "  library(vegan)"
      )
    }
  }
}


.anosim <- function(x, grouping, permutations = 999, distance = "bray", strata = NULL, parallel = getOption("mc.cores"), ...) {
  vegan::anosim(
    x = x,
    grouping = grouping,
    permutations = permutations,
    distance = distance,
    strata = strata,
    parallel = parallel
  )
}

#' Non-Metric Multidimensional Scaling (NMDS) Analysis
#'
#' @description
#' This function performs NMDS analysis on a dataset using the specified distance metric,
#' and optionally runs PERMANOVA (adonis2) and ANOSIM tests for group differences.
#' It supports both data.frame and matrix inputs.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names.
#' @param group \code{data.frame}. A data frame with sample names as row names and a single column containing group information for each sample.
#' @param distance Distance metric to use (default: "bray"). See \code{\link[vegan]{metaMDS}} for all available options.
#' @param decostand.method Standardization methods for community ecology data (default: "hellinger"). Set to `NULL` for no transformation. See \code{\link[vegan]{decostand}} for all available options.
#' @param k Number of dimensions for NMDS (default: 2).
#' @param autotransform Logical; whether to automatically transform the data (default: TRUE).
#'                      See \code{\link[vegan]{metaMDS}} for details.
#' @param adonis2 Logical; whether to perform PERMANOVA test using \code{\link[vegan]{adonis2}} (default: TRUE).
#' @param anosim Logical; whether to perform ANOSIM test using \code{\link[vegan]{anosim}} (default: TRUE).
#' @param simper Logical; whether to perform SIMPER test using \code{\link[vegan]{simper}} (default: TRUE).
#' @param ... Additional arguments passed to \code{\link[vegan]{metaMDS}}, \code{\link[vegan]{decostand}}, \code{\link[vegan]{adonis2}}, \code{\link[vegan]{anosim}}, or \code{\link[vegan]{simper}}.
#'
#' @returns An object of class "NMDS" containing:
#' \itemize{
#'   \item data - List containing the input data and group information
#'   \item call - The function call
#'   \item NMDS - NMDS results from \code{\link[vegan]{metaMDS}}
#'   \item adonis2 - PERMANOVA results (if adonis2 = TRUE)
#'   \item anosim - ANOSIM results (if anosim = TRUE)
#'   \item SIMPER - SIMPER results (if simper = TRUE)
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[vegan]{metaMDS}} for details on NMDS implementation and distance measures
#'   \item \code{\link[vegan]{decostand}} for details on standardization methods
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
#' nmds <- runNMDS(data, group = group_df)
#'
#' @export
runNMDS <- function(data, group, distance = "bray", k = 2, decostand.method = "hellinger", autotransform = TRUE, adonis2 = TRUE, anosim = TRUE, simper = TRUE, ...) {
  UseMethod("runNMDS")
}

#' @rdname runNMDS
#' @method runNMDS data.frame
#' @exportS3Method Nematode::runNMDS
runNMDS.data.frame <- function(data, group, distance = "bray", k = 2, decostand.method = "hellinger", autotransform = TRUE, adonis2 = TRUE, anosim = TRUE, simper = TRUE, ...) {
  if (is.null(rownames(data))) {
    stop("Data must have row names")
  }

  if (!all(row.names(data) %in% row.names(group)) || any(is.na(group[, 1]))) {
    stop("Please provide the group information for all samples!")
  }

  rawdata <- data
  data[is.na(data)] <- 0
  if (!is.null(decostand.method)) {
    data <- vegan::decostand(data, method = decostand.method, ...)
  }
  group <- group[row.names(data), , drop = F]

  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }

  result <- structure(
    list(
      data = list(
        data = rawdata,
        decostand.method = decostand.method,
        data.std = data,
        group = group
      ),
      call = match.call()
    ),
    class = c("Ordination", "NMDS")
  )
  result$adonis2 <- NULL
  result$anosim <- NULL
  data.nmds <- vegan::metaMDS(data, distance = distance, k = k, autotransform = autotransform, ...)
  result$NMDS <- data.nmds
  plot(data.nmds, type = "t")
  vegan::stressplot(data.nmds)
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

#' @rdname runNMDS
#' @method runNMDS matrix
#' @exportS3Method Nematode::runNMDS
runNMDS.matrix <- function(data, group, distance = "bray", k = 2, decostand.method = "hellinger", autotransform = TRUE, adonis2 = TRUE, anosim = TRUE, simper = TRUE, ...) {
  if (is.null(rownames(data))) {
    stop("Data must have row names")
  }

  if (!all(row.names(data) %in% row.names(group)) || any(is.na(group[, 1]))) {
    stop("Please provide the group information for all samples!")
  }

  rawdata <- data
  data[is.na(data)] <- 0
  if (!is.null(decostand.method)) {
    data <- vegan::decostand(data, method = decostand.method, ...)
  }
  group <- group[row.names(data), , drop = F]

  if (!is.numeric(data)) {
    stop("Input matrix contains non-numeric values")
  }

  result <- structure(
    list(
      data = list(
        data = rawdata,
        decostand.method = decostand.method,
        data.std = data,
        group = group
      ),
      call = match.call()
    ),
    class = c("Ordination", "NMDS")
  )
  result$adonis2 <- NULL
  result$anosim <- NULL
  data.nmds <- vegan::metaMDS(data, distance = distance, k = k, autotransform = autotransform, ...)
  result$NMDS <- data.nmds
  plot(data.nmds, type = "t")
  vegan::stressplot(data.nmds)
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

#' @rdname runNMDS
#' @method runNMDS default
#' @exportS3Method Nematode::runNMDS
runNMDS.default <- function(data, group, distance = "bray", k = 2, decostand.method = "hellinger", autotransform = TRUE, adonis2 = TRUE, anosim = TRUE, simper = TRUE, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}

#' Summarize NMDS Results
#'
#' Provides a concise summary of Non-Metric Multidimensional Scaling (NMDS) analysis results,
#' including stress value, PERMANOVA (adonis2) and ANOSIM test statistics.
#'
#' @param object An object of class "NMDS" produced by \code{\link{runNMDS}} function.
#' @param ... Additional arguments (currently not used).
#'
#' @return A list containing:
#' \itemize{
#'   \item stress - NMDS stress value
#'   \item points - Sample coordinates in the reduced space
#'   \item adonis2 - PERMANOVA results (R2, p-value, significance)
#'   \item anosim - ANOSIM results (R statistic, p-value, significance)
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
#' # Example for summary.NMDS
#' nmds <- runNMDS(data, group = group_df)
#' summary(nmds)
#'
#' @rdname summary
#' @method summary NMDS
#' @export
summary.NMDS <- function(object, ...) {
  result <- list(
    call = object$call,
    stress = object$NMDS$stress,
    points = object$NMDS$points,
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

#' Similarity Percentages Analysis
#'
#' @description
#' Discriminating species between two groups using Bray-Curtis dissimilarities
#'
#'
#' @param object An object of class "Ordination".
#' @param ... Additional arguments passed to \code{\link[vegan]{simper}}.
#'
#' @returns The object of class "Ordination" containing (See \code{\link{runNMDS}} for details):
#' \itemize{
#'   \item data - List containing the input data and group information
#'   \item call - The function call
#'   \item NMDS - NMDS results from \code{\link[vegan]{metaMDS}}
#'   \item SIMPER - SIMPER results
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
#' nmds <- runNMDS(data, group = group_df, simper = FALSE)
#'
#' # Example
#' nmds_simper <- runSimper(nmds)
#' print(nmds_simper$SIMPER)
#'
#' @export
runSimper <- function(object, ...) {
  UseMethod("runSimper")
}

#' @rdname runSimper
#' @method runSimper Ordination
#' @exportS3Method Nematode::runSimper
runSimper.Ordination <- function(object, ...) {
  sim_results <- vegan::simper(object$data$data.std, group = object$data$group[,1], ...)
  sim_summary <- summary(sim_results)
  pair_name <- utils::combn(unique(object$data$group[,1]), 2, FUN = function(x) paste(x, collapse = "_"))
  sim_summary_list <- lapply(sim_summary, function(df) {
    df <- data.frame(
      Genus = rownames(df),
      Contribution.Std = df$average / sum(df$average),
      df,
      Signif = ifelse(
        df$p < 0.001, "***", ifelse(
          df$p < 0.01, "**", ifelse(
            df$p < 0.05, "*", ""
          )
        )
      )
    )
    return(df)
  })
  object$SIMPER <- sim_summary_list
  return(object)
}

#' @rdname runSimper
#' @method runSimper default
#' @exportS3Method Nematode::runSimper
runSimper.default <- function(object, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}
