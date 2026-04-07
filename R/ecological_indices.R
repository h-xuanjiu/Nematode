#' @importFrom dplyr %>% mutate across everything select all_of any_of filter pull
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom purrr reduce
NULL
# ======relative abundance=====
#' Calculate the Relative Abundance of Nematodes
#'
#' @description
#' This function calculates the relative abundance of nematodes for each sample.
#' The relative abundance is defined as the proportion of each nematode's count to the total count of all nematodes in a sample.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematodes.
#' Each element indicates the count of a specific nematode in the corresponding sample.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{data.frame} or \code{matrix} (matching the input type) containing the relative abundance of each nematode in each sample.
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Species1 = c(10, NA, 15),
#'   Species2 = c(5, 10, NA),
#'   Species3 = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' rel_abundance(df)
#'
#' # Example with a matrix
#' mat <- matrix(c(10, NA, 15, 5, 10, NA, 8, 12, 10), nrow = 3, byrow = TRUE)
#' colnames(mat) <- c("Species1", "Species2", "Species3")
#' row.names(mat) <- c("A", "B", "C")
#' rel_abundance(mat)
#'
#' @export
rel_abundance <- function(data, ...) {
  UseMethod("rel_abundance")
}

#' @rdname rel_abundance
#' @method rel_abundance data.frame
#' @exportS3Method Nematode::rel_abundance
rel_abundance.data.frame <- function(data, ...) {
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  data[is.na(data)] <- 0
  total_abundance <- rowSums(data, na.rm = TRUE)
  if (any(total_abundance == 0)) {
    warning("Some rows sum to zero.")
  }
  result <- data / total_abundance
  result[] <- lapply(result, function(x) {
    if (is.numeric(x)) replace(x, is.nan(x), 0) else x
  })
  return(result)
}

#' @rdname rel_abundance
#' @method rel_abundance matrix
#' @exportS3Method Nematode::rel_abundance
rel_abundance.matrix <- function(data, ...) {
  if (!is.numeric(data)) {
    stop("The data contains non-numeric characters!")
  }
  data[is.na(data)] <- 0
  total_abundance <- rowSums(data, na.rm = TRUE)
  if (any(total_abundance == 0)) {
    warning("Some rows sum to zero.")
  }
  result <- data / total_abundance
  result[is.nan(result)] <- 0
  return(result)
}

#' @rdname rel_abundance
#' @method rel_abundance default
#' @exportS3Method Nematode::rel_abundance
rel_abundance.default <- function(data, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"),
    call. = FALSE
  )
}

# =====Number of species=====
#' Calculate Number of Species
#'
#' @description
#' This function calculates the number of nematode species present in each sample.
#' It counts the number of non-zero and non-empty nematode species for each sample.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance
#' table where rows represent samples and columns represent nematodes.
#' Each element indicates the count of a specific nematode in the
#' corresponding sample. Row names must be sample names.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{data.frame} with two columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names)}
#'   \item{NumSpecies}{Number of non-zero nematode species in each sample}
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Species1 = c(10, NA, 15),
#'   Species2 = c(5, 10, NA),
#'   Species3 = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' num_species(df)
#'
#' # Example with a matrix
#' mat <- matrix(c(10, NA, 15, 5, 10, NA, 8, 12, 10), nrow = 3, byrow = TRUE)
#' colnames(mat) <- c("Species1", "Species2", "Species3")
#' row.names(mat) <- c("A", "B", "C")
#' num_species(mat)
#'
#' @export
num_species <- function(data, ...) {
  UseMethod("num_species")
}

#' @rdname num_species
#' @method num_species data.frame
#' @exportS3Method Nematode::num_species
num_species.data.frame <- function(data, ...) {
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  if (is.null(rownames(data))) {
    stop("Data must have row names")
  }
  count <- apply(data, 1, function(row) {
    sum(row > 0, na.rm = TRUE)
  })
  result <- data.frame(
    Sample.ID = row.names(data),
    NumSpecies = count,
    stringsAsFactors = FALSE
  )
  return(result)
}

#' @rdname num_species
#' @method num_species matrix
#' @exportS3Method Nematode::num_species
num_species.matrix <- function(data, ...) {
  if (!is.numeric(data)) {
    stop("The data contains non-numeric characters!")
  }
  if (is.null(rownames(data))) {
    stop("Data must have row names")
  }
  count <- apply(data, 1, function(row) {
    sum(row > 0, na.rm = TRUE)
  })
  result <- data.frame(
    Sample.ID = row.names(data),
    NumSpecies = count,
    stringsAsFactors = FALSE
  )
  return(result)
}

#' @rdname num_species
#' @method num_species default
#' @exportS3Method Nematode::num_species
num_species.default <- function(data, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}

# =====diet_rel_abundance=====
#' Calculate Diet Relative or Absolute Abundance
#'
#' This function calculates the relative or absolute abundance of four feeding types of nematodes in each sample.
#' The feeding types include bacterial feeders (Ba), fungus feeders (Fu), plant feeders (Pp), and omnivores/predators (Op).
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param total.abundance \code{data.frame}. A data frame with sample names as row names and a single column containing the total nematode abundance for each sample.
#' This parameter is required when \code{relative} is set to \code{FALSE}. Default is \code{NULL}.
#' @param relative \code{Logical}. If \code{TRUE} (default), the function calculates relative abundance (does not require \code{total.abundance}).
#' If \code{FALSE}, the function calculates absolute abundance (requires \code{total.abundance}).
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame with five columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names of \code{data})}
#'   \item{Ba}{Relative or absolute abundance of bacterial feeders}
#'   \item{Fu}{Relative or absolute abundance of fungus feeders}
#'   \item{Pp}{Relative or absolute abundance of plant feeders}
#'   \item{Op}{Relative or absolute abundance of omnivores/predators}
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' abundance <- data.frame(
#'   abundance = c(100, 150, 120),
#'   row.names = c("A", "B", "C")
#' )
#' diet_rel_abundance(df, abundance, relative = FALSE)
#'
#' # Example with a matrix
#' mat <- matrix(c(10, NA, 15, 5, 10, NA, 8, 12, 10), nrow = 3, byrow = TRUE)
#' colnames(mat) <- c("Cephalobus", "Caenorhabditis", "Pratylenchus")
#' row.names(mat) <- c("A", "B", "C")
#' diet_rel_abundance(mat)
#' @export
diet_rel_abundance <- function(data, total.abundance = NULL, relative = TRUE, ...) {
  UseMethod("diet_rel_abundance")
}

#' @rdname diet_rel_abundance
#' @method diet_rel_abundance data.frame
#' @exportS3Method Nematode::diet_rel_abundance
diet_rel_abundance.data.frame <- function(data, total.abundance = NULL, relative = TRUE, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  if (!relative) {
    if (is.null(total.abundance)) {
      stop("total.abundance must be provided when relative = FALSE.")
    }
    if (!all(row.names(data) %in% row.names(total.abundance))) {
      stop("Please provide the total.abundance information for all samples!")
    }
    if (!all(sapply(total.abundance, is.numeric))) {
      stop("The total.abundance contains non-numeric characters!")
    }
  }
  genus_info <- check_nematode_genus(colnames(data))
  if (!all(genus_info$Exist)) {
    missing <- genus_info$Query.genus[!genus_info$Exist]
    stop(
      length(missing), " genus names not found: ", paste(head(missing), collapse = ", "),
      ifelse(length(missing) > 6, "...", "")
    )
  }
  valid_habits <- c("Bacterial feeders", "Fungus feeders", "Plant feeders", "Omnivores", "Predators")
  invalid_habits <- setdiff(unique(genus_info$Feeding_habit), valid_habits)
  if (length(invalid_habits) > 0) {
    stop("Unsupported feeding habit(s): ", paste(invalid_habits, collapse = ", "))
  }
  rel_data <- rel_abundance(data)
  if (!relative) {
    abundance_value <- total.abundance[
      match(row.names(data), row.names(total.abundance)),
      colnames(total.abundance)[1]
    ]
    rel_data <- rel_data %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x * abundance_value))
  }

  type_map <- c(
    "Bacterial feeders" = "Ba",
    "Fungus feeders" = "Fu",
    "Plant feeders" = "Pp",
    "Omnivores" = "Op",
    "Predators" = "Op"
  )
  genus_info$type <- type_map[genus_info$Feeding_habit]
  conditions <- list(
    Ba = genus_info$type == "Ba",
    Fu = genus_info$type == "Fu",
    Pp = genus_info$type == "Pp",
    Op = genus_info$type == "Op"
  )

  cal_fun <- function(data, genus_info, condition) {
    sel_genus <- genus_info %>%
      dplyr::filter(!!condition) %>%
      dplyr::pull("Genus")
    data %>%
      dplyr::select(dplyr::all_of(sel_genus)) %>%
      rowSums()
  }
  result <- lapply(names(conditions), function(names) {
    cal_fun(rel_data, genus_info, conditions[[names]])
  })
  result_df <- as.data.frame(result)
  colnames(result_df) <- names(conditions)

  result_df <- data.frame(
    Sample.ID = row.names(result_df),
    result_df
  )
  row.names(result_df) <- NULL
  return(result_df)
}

#' @rdname diet_rel_abundance
#' @method diet_rel_abundance matrix
#' @exportS3Method Nematode::diet_rel_abundance
diet_rel_abundance.matrix <- function(data, total.abundance = NULL, relative = TRUE, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }
  if (!relative) {
    if (is.null(total.abundance)) {
      stop("total.abundance must be provided when relative = FALSE.")
    }
    if (!all(row.names(data) %in% row.names(total.abundance))) {
      stop("Please provide the total.abundance information for all samples!")
    }
    if (!all(sapply(total.abundance, is.numeric))) {
      stop("The total.abundance contains non-numeric characters!")
    }
  }
  genus_info <- check_nematode_genus(colnames(data))
  if (!all(genus_info$Exist)) {
    missing <- genus_info$Query.genus[!genus_info$Exist]
    stop(
      length(missing), " genus names not found: ", paste(head(missing), collapse = ", "),
      ifelse(length(missing) > 6, "...", "")
    )
  }
  valid_habits <- c("Bacterial feeders", "Fungus feeders", "Plant feeders", "Omnivores", "Predators")
  invalid_habits <- setdiff(unique(genus_info$Feeding_habit), valid_habits)
  if (length(invalid_habits) > 0) {
    stop("Unsupported feeding habit(s): ", paste(invalid_habits, collapse = ", "))
  }
  rel_data <- rel_abundance(data)
  if (!relative) {
    abundance_value <- total.abundance[match(row.names(data), row.names(total.abundance)), colnames(total.abundance)[1]]
    rel_data <- rel_data * abundance_value
  }
  type_map <- c(
    "Bacterial feeders" = "Ba",
    "Fungus feeders" = "Fu",
    "Plant feeders" = "Pp",
    "Omnivores" = "Op",
    "Predators" = "Op"
  )
  genus_info$type <- type_map[genus_info$Feeding_habit]
  conditions <- list(
    Ba = genus_info$type == "Ba",
    Fu = genus_info$type == "Fu",
    Pp = genus_info$type == "Pp",
    Op = genus_info$type == "Op"
  )
  cal_fun <- function(data_mat, genus_info, condition) {
    sel_genus <- genus_info[["Genus"]][condition]
    rowSums(data_mat[, sel_genus, drop = FALSE])
  }
  result <- lapply(names(conditions), function(names) {
    cal_fun(rel_data, genus_info, conditions[[names]])
  })
  result_mat <- do.call(cbind, result)
  colnames(result_mat) <- names(conditions)
  result_df <- data.frame(
    Sample.ID = row.names(result_mat),
    result_mat
  )
  row.names(result_df) <- NULL
  return(result_df)
}

#' @rdname diet_rel_abundance
#' @method diet_rel_abundance default
#' @exportS3Method Nematode::diet_rel_abundance
diet_rel_abundance.default <- function(data, total.abundance = NULL, relative = TRUE, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}


# =====cp_rel_abundance=====
#' Calculate CP Relative or Absolute Abundance
#'
#' This function calculates the relative or absolute abundance of nematodes in different CP (Colonizer-Persister) groups for each sample.
#' The CP groups range from CP1 (colonizers, r-strategists) to CP5 (persisters, K-strategists). Genera without CP classification are grouped as No_CP.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param total.abundance \code{data.frame}. A data frame with sample names as row names and a single column containing the total nematode abundance for each sample.
#' This parameter is required when \code{relative} is set to \code{FALSE}. Default is \code{NULL}.
#' @param relative \code{Logical}. If \code{TRUE} (default), the function calculates relative abundance (does not require \code{total.abundance}).
#' If \code{FALSE}, the function calculates absolute abundance (requires \code{total.abundance}).
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame with seven columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names of \code{data})}
#'   \item{CP1}{Relative or absolute abundance of CP1 group (colonizers, r-strategists)}
#'   \item{CP2}{Relative or absolute abundance of CP2 group}
#'   \item{CP3}{Relative or absolute abundance of CP3 group}
#'   \item{CP4}{Relative or absolute abundance of CP4 group}
#'   \item{CP5}{Relative or absolute abundance of CP5 group (persisters, K-strategists)}
#'   \item{No_CP}{Relative or absolute abundance of genera without CP classification (only present if such genera exist)}
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' abundance <- data.frame(
#'   abundance = c(100, 150, 120),
#'   row.names = c("A", "B", "C")
#' )
#' cp_rel_abundance(df, abundance, relative = FALSE)
#'
#' # Example with a matrix
#' mat <- matrix(c(10, NA, 15, 5, 10, NA, 8, 12, 10), nrow = 3, byrow = TRUE)
#' colnames(mat) <- c("Cephalobus", "Caenorhabditis", "Pratylenchus")
#' row.names(mat) <- c("A", "B", "C")
#' cp_rel_abundance(mat)
#' @export
cp_rel_abundance <- function(data, total.abundance = NULL, relative = TRUE, ...) {
  UseMethod("cp_rel_abundance")
}

#' @rdname cp_rel_abundance
#' @method cp_rel_abundance data.frame
#' @exportS3Method Nematode::cp_rel_abundance
cp_rel_abundance.data.frame <- function(data, total.abundance = NULL, relative = TRUE, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  if (!relative) {
    if (is.null(total.abundance)) {
      stop("total.abundance must be provided when relative = FALSE.")
    }
    if (!all(row.names(data) %in% row.names(total.abundance))) {
      stop("Please provide the total.abundance information for all samples!")
    }
    if (!all(sapply(total.abundance, is.numeric))) {
      stop("The total.abundance contains non-numeric characters!")
    }
  }
  genus_info <- check_nematode_genus(colnames(data))
  if (!all(genus_info$Exist)) {
    missing <- genus_info$Query.genus[!genus_info$Exist]
    stop(
      length(missing), " genus names not found: ", paste(head(missing), collapse = ", "),
      ifelse(length(missing) > 6, "...", "")
    )
  }
  valid_habits <- c("Bacterial feeders", "Fungus feeders", "Plant feeders", "Omnivores", "Predators")
  invalid_habits <- setdiff(unique(genus_info$Feeding_habit), valid_habits)
  if (length(invalid_habits) > 0) {
    stop("Unsupported feeding habit(s): ", paste(invalid_habits, collapse = ", "))
  }
  rel_data <- rel_abundance(data)
  if (!relative) {
    abundance_value <- total.abundance[
      match(row.names(data), row.names(total.abundance)),
      colnames(total.abundance)[1]
    ]
    rel_data <- rel_data %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x * abundance_value))
  }

  conditions <- list(
    CP1 = genus_info$CP_group == 1,
    CP2 = genus_info$CP_group == 2,
    CP3 = genus_info$CP_group == 3,
    CP4 = genus_info$CP_group == 4,
    CP5 = genus_info$CP_group == 5
  )

  if (any(is.na(genus_info$CP_group))) {
    conditions$No_CP <- is.na(genus_info$CP_group)
  }

  cal_fun <- function(data, genus_info, condition) {
    sel_genus <- genus_info %>%
      dplyr::filter(!!condition) %>%
      dplyr::pull("Genus")
    data %>%
      dplyr::select(dplyr::all_of(sel_genus)) %>%
      rowSums()
  }
  result <- lapply(names(conditions), function(names) {
    cal_fun(rel_data, genus_info, conditions[[names]])
  })
  result_df <- as.data.frame(result)
  colnames(result_df) <- names(conditions)
  result_df <- data.frame(
    Sample.ID = row.names(result_df),
    result_df
  )
  row.names(result_df) <- NULL
  return(result_df)
}

#' @rdname cp_rel_abundance
#' @method cp_rel_abundance matrix
#' @exportS3Method Nematode::cp_rel_abundance
cp_rel_abundance.matrix <- function(data, total.abundance = NULL, relative = TRUE, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }
  if (!relative) {
    if (is.null(total.abundance)) {
      stop("total.abundance must be provided when relative = FALSE.")
    }
    if (!all(row.names(data) %in% row.names(total.abundance))) {
      stop("Please provide the total.abundance information for all samples!")
    }
    if (!all(sapply(total.abundance, is.numeric))) {
      stop("The total.abundance contains non-numeric characters!")
    }
  }
  genus_info <- check_nematode_genus(colnames(data))
  if (!all(genus_info$Exist)) {
    missing <- genus_info$Query.genus[!genus_info$Exist]
    stop(
      length(missing), " genus names not found: ", paste(head(missing), collapse = ", "),
      ifelse(length(missing) > 6, "...", "")
    )
  }
  valid_habits <- c("Bacterial feeders", "Fungus feeders", "Plant feeders", "Omnivores", "Predators")
  invalid_habits <- setdiff(unique(genus_info$Feeding_habit), valid_habits)
  if (length(invalid_habits) > 0) {
    stop("Unsupported feeding habit(s): ", paste(invalid_habits, collapse = ", "))
  }
  rel_data <- rel_abundance(data)
  if (!relative) {
    abundance_value <- total.abundance[match(row.names(data), row.names(total.abundance)), colnames(total.abundance)[1]]
    rel_data <- rel_data * abundance_value
  }

  conditions <- list(
    CP1 = genus_info$CP_group == 1,
    CP2 = genus_info$CP_group == 2,
    CP3 = genus_info$CP_group == 3,
    CP4 = genus_info$CP_group == 4,
    CP5 = genus_info$CP_group == 5
  )

  if (any(is.na(genus_info$CP_group))) {
    conditions$No_CP <- is.na(genus_info$CP_group)
  }
  cal_fun <- function(data_mat, genus_info, condition) {
    sel_genus <- genus_info[["Genus"]][condition]
    rowSums(data_mat[, sel_genus, drop = FALSE])
  }
  result <- lapply(names(conditions), function(names) {
    cal_fun(rel_data, genus_info, conditions[[names]])
  })
  result_mat <- do.call(cbind, result)
  colnames(result_mat) <- names(conditions)
  result_df <- data.frame(
    Sample.ID = row.names(result_mat),
    result_mat
  )
  row.names(result_df) <- NULL
  return(result_df)
}

#' @rdname cp_rel_abundance
#' @method cp_rel_abundance default
#' @exportS3Method Nematode::cp_rel_abundance
cp_rel_abundance.default <- function(data, total.abundance = NULL, relative = TRUE, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}



# =====Trophic diversity=====
#' Calculate Trophic Diversity (TD) Index
#'
#' @description
#' This function calculates the Trophic Diversity (TD) Index for ecological communities.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame with two columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names of \code{data})}
#'   \item{TD}{Trophic Diversity index for each sample}
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' cal.TD(data = df)
#' @export
cal.TD <- function(data, ...) {
  UseMethod("cal.TD")
}

#' @rdname cal.TD
#' @method cal.TD data.frame
#' @exportS3Method Nematode::cal.TD
cal.TD.data.frame <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  diet_rel <- diet_rel_abundance(data = data, total.abundance = NULL, relative = T)
  row.names(diet_rel) <- diet_rel$Sample.ID
  diet_rel <- diet_rel[, -1, drop = FALSE]^2
  TD <- data.frame(
    Sample.ID = rownames(diet_rel),
    TD = 1 / rowSums(diet_rel)
  )
  return(TD)
}

#' @rdname cal.TD
#' @method cal.TD matrix
#' @exportS3Method Nematode::cal.TD
cal.TD.matrix <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }
  diet_rel <- diet_rel_abundance(data = data, total.abundance = NULL, relative = T)
  row.names(diet_rel) <- diet_rel$Sample.ID
  diet_rel <- diet_rel[, -1, drop = FALSE]^2
  TD <- data.frame(
    Sample.ID = rownames(diet_rel),
    TD = 1 / rowSums(diet_rel)
  )
  return(TD)
}

#' @rdname cal.TD
#' @method cal.TD default
#' @exportS3Method Nematode::cal.TD
cal.TD.default <- function(data, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}

# =====Simpson index=====
#' Calculate Simpson Index
#'
#' @description
#' This function calculates the Simpson Index for ecological communities.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame with two columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names of \code{data})}
#'   \item{Simpson}{Simpson's Index for each sample}
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' cal.Simpson(data = df)
#' @export
cal.Simpson <- function(data, ...) {
  UseMethod("cal.Simpson")
}

#' @rdname cal.Simpson
#' @method cal.Simpson data.frame
#' @exportS3Method Nematode::cal.Simpson
cal.Simpson.data.frame <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  rel_data <- rel_abundance(data = data)^2
  Simpson <- data.frame(
    Sample.ID = row.names(rel_data),
    Simpson = 1 - rowSums(rel_data)
  )
  return(Simpson)
}

#' @rdname cal.Simpson
#' @method cal.Simpson matrix
#' @exportS3Method Nematode::cal.Simpson
cal.Simpson.matrix <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }
  rel_data <- rel_abundance(data = data)^2
  Simpson <- data.frame(
    Sample.ID = row.names(rel_data),
    Simpson = 1 - rowSums(rel_data)
  )
  return(Simpson)
}

#' @rdname cal.Simpson
#' @method cal.Simpson default
#' @exportS3Method Nematode::cal.Simpson
cal.Simpson.default <- function(data, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}


# =====Shannon-Wiener index======
#' Calculate Shannon-Wiener Index (H)
#'
#' @description
#' This function calculates the Shannon-Wiener Index (H) for ecological communities.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame with two columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names of \code{data})}
#'   \item{H}{Shannon-Wiener Index for each sample}
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' cal.H(data = df)
#' @export
cal.H <- function(data, ...) {
  UseMethod("cal.H")
}

#' @rdname cal.H
#' @method cal.H data.frame
#' @exportS3Method Nematode::cal.H
cal.H.data.frame <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  rel_data <- rel_abundance(data = data)
  log_rel <- log(rel_data) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ ifelse(. == -Inf, 0, .)))

  metaH <- rel_data * log_rel
  H <- data.frame(
    Sample.ID = row.names(metaH),
    H = -rowSums(metaH)
  )
  return(H)
}

#' @rdname cal.H
#' @method cal.H matrix
#' @exportS3Method Nematode::cal.H
cal.H.matrix <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }
  rel_data <- rel_abundance(data = data)
  log_rel <- log(rel_data)
  log_rel[log_rel == -Inf] <- 0

  metaH <- rel_data * log_rel
  H <- data.frame(
    Sample.ID = row.names(metaH),
    H = -rowSums(metaH)
  )
  return(H)
}

#' @rdname cal.H
#' @method cal.H default
#' @exportS3Method Nematode::cal.H
cal.H.default <- function(data, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}

# =====Pielou's Evenness Index=====
#' Calculate Pielou's Evenness Index (J)
#'
#' @description
#' This function calculates the Pielou's Evenness Index (J) for ecological communities.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame with two columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names of \code{data})}
#'   \item{J}{Pielou's Evenness Index for each sample}
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' cal.J(data = df)
#' @export
cal.J <- function(data, ...) {
  UseMethod("cal.J")
}

#' @rdname cal.J
#' @method cal.J data.frame
#' @exportS3Method Nematode::cal.J
cal.J.data.frame <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  H <- cal.H(data)
  N <- num_species(data)
  H_N <- merge(H, N, by = "Sample.ID")
  J <- data.frame(
    Sample.ID = H_N$Sample.ID,
    J = H_N$H / log(H_N$NumSpecies)
  )
  return(J)
}

#' @rdname cal.J
#' @method cal.J matrix
#' @exportS3Method Nematode::cal.J
cal.J.matrix <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }
  H <- cal.H(data)
  N <- num_species(data)
  H_N <- merge(H, N, by = "Sample.ID")
  J <- data.frame(
    Sample.ID = H_N$Sample.ID,
    J = H_N$H / log(H_N$NumSpecies)
  )
  return(J)
}

#' @rdname cal.J
#' @method cal.J default
#' @exportS3Method Nematode::cal.J
cal.J.default <- function(data, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}

# =====Species Richness Index=====
#' Calculate Species Richness Index (SRI)
#'
#' @description
#' This function calculates the Species Richness Index (SRI) for ecological communities.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param method The method used to calculate the Species Richness Index. Default is \code{"Margalef"}.
#' Supported methods are \code{"Margalef"} and \code{"Menhinick"}. Only one method can be specified.
#'   \itemize{
#'     \item \code{"Margalef"}: Margalef's Richness Index, calculated as \eqn{(S - 1) / \ln(N)}.
#'     \item \code{"Menhinick"}: Menhinick's Richness Index, calculated as \eqn{S / \sqrt{N}}.
#'   }
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame with two columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names of \code{data})}
#'   \item{SRI}{Species Richness Index for each sample}
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' cal.SRI(data = df, method = "Margalef")
#' @export
cal.SRI <- function(data, method = "Margalef", ...) {
  UseMethod("cal.SRI")
}

#' @rdname cal.SRI
#' @method cal.SRI data.frame
#' @exportS3Method Nematode::cal.SRI
cal.SRI.data.frame <- function(data, method = "Margalef", ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  if (!method %in% c("Margalef", "Menhinick")) {
    stop("Invalid method. Choose from 'Margalef' or 'Menhinick'.")
  }
  abundance <- data.frame(Sample.ID = rownames(data), abundance = rowSums(data, na.rm = T))
  N <- num_species(data)
  N_abundance <- merge(N, abundance, by = "Sample.ID")
  if (method == "Margalef") {
    SRI <- data.frame(
      Sample.ID = N_abundance$Sample.ID,
      SRI = (N_abundance$NumSpecies - 1) / log(N_abundance$abundance)
    )
  } else if (method == "Menhinick") {
    SRI <- data.frame(
      Sample.ID = N_abundance$Sample.ID,
      SRI = N_abundance$NumSpecies / sqrt(N_abundance$abundance)
    )
  }
  return(SRI)
}

#' @rdname cal.SRI
#' @method cal.SRI matrix
#' @exportS3Method Nematode::cal.SRI
cal.SRI.matrix <- function(data, method = "Margalef", ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }
  if (!method %in% c("Margalef", "Menhinick")) {
    stop("Invalid method. Choose from 'Margalef' or 'Menhinick'.")
  }
  abundance <- data.frame(Sample.ID = rownames(data), abundance = rowSums(data, na.rm = T))
  N <- num_species(data)
  N_abundance <- merge(N, abundance, by = "Sample.ID")
  if (method == "Margalef") {
    SRI <- data.frame(
      Sample.ID = N_abundance$Sample.ID,
      SRI = (N_abundance$NumSpecies - 1) / log(N_abundance$abundance)
    )
  } else if (method == "Menhinick") {
    SRI <- data.frame(
      Sample.ID = N_abundance$Sample.ID,
      SRI = N_abundance$NumSpecies / sqrt(N_abundance$abundance)
    )
  }
  return(SRI)
}

#' @rdname cal.SRI
#' @method cal.SRI default
#' @exportS3Method Nematode::cal.SRI
cal.SRI.default <- function(data, method = "Margalef", ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}

# =====Nematode Channel Ratio=====
#' Calculate Nematode Channel Ratio (NCR)
#'
#' @description
#' This function calculates the Nematode Channel Ratio (NCR) for ecological communities.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame with two columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names of \code{data})}
#'   \item{NCR}{Nematode Channel Ratio for each sample}
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' cal.NCR(data = df)
#' @export
cal.NCR <- function(data, ...) {
  UseMethod("cal.NCR")
}

#' @rdname cal.NCR
#' @method cal.NCR data.frame
#' @exportS3Method Nematode::cal.NCR
cal.NCR.data.frame <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  diet_rel <- diet_rel_abundance(data = data, relative = T)
  NCR <- data.frame(
    Sample.ID = diet_rel$Sample.ID,
    NCR = diet_rel$Ba / (diet_rel$Ba + diet_rel$Fu)
  )
  return(NCR)
}

#' @rdname cal.NCR
#' @method cal.NCR matrix
#' @exportS3Method Nematode::cal.NCR
cal.NCR.matrix <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }
  diet_rel <- diet_rel_abundance(data = data, relative = T)
  NCR <- data.frame(
    Sample.ID = diet_rel$Sample.ID,
    NCR = diet_rel$Ba / (diet_rel$Ba + diet_rel$Fu)
  )
  return(NCR)
}

#' @rdname cal.NCR
#' @method cal.NCR default
#' @exportS3Method Nematode::cal.NCR
cal.NCR.default <- function(data, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}


# =====Maturity Index====
#' Calculate Maturity Index (MI)
#'
#' @description
#' This function calculates the Maturity Index (MI) for ecological communities.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame with two columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names of \code{data})}
#'   \item{MI}{Maturity Index for each sample}
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' cal.MI(data = df)
#' @export
cal.MI <- function(data, ...) {
  UseMethod("cal.MI")
}

#' @rdname cal.MI
#' @method cal.MI data.frame
#' @exportS3Method Nematode::cal.MI
cal.MI.data.frame <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }

  genus_info <- check_nematode_genus(colnames(data), show.details = T)
  if (!all(genus_info$Exist)) {
    missing <- genus_info$Query.genus[!genus_info$Exist]
    stop(
      length(missing), " genus names not found: ", paste(head(missing), collapse = ", "),
      ifelse(length(missing) > 6, "...", "")
    )
  }

  valid_habits <- c("Bacterial feeders", "Fungus feeders", "Plant feeders", "Omnivores", "Predators")
  invalid_habits <- setdiff(unique(genus_info$Feeding_habit), valid_habits)
  if (length(invalid_habits) > 0) {
    stop("Unsupported feeding habit(s): ", paste(invalid_habits, collapse = ", "))
  }

  if (any(is.na(genus_info$CP_group))) {
    missing_cp <- genus_info$Query.genus[is.na(genus_info$CP_group)]
    stop(
      length(missing_cp), " genus lack CP values: ", paste(head(missing_cp), collapse = ", "),
      ifelse(length(missing_cp) > 6, "...", "")
    )
  }
  colnames(data) <- genus_info$Genus[match(colnames(data), genus_info$Genus)]
  data <- data %>% dplyr::select(dplyr::all_of(genus_info$Genus[genus_info$Feeding_habit != "Plant feeders"]))
  rel_data <- rel_abundance(data)
  cp_map <- stats::setNames(genus_info$CP_group, genus_info$Genus)
  cp_rel <- rel_data %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x * cp_map[dplyr::cur_column()]))
  MI <- data.frame(
    Sample.ID = rownames(cp_rel),
    MI = rowSums(cp_rel)
  )
  return(MI)
}

#' @rdname cal.MI
#' @method cal.MI matrix
#' @exportS3Method Nematode::cal.MI
cal.MI.matrix <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }

  genus_info <- check_nematode_genus(colnames(data), show.details = T)
  if (!all(genus_info$Exist)) {
    missing <- genus_info$Query.genus[!genus_info$Exist]
    stop(
      length(missing), " genus names not found: ", paste(head(missing), collapse = ", "),
      ifelse(length(missing) > 6, "...", "")
    )
  }

  valid_habits <- c("Bacterial feeders", "Fungus feeders", "Plant feeders", "Omnivores", "Predators")
  invalid_habits <- setdiff(unique(genus_info$Feeding_habit), valid_habits)
  if (length(invalid_habits) > 0) {
    stop("Unsupported feeding habit(s): ", paste(invalid_habits, collapse = ", "))
  }

  if (any(is.na(genus_info$CP_group))) {
    missing_cp <- genus_info$Query.genus[is.na(genus_info$CP_group)]
    stop(
      length(missing_cp), " genus lack CP values: ", paste(head(missing_cp), collapse = ", "),
      ifelse(length(missing_cp) > 6, "...", "")
    )
  }
  colnames(data) <- genus_info$Genus[match(colnames(data), genus_info$Genus)]
  sel_col <- genus_info$Genus[genus_info$Feeding_habit != "Plant feeders"]
  data <- data[, sel_col, drop = F]
  rel_data <- rel_abundance(data)
  cp_map <- stats::setNames(genus_info$CP_group, genus_info$Genus)
  cp_rel <- sweep(rel_data, 2, cp_map[colnames(rel_data)], "*")
  MI <- data.frame(
    Sample.ID = rownames(cp_rel),
    MI = rowSums(cp_rel)
  )
  return(MI)
}

#' @rdname cal.MI
#' @method cal.MI default
#' @exportS3Method Nematode::cal.MI
cal.MI.default <- function(data, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}


# =====Plant Parasite Index====
#' Calculate Plant Parasite Index (PPI)
#'
#' @description
#' This function calculates the Plant Parasite Index (PPI) for ecological communities.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame with two columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names of \code{data})}
#'   \item{PPI}{Plant Parasite Index for each sample}
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' cal.PPI(data = df)
#' @export
cal.PPI <- function(data, ...) {
  UseMethod("cal.PPI")
}

#' @rdname cal.PPI
#' @method cal.PPI data.frame
#' @exportS3Method Nematode::cal.PPI
cal.PPI.data.frame <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }

  genus_info <- check_nematode_genus(colnames(data), show.details = T)
  if (!all(genus_info$Exist)) {
    missing <- genus_info$Query.genus[!genus_info$Exist]
    stop(
      length(missing), " genus names not found: ", paste(head(missing), collapse = ", "),
      ifelse(length(missing) > 6, "...", "")
    )
  }

  valid_habits <- c("Bacterial feeders", "Fungus feeders", "Plant feeders", "Omnivores", "Predators")
  invalid_habits <- setdiff(unique(genus_info$Feeding_habit), valid_habits)
  if (length(invalid_habits) > 0) {
    stop("Unsupported feeding habit(s): ", paste(invalid_habits, collapse = ", "))
  }

  if (any(is.na(genus_info$CP_group))) {
    missing_cp <- genus_info$Query.genus[is.na(genus_info$CP_group)]
    stop(
      length(missing_cp), " genus lack CP values: ", paste(head(missing_cp), collapse = ", "),
      ifelse(length(missing_cp) > 6, "...", "")
    )
  }
  colnames(data) <- genus_info$Genus[match(colnames(data), genus_info$Genus)]
  data <- data %>% dplyr::select(dplyr::all_of(genus_info$Genus[genus_info$Feeding_habit == "Plant feeders"]))
  rel_data <- rel_abundance(data)
  cp_map <- stats::setNames(genus_info$CP_group, genus_info$Genus)
  cp_rel <- rel_data %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x * cp_map[dplyr::cur_column()]))
  PPI <- data.frame(
    Sample.ID = rownames(cp_rel),
    PPI = rowSums(cp_rel)
  )
  return(PPI)
}

#' @rdname cal.PPI
#' @method cal.PPI matrix
#' @exportS3Method Nematode::cal.PPI
cal.PPI.matrix <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }

  genus_info <- check_nematode_genus(colnames(data), show.details = T)
  if (!all(genus_info$Exist)) {
    missing <- genus_info$Query.genus[!genus_info$Exist]
    stop(
      length(missing), " genus names not found: ", paste(head(missing), collapse = ", "),
      ifelse(length(missing) > 6, "...", "")
    )
  }

  valid_habits <- c("Bacterial feeders", "Fungus feeders", "Plant feeders", "Omnivores", "Predators")
  invalid_habits <- setdiff(unique(genus_info$Feeding_habit), valid_habits)
  if (length(invalid_habits) > 0) {
    stop("Unsupported feeding habit(s): ", paste(invalid_habits, collapse = ", "))
  }

  if (any(is.na(genus_info$CP_group))) {
    missing_cp <- genus_info$Query.genus[is.na(genus_info$CP_group)]
    stop(
      length(missing_cp), " genus lack CP values: ", paste(head(missing_cp), collapse = ", "),
      ifelse(length(missing_cp) > 6, "...", "")
    )
  }
  colnames(data) <- genus_info$Genus[match(colnames(data), genus_info$Genus)]
  sel_col <- genus_info$Genus[genus_info$Feeding_habit == "Plant feeders"]
  data <- data[, sel_col, drop = F]
  rel_data <- rel_abundance(data)
  cp_map <- stats::setNames(genus_info$CP_group, genus_info$Genus)
  cp_rel <- sweep(rel_data, 2, cp_map[colnames(rel_data)], "*")
  PPI <- data.frame(
    Sample.ID = rownames(cp_rel),
    PPI = rowSums(cp_rel)
  )
  return(PPI)
}

#' @rdname cal.PPI
#' @method cal.PPI default
#' @exportS3Method Nematode::cal.PPI
cal.PPI.default <- function(data, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}


# =====Wasilewska Index=====
#' Calculate Wasilewska Index (WI)
#'
#' @description
#' This function calculates the Wasilewska Index (WI) for ecological communities.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame with two columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names of \code{data})}
#'   \item{WI}{Wasilewska Index for each sample}
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' cal.WI(data = df)
#' @export
cal.WI <- function(data, ...) {
  UseMethod("cal.WI")
}

#' @rdname cal.WI
#' @method cal.WI data.frame
#' @exportS3Method Nematode::cal.WI
cal.WI.data.frame <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  diet_rel <- diet_rel_abundance(data = data, total.abundance = NULL, relative = T)
  WI <- data.frame(
    Sample.ID = diet_rel$Sample.ID,
    WI = (diet_rel$Ba + diet_rel$Fu) / diet_rel$Pp
  )
  return(WI)
}

#' @rdname cal.WI
#' @method cal.WI matrix
#' @exportS3Method Nematode::cal.WI
cal.WI.matrix <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }
  diet_rel <- diet_rel_abundance(data = data, total.abundance = NULL, relative = T)
  WI <- data.frame(
    Sample.ID = diet_rel$Sample.ID,
    WI = (diet_rel$Ba + diet_rel$Fu) / diet_rel$Pp
  )
  return(WI)
}

#' @rdname cal.WI
#' @method cal.WI default
#' @exportS3Method Nematode::cal.WI
cal.WI.default <- function(data, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}


# =====Contributions=====
.contribution <- function(data, ...) {
  UseMethod(".contribution")
}

.contribution.data.frame <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  rel_data <- rel_abundance(data)
  genus_info <- check_nematode_genus(colnames(rel_data), show.details = T)
  if (!all(genus_info$Exist)) {
    missing <- genus_info$Query.genus[!genus_info$Exist]
    stop(
      length(missing), " genus names not found: ", paste(head(missing), collapse = ", "),
      ifelse(length(missing) > 6, "...", "")
    )
  }

  valid_habits <- c("Bacterial feeders", "Fungus feeders", "Plant feeders", "Omnivores", "Predators")
  invalid_habits <- setdiff(unique(genus_info$Feeding_habit), valid_habits)
  if (length(invalid_habits) > 0) {
    stop("Unsupported feeding habit(s): ", paste(invalid_habits, collapse = ", "))
  }

  if (any(is.na(genus_info$CP_group))) {
    missing_cp <- genus_info$Query.genus[is.na(genus_info$CP_group)]
    stop(
      length(missing_cp), " genus lack CP values: ", paste(head(missing_cp), collapse = ", "),
      ifelse(length(missing_cp) > 6, "...", "")
    )
  }
  colnames(rel_data) <- genus_info$Genus[match(colnames(rel_data), genus_info$Genus)]
  type_map <- c(
    "Bacterial feeders" = "Ba",
    "Fungus feeders" = "Fu",
    "Plant feeders" = "Pp",
    "Omnivores" = "Op",
    "Predators" = "Op"
  )
  genus_info$type <- type_map[genus_info$Feeding_habit]
  conditions <- list(
    Ba1 = genus_info$type == "Ba" & genus_info$CP_group == 1,
    Ba2 = genus_info$type == "Ba" & genus_info$CP_group == 2,
    Ba3 = genus_info$type == "Ba" & genus_info$CP_group == 3,
    Ba4 = genus_info$type == "Ba" & genus_info$CP_group == 4,
    Fu2 = genus_info$type == "Fu" & genus_info$CP_group == 2,
    Fu3 = genus_info$type == "Fu" & genus_info$CP_group == 3,
    Fu4 = genus_info$type == "Fu" & genus_info$CP_group == 4,
    Op4 = genus_info$type == "Op" & genus_info$CP_group == 4,
    Op5 = genus_info$type == "Op" & genus_info$CP_group == 5
  )
  cal_fun <- function(data, genus_info, condition) {
    sel_genus <- genus_info %>%
      dplyr::filter(!!condition) %>%
      dplyr::pull("Genus")
    data %>%
      dplyr::select(dplyr::all_of(sel_genus)) %>%
      rowSums()
  }
  rel_sum <- lapply(names(conditions), function(names) {
    cal_fun(rel_data, genus_info, conditions[[names]])
  })
  rel_sum <- as.data.frame(rel_sum)
  colnames(rel_sum) <- names(conditions)
  coefs <- c(
    Ba1 = 3.2,
    Ba2 = 0.8,
    Ba3 = 1.8,
    Ba4 = 3.2,
    Fu2 = 0.8,
    Fu3 = 1.8,
    Fu4 = 3.2,
    Op4 = 3.2,
    Op5 = 5
  )
  contributions <- rel_sum %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x * coefs[dplyr::cur_column()]))
  return(contributions)
}

.contribution.matrix <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }
  rel_data <- rel_abundance(data)
  genus_info <- check_nematode_genus(colnames(rel_data), show.details = T)
  if (!all(genus_info$Exist)) {
    missing <- genus_info$Query.genus[!genus_info$Exist]
    stop(
      length(missing), " genus names not found: ", paste(head(missing), collapse = ", "),
      ifelse(length(missing) > 6, "...", "")
    )
  }

  valid_habits <- c("Bacterial feeders", "Fungus feeders", "Plant feeders", "Omnivores", "Predators")
  invalid_habits <- setdiff(unique(genus_info$Feeding_habit), valid_habits)
  if (length(invalid_habits) > 0) {
    stop("Unsupported feeding habit(s): ", paste(invalid_habits, collapse = ", "))
  }

  if (any(is.na(genus_info$CP_group))) {
    missing_cp <- genus_info$Query.genus[is.na(genus_info$CP_group)]
    stop(
      length(missing_cp), " genus lack CP values: ", paste(head(missing_cp), collapse = ", "),
      ifelse(length(missing_cp) > 6, "...", "")
    )
  }
  colnames(rel_data) <- genus_info$Genus[match(colnames(rel_data), genus_info$Genus)]
  type_map <- c(
    "Bacterial feeders" = "Ba",
    "Fungus feeders" = "Fu",
    "Plant feeders" = "Pp",
    "Omnivores" = "Op",
    "Predators" = "Op"
  )
  genus_info$type <- type_map[genus_info$Feeding_habit]
  conditions <- list(
    Ba1 = genus_info$type == "Ba" & genus_info$CP_group == 1,
    Ba2 = genus_info$type == "Ba" & genus_info$CP_group == 2,
    Ba3 = genus_info$type == "Ba" & genus_info$CP_group == 3,
    Ba4 = genus_info$type == "Ba" & genus_info$CP_group == 4,
    Fu2 = genus_info$type == "Fu" & genus_info$CP_group == 2,
    Fu3 = genus_info$type == "Fu" & genus_info$CP_group == 3,
    Fu4 = genus_info$type == "Fu" & genus_info$CP_group == 4,
    Op4 = genus_info$type == "Op" & genus_info$CP_group == 4,
    Op5 = genus_info$type == "Op" & genus_info$CP_group == 5
  )

  rel_sum_list <- lapply(names(conditions), function(t) {
    selected_cols <- colnames(rel_data) %in% genus_info$Genus[conditions[[t]]]
    if (sum(selected_cols) == 0) {
      stats::setNames(numeric(nrow(rel_data)), rownames(rel_data))
    } else {
      rowSums(rel_data[, selected_cols, drop = FALSE])
    }
  })

  rel_sum <- as.data.frame(do.call(cbind, rel_sum_list))
  colnames(rel_sum) <- names(conditions)
  rownames(rel_sum) <- rownames(rel_data)

  coefs <- c(
    Ba1 = 3.2,
    Ba2 = 0.8,
    Ba3 = 1.8,
    Ba4 = 3.2,
    Fu2 = 0.8,
    Fu3 = 1.8,
    Fu4 = 3.2,
    Op4 = 3.2,
    Op5 = 5
  )
  contributions <- rel_sum %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x * coefs[dplyr::cur_column()]))
  return(contributions)
}

.contribution.default <- function(data, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}

# =====Channel Index=====
#' Calculate Channel Index (CI)
#'
#' @description
#' This function calculates the Channel Index (CI) for ecological communities.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame with two columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names of \code{data})}
#'   \item{CI}{Channel Index for each sample}
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' cal.CI(data = df)
#' @export
cal.CI <- function(data, ...) {
  UseMethod("cal.CI")
}

#' @rdname cal.CI
#' @method cal.CI data.frame
#' @exportS3Method Nematode::cal.CI
cal.CI.data.frame <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  contributions <- .contribution(data)
  CI <- data.frame(
    Sample.ID = row.names(contributions),
    CI = contributions$Fu2 / (contributions$Ba1 + contributions$Fu2) * 100
  )
  return(CI)
}

#' @rdname cal.CI
#' @method cal.CI matrix
#' @exportS3Method Nematode::cal.CI
cal.CI.matrix <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }
  contributions <- .contribution(data)
  CI <- data.frame(
    Sample.ID = row.names(contributions),
    CI = contributions$Fu2 / (contributions$Ba1 + contributions$Fu2) * 100
  )
  return(CI)
}

#' @rdname cal.CI
#' @method cal.CI default
#' @exportS3Method Nematode::cal.CI
cal.CI.default <- function(data, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}


# =====Enrichment Index=====
#' Calculate Enrichment Index (EI)
#'
#' @description
#' This function calculates the Enrichment Index (EI) for ecological communities.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame with two columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names of \code{data})}
#'   \item{EI}{Enrichment Index for each sample}
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' cal.EI(data = df)
#' @export
cal.EI <- function(data, ...) {
  UseMethod("cal.EI")
}

#' @rdname cal.EI
#' @method cal.EI data.frame
#' @exportS3Method Nematode::cal.EI
cal.EI.data.frame <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  contributions <- .contribution(data)
  b <- contributions %>%
    dplyr::select(dplyr::all_of(c("Ba2", "Fu2"))) %>%
    rowSums()
  e <- contributions %>%
    dplyr::select(dplyr::all_of(c("Ba1", "Fu2"))) %>%
    rowSums()
  Sample.ID <- row.names(contributions)
  b <- b[Sample.ID]
  e <- e[Sample.ID]
  EI <- data.frame(
    Sample.ID = Sample.ID,
    EI = e / (b + e) * 100
  )
  return(EI)
}

#' @rdname cal.EI
#' @method cal.EI matrix
#' @exportS3Method Nematode::cal.EI
cal.EI.matrix <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }
  contributions <- .contribution(data)
  b <- contributions %>%
    dplyr::select(dplyr::all_of(c("Ba2", "Fu2"))) %>%
    rowSums()
  e <- contributions %>%
    dplyr::select(dplyr::all_of(c("Ba1", "Fu2"))) %>%
    rowSums()
  Sample.ID <- row.names(contributions)
  b <- b[Sample.ID]
  e <- e[Sample.ID]
  EI <- data.frame(
    Sample.ID = Sample.ID,
    EI = e / (b + e) * 100
  )
  return(EI)
}

#' @rdname cal.EI
#' @method cal.EI default
#' @exportS3Method Nematode::cal.EI
cal.EI.default <- function(data, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}


# =====Structure Index=====
#' Calculate Structure Index (SI)
#'
#' @description
#' This function calculates the Structure Index (SI) for ecological communities.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame with two columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names of \code{data})}
#'   \item{SI}{Structure Index for each sample}
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' cal.SI(data = df)
#' @export
cal.SI <- function(data, ...) {
  UseMethod("cal.SI")
}

#' @rdname cal.SI
#' @method cal.SI data.frame
#' @exportS3Method Nematode::cal.SI
cal.SI.data.frame <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  contributions <- .contribution(data)
  b <- contributions %>%
    dplyr::select(dplyr::all_of(c("Ba2", "Fu2"))) %>%
    rowSums()
  s <- contributions %>%
    dplyr::select(dplyr::all_of(c("Ba3", "Ba4", "Fu3", "Fu4", "Op4", "Op5"))) %>%
    rowSums()
  Sample.ID <- row.names(contributions)
  b <- b[Sample.ID]
  s <- s[Sample.ID]
  SI <- data.frame(
    Sample.ID = Sample.ID,
    SI = s / (b + s) * 100
  )
  return(SI)
}

#' @rdname cal.SI
#' @method cal.SI matrix
#' @exportS3Method Nematode::cal.SI
cal.SI.matrix <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }
  contributions <- .contribution(data)
  b <- contributions %>%
    dplyr::select(dplyr::all_of(c("Ba2", "Fu2"))) %>%
    rowSums()
  s <- contributions %>%
    dplyr::select(dplyr::all_of(c("Ba3", "Ba4", "Fu3", "Fu4", "Op4", "Op5"))) %>%
    rowSums()
  Sample.ID <- row.names(contributions)
  b <- b[Sample.ID]
  s <- s[Sample.ID]
  SI <- data.frame(
    Sample.ID = Sample.ID,
    SI = s / (b + s) * 100
  )
  return(SI)
}

#' @rdname cal.SI
#' @method cal.SI default
#' @exportS3Method Nematode::cal.SI
cal.SI.default <- function(data, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}


# =====Basic Index=====
#' Calculate Basic Index (BI)
#'
#' @description
#' This function calculates the Basic Index (BI) for ecological communities.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame with two columns:
#'   \item{Sample.ID}{Character vector of sample identifiers (from row names of \code{data})}
#'   \item{BI}{Basic Index for each sample}
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' cal.BI(data = df)
#' @export
cal.BI <- function(data, ...) {
  UseMethod("cal.BI")
}

#' @rdname cal.BI
#' @method cal.BI data.frame
#' @exportS3Method Nematode::cal.BI
cal.BI.data.frame <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }
  contributions <- .contribution(data)
  b <- contributions %>%
    dplyr::select(dplyr::all_of(c("Ba2", "Fu2"))) %>%
    rowSums()
  s <- contributions %>%
    dplyr::select(dplyr::all_of(c("Ba3", "Ba4", "Fu3", "Fu4", "Op4", "Op5"))) %>%
    rowSums()
  e <- contributions %>%
    dplyr::select(dplyr::all_of(c("Ba1", "Fu2"))) %>%
    rowSums()
  Sample.ID <- row.names(contributions)
  b <- b[Sample.ID]
  s <- s[Sample.ID]
  e <- e[Sample.ID]
  BI <- data.frame(
    Sample.ID = Sample.ID,
    BI = b / (b + s + e) * 100
  )
  return(BI)
}

#' @rdname cal.BI
#' @method cal.BI matrix
#' @exportS3Method Nematode::cal.BI
cal.BI.matrix <- function(data, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }
  contributions <- .contribution(data)
  b <- contributions %>%
    dplyr::select(dplyr::all_of(c("Ba2", "Fu2"))) %>%
    rowSums()
  s <- contributions %>%
    dplyr::select(dplyr::all_of(c("Ba3", "Ba4", "Fu3", "Fu4", "Op4", "Op5"))) %>%
    rowSums()
  e <- contributions %>%
    dplyr::select(dplyr::all_of(c("Ba1", "Fu2"))) %>%
    rowSums()
  Sample.ID <- row.names(contributions)
  b <- b[Sample.ID]
  s <- s[Sample.ID]
  e <- e[Sample.ID]
  BI <- data.frame(
    Sample.ID = Sample.ID,
    BI = b / (b + s + e) * 100
  )
  return(BI)
}

#' @rdname cal.BI
#' @method cal.BI default
#' @exportS3Method Nematode::cal.BI
cal.BI.default <- function(data, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}


# =====Ecological Indices=====
#' Calculate Ecological Indices of Nematodes
#'
#' @description
#' This function calculates various ecological indices based on the provided nematode genus abundance data.
#' It supports a range of indices, including taxonomic diversity, Shannon diversity index, Pielou's evenness index,
#' Simpson's index, and more. Users can specify which indices to calculate or use the default option to calculate all supported indices.
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param indices A character vector specifying the ecological indices to be calculated.
#'   The following indices are supported:
#'   \itemize{
#'     \item "TD" - Trophic Diversity
#'     \item "H" - Shannon-Wiener Index
#'     \item "J" - Pielou's Evenness Index
#'     \item "Simpson" - Simpson Index
#'     \item "WI" - Wasilewska Index
#'     \item "MI" - Maturity Index
#'     \item "PPI" - Plant Parasite Index
#'     \item "SRI" - Species Richness Index
#'     \item "NCR" - Nematode Channel Ratio
#'     \item "CI" - Channel Index
#'     \item "BI" - Basic Index
#'     \item "EI" - Enrichment Index
#'     \item "SI" - Structure Index
#'   }
#'   Additionally, specifying \code{All} will calculate all supported indices. \code{All} is the default value.
#' @param method The method to use for calculating the Species Richness Index. Default is \code{NULL}, which uses the default method \code{Margalef}. Options include:
#'   \itemize{
#'     \item \code{"Margalef"}: Margalef's Richness Index, calculated as \eqn{(S - 1) / \ln(N)}.
#'     \item \code{"Menhinick"}: Menhinick's Richness Index, calculated as \eqn{S / \sqrt{N}}.
#'   }
#' @param ... Additional arguments (currently unused).
#'
#' @returns A data frame containing the calculated indices. The data frame includes a \code{Sample.ID} column and additional columns for each requested index.
#'
#' @examples
#' # Example with a data frame
#' df <- data.frame(
#'   Cephalobus = c(10, NA, 15),
#'   Caenorhabditis = c(5, 10, NA),
#'   Pratylenchus = c(8, 12, 10),
#'   row.names = c("A", "B", "C")
#' )
#' Ecological.Indices(data = df, indices = "All", method = "Menhinick")
#' @export
Ecological.Indices <- function(data, indices = "All", method = NULL, ...) {
  UseMethod("Ecological.Indices")
}

#' @rdname Ecological.Indices
#' @method Ecological.Indices data.frame
#' @exportS3Method Nematode::Ecological.Indices
Ecological.Indices.data.frame <- function(data, indices = "All", method = NULL, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }

  support.indices <- c(
    "TD", "H", "J", "Simpson", "WI", "MI", "PPI",
    "SRI", "NCR", "CI", "BI", "EI", "SI"
  )
  unsupported_indices <- setdiff(indices, c(support.indices, "All"))

  if (length(unsupported_indices) > 0) {
    stop(paste(
      "The following indices are not supported:",
      paste(unsupported_indices, collapse = ", ")
    ))
  }

  if ("All" %in% indices) {
    indices <- support.indices
  }

  result_list <- list()

  if ("SRI" %in% indices) {
    if (is.null(method)) {
      warning("Using default method 'Margalef' for 'SRI'. To use 'Menhinick' method, specify method = 'Menhinick'.")
      method <- "Margalef"
    }
    result_list[["SRI"]] <- cal.SRI(data = data, method = method)
    indices <- setdiff(indices, "SRI")
  }

  if (length(indices) > 0) {
    result_list <- c(
      result_list,
      lapply(indices, function(index) {
        func <- match.fun(paste0("cal.", index))
        func(data)
      })
    )
  }

  result_df <- purrr::reduce(result_list, dplyr::left_join, by = "Sample.ID")
  result_df <- result_df %>% dplyr::select(dplyr::any_of(c("Sample.ID", support.indices)))
  return(result_df)
}

#' @rdname Ecological.Indices
#' @method Ecological.Indices matrix
#' @exportS3Method Nematode::Ecological.Indices
Ecological.Indices.matrix <- function(data, indices = "All", method = NULL, ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }
  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }

  support.indices <- c(
    "TD", "H", "J", "Simpson", "WI", "MI", "PPI",
    "SRI", "NCR", "CI", "BI", "EI", "SI"
  )
  unsupported_indices <- setdiff(indices, c(support.indices, "All"))

  if (length(unsupported_indices) > 0) {
    stop(paste(
      "The following indices are not supported:",
      paste(unsupported_indices, collapse = ", ")
    ))
  }

  if ("All" %in% indices) {
    indices <- support.indices
  }

  result_list <- list()

  if ("SRI" %in% indices) {
    if (is.null(method)) {
      warning("Using default method 'Margalef' for 'SRI'. To use 'Menhinick' method, specify method = 'Menhinick'.")
      method <- "Margalef"
    }
    result_list[["SRI"]] <- cal.SRI(data = data, method = method)
    indices <- setdiff(indices, "SRI")
  }

  if (length(indices) > 0) {
    result_list <- c(
      result_list,
      lapply(indices, function(index) {
        func <- match.fun(paste0("cal.", index))
        func(data)
      })
    )
  }

  result_df <- purrr::reduce(result_list, dplyr::left_join, by = "Sample.ID")
  result_df <- result_df %>% dplyr::select(dplyr::any_of(c("Sample.ID", support.indices)))
  return(result_df)
}

#' @rdname Ecological.Indices
#' @method Ecological.Indices default
#' @exportS3Method Nematode::Ecological.Indices
Ecological.Indices.default <- function(data, indices = "All", method = NULL, ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}
