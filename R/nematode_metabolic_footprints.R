#' @importFrom dplyr %>% mutate across everything select all_of any_of filter pull
#' @importFrom stats setNames
#' @importFrom utils data head
NULL
#' Nematode Metabolic Footprints (NMF) Calculation
#'
#' @description
#' This function calculates various Nematode Metabolic Footprints (NMF) based on the input data and abundance information.
#' It supports multiple types of NMF calculations and can handle data in different formats (data.frame or matrix).
#'
#' @param data \code{data.frame} or \code{matrix}. The nematode abundance table where rows represent samples and columns represent nematode genera.
#' Each element indicates the count of a specific nematode genus in the corresponding sample. Row names must be sample names, and column names must be nematode genus names.
#' @param abundance \code{data.frame}. A data frame with sample names as row names and a single column containing the total nematode abundance for each sample.
#' @param type Character vector specifying the type(s) of NMF to calculate. \code{All} is the default value. Valid options include:
#'   \itemize{
#'     \item "BaMF" - Bacterial Feeders Nematode Metabolic Footprints
#'     \item "FuMF" - Fungal Feeders Nematode Metabolic Footprints
#'     \item "PpMF" - Plant Feeders Nematode Metabolic Footprints
#'     \item "OpMF" - Omnivore/Predator Nematode Metabolic Footprints
#'     \item "Fe" - Enrichment Footprints (CP group <= 2)
#'     \item "Fs" - Structure Footprints (CP group > 2)
#'     \item "TNMF" - Total Nematode Metabolic Footprints
#'     \item "FMF" - Functional Metabolic Footprints (product of Fe and Fs)
#'     \item "All" - Calculate all the above types
#'   }
#' @param ... Additional arguments (currently unused).
#'
#' @return A data.frame containing the calculated NMF values for each sample. The columns represent different NMF types, and the rows correspond to samples.
#'
#' @examples
#' data <- data.frame(
#'   Cephalobus = c(10, 20, 30),
#'   Eucephalobus = c(5, 10, 12),
#'   Acrobeloides = c(1, 2, 3),
#'   Caenorhabditis = c(5, 8, 15),
#'   Aphelenchus = c(5, 13, 11),
#'   Leptonchus = c(3, 10, 15),
#'   Pratylenchus = c(9, 2, 15),
#'   Tylenchus = c(5, 0, 15),
#'   Mesodorylaimus = c(7, 10, 18),
#'   Discolaimus = c(1, 10, 25),
#'   row.names = c("Sample1", "Sample2", "Sample3")
#' )
#' abundance <- data.frame(
#'   Abundance = c(100, 200, 300),
#'   row.names = c("Sample1", "Sample2", "Sample3")
#' )
#' result <- NMF(data, abundance, type = "All")
#' print(result)
#' @export
NMF <- function(data, abundance, type = "All", ...) {
  UseMethod("NMF")
}

#' @rdname NMF
#' @method NMF data.frame
#' @exportS3Method Nematode::NMF
NMF.data.frame <- function(data, abundance, type = "All", ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }

  if (!all(row.names(data) %in% row.names(abundance))) {
    stop("Please provide the abundance information for all samples!")
  }

  data[is.na(data)] <- 0

  if (!all(sapply(data, is.numeric))) {
    stop("The data contains non-numeric characters!")
  }

  if (!all(sapply(abundance, is.numeric))) {
    stop("The abundance contains non-numeric characters!")
  }

  valid_types <- c("BaMF", "FuMF", "PpMF", "OpMF", "Fe", "Fs", "TNMF", "FMF", "All")
  invalid_types <- setdiff(type, valid_types)
  if (length(invalid_types) > 0) {
    stop("Invalid type(s): ", paste(invalid_types, collapse = ", "),
         "\nValid types are: ", paste(valid_types, collapse = ", "))
  }

  genus_info <- check_nematode_genus(colnames(data))

  if (!all(genus_info$Exist)) {
    missing <- genus_info$Query.genus[!genus_info$Exist]
    stop(length(missing), " genus names not found: ", paste(head(missing), collapse = ", "),
         ifelse(length(missing) > 6, "...", ""))
  }

  valid_habits <- c("Bacterial feeders", "Fungus feeders", "Plant feeders", "Omnivores", "Predators")
  invalid_habits <- setdiff(unique(genus_info$Feeding_habit), valid_habits)
  if (length(invalid_habits) > 0) {
    stop("Unsupported feeding habit(s): ", paste(invalid_habits, collapse = ", "))
  }

  if (any(is.na(genus_info$CP_group))) {
    missing_cp <- genus_info$Query.genus[is.na(genus_info$CP_group)]
    stop(length(missing_cp), " genus lack CP values: ", paste(head(missing_cp), collapse = ", "),
         ifelse(length(missing_cp) > 6, "...", ""))
  }

  nematode.ave.mass <- Nematode::nematode.ave.mass
  colnames(data) <- genus_info$Genus[match(colnames(data), genus_info$Genus)]
  genus_info$mass <- nematode.ave.mass$Genus.Average.Mass[match(genus_info$Genus, nematode.ave.mass$Genus)]
  genus_info[is.na(genus_info$mass), ]$mass <- nematode.ave.mass$Family.Average.Mass[match(genus_info[is.na(genus_info$mass), ]$Family, nematode.ave.mass$Family)]

  if (any(is.na(genus_info$mass))) {
    no_mass_genus_str <- paste0("*", genus_info$Query.genus[is.na(genus_info$mass)], collapse = ", ")
    stop(paste("The following genus lack mass values:", no_mass_genus_str))
  }

  type_map <- c(
    "Bacterial feeders" = "BaMF",
    "Fungus feeders" = "FuMF",
    "Plant feeders" = "PpMF",
    "Omnivores" = "OpMF",
    "Predators" = "OpMF"
  )
  genus_info$type <- type_map[genus_info$Feeding_habit]
  genus_info$MF <- (0.1 * (genus_info$mass / genus_info$CP_group)) + (0.273 * (genus_info$mass^0.75))
  total_ident <- rowSums(data)
  abundance_value <- abundance[match(row.names(data), row.names(abundance)), colnames(abundance)[1]]
  rel_data <- data %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ (.x / total_ident) * abundance_value))
  mul_map <- stats::setNames(genus_info$MF, genus_info$Genus)
  MF_df <- rel_data %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x * mul_map[dplyr::cur_column()]))

  conditions <- list(
    BaMF = genus_info$type == "BaMF",
    FuMF = genus_info$type == "FuMF",
    PpMF = genus_info$type == "PpMF",
    OpMF = genus_info$type == "OpMF",
    TNMF = rep(TRUE, nrow(genus_info)),
    Fe = genus_info$CP_group <= 2,
    Fs = genus_info$CP_group > 2
  )

  cal_fun <- function(data, genus_info, condition) {
    sel_genus <- genus_info %>%
      dplyr::filter(!!condition) %>%
      dplyr::pull('Genus')
    data %>%
      dplyr::select(dplyr::all_of(sel_genus)) %>%
      rowSums()
  }

  adj_type <- if ("All" %in% type) {
    c("BaMF", "FuMF", "PpMF", "OpMF", "Fe", "Fs", "TNMF")
  } else if ("FMF" %in% type) {
    setdiff(unique(c(type, "Fe", "Fs")), 'FMF')
  } else {
    type
  }

  result <- lapply(adj_type, function(names) {
    cal_fun(MF_df, genus_info, conditions[[names]])
  })
  result_df <- as.data.frame(result)
  colnames(result_df) <- adj_type

  if ("All" %in% type) {
    result_df$FMF <- result_df$Fe * result_df$Fs
  } else if ("FMF" %in% type) {
    result_df$FMF <- result_df$Fe * result_df$Fs
    result_df <- result_df %>% dplyr::select(dplyr::all_of(type))
  }
  col_order <- c("BaMF", "FuMF", "PpMF", "OpMF", "Fe", "Fs", "TNMF", "FMF")
  result_df <- result_df %>% dplyr::select(dplyr::any_of(col_order))
  result_df <- data.frame(Sample.ID = row.names(result_df), result_df)
  return(result_df)
}


#' @rdname NMF
#' @method NMF matrix
#' @exportS3Method Nematode::NMF
NMF.matrix <- function(data, abundance, type = "All", ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }

  if (!all(rownames(data) %in% rownames(abundance))) {
    stop("Abundance data missing for some samples in the input matrix")
  }

  data[is.na(data)] <- 0

  if (!is.numeric(data)) stop("Input matrix contains non-numeric values")
  if (!all(vapply(abundance, is.numeric, logical(1)))) {
    stop("Abundance data contains non-numeric columns")
  }

  valid_types <- c("BaMF", "FuMF", "PpMF", "OpMF", "Fe", "Fs", "TNMF", "FMF", "All")
  invalid_types <- setdiff(type, valid_types)
  if (length(invalid_types) > 0) {
    stop("Invalid type(s): ", paste(invalid_types, collapse = ", "),
         "\nValid types are: ", paste(valid_types, collapse = ", "))
  }

  genus_info <- check_nematode_genus(colnames(data))

  if (!all(genus_info$Exist)) {
    missing <- genus_info$Query.genus[!genus_info$Exist]
    stop(length(missing), " genus names not found: ", paste(head(missing), collapse = ", "),
         ifelse(length(missing) > 6, "...", ""))
  }

  valid_habits <- c("Bacterial feeders", "Fungus feeders", "Plant feeders", "Omnivores", "Predators")
  invalid_habits <- setdiff(unique(genus_info$Feeding_habit), valid_habits)
  if (length(invalid_habits) > 0) {
    stop("Unsupported feeding habit(s): ", paste(invalid_habits, collapse = ", "))
  }

  if (any(is.na(genus_info$CP_group))) {
    missing_cp <- genus_info$Query.genus[is.na(genus_info$CP_group)]
    stop(length(missing_cp), " genus lack CP values: ", paste(head(missing_cp), collapse = ", "),
         ifelse(length(missing_cp) > 6, "...", ""))
  }

  nematode.ave.mass <- Nematode::nematode.ave.mass
  colnames(data) <- genus_info$Genus[match(colnames(data), genus_info$Query.genus)]

  genus_mass <- nematode.ave.mass$Genus.Average.Mass[match(genus_info$Genus, nematode.ave.mass$Genus)]
  family_mass <- nematode.ave.mass$Family.Average.Mass[match(genus_info$Family, nematode.ave.mass$Family)]
  genus_info$mass <- ifelse(is.na(genus_mass), family_mass, genus_mass)

  if (any(is.na(genus_info$mass))) {
    missing_mass <- genus_info$Query.genus[is.na(genus_info$mass)]
    stop(length(missing_mass), " genus lack mass values: ", paste(head(missing_mass), collapse = ", "),
         ifelse(length(missing_mass) > 6, "...", ""))
  }

  MF_params <- (0.1 * (genus_info$mass / genus_info$CP_group)) + (0.273 * (genus_info$mass^0.75))
  names(MF_params) <- genus_info$Genus

  total_ident <- rowSums(data)
  abundance_vec <- abundance[match(rownames(data), rownames(abundance)), 1]
  rel_data <- (data / total_ident) * abundance_vec

  MF_matrix <- t(t(rel_data) * MF_params[colnames(rel_data)])

  type_map <- c(
    "Bacterial feeders" = "BaMF",
    "Fungus feeders" = "FuMF",
    "Plant feeders" = "PpMF",
    "Omnivores" = "OpMF",
    "Predators" = "OpMF"
  )
  genus_info$type <- type_map[genus_info$Feeding_habit]

  conditions <- list(
    BaMF = genus_info$type == "BaMF",
    FuMF = genus_info$type == "FuMF",
    PpMF = genus_info$type == "PpMF",
    OpMF = genus_info$type == "OpMF",
    TNMF = rep(TRUE, nrow(genus_info)),
    Fe = genus_info$CP_group <= 2,
    Fs = genus_info$CP_group > 2
  )

  adj_type <- if ("All" %in% type) {
    c("BaMF", "FuMF", "PpMF", "OpMF", "Fe", "Fs", "TNMF")
  } else if ("FMF" %in% type) {
    setdiff(unique(c(type, "Fe", "Fs")), 'FMF')
  } else {
    type
  }

  result_list <- lapply(adj_type, function(t) {
    selected_cols <- colnames(MF_matrix) %in% genus_info$Genus[conditions[[t]]]
    if (sum(selected_cols) == 0) {
      stats::setNames(numeric(nrow(MF_matrix)), rownames(MF_matrix))
    } else {
      rowSums(MF_matrix[, selected_cols, drop = FALSE])
    }
  })

  result_df <- as.data.frame(do.call(cbind, result_list))
  colnames(result_df) <- adj_type
  rownames(result_df) <- rownames(data)

  if ("All" %in% type) {
    result_df$FMF <- result_df$Fe * result_df$Fs
  } else if ("FMF" %in% type) {
    result_df$FMF <- result_df$Fe * result_df$Fs
    result_df <- result_df %>% dplyr::select(dplyr::all_of(type))
  }

  col_order <- c("BaMF", "FuMF", "PpMF", "OpMF", "Fe", "Fs", "TNMF", "FMF")
  result_df <- result_df[, intersect(col_order, colnames(result_df)), drop = FALSE]

  data.frame(Sample.ID = rownames(data), result_df, row.names = NULL, check.names = FALSE)
}

#' @rdname NMF
#' @method NMF default
#' @exportS3Method Nematode::NMF
NMF.default <- function(data, abundance, type = "all", ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}
