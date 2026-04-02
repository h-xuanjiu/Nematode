#' @importFrom dplyr %>% mutate across everything select all_of any_of filter pull
#' @importFrom stats setNames sd
#' @importFrom utils data head
NULL

#' Nematode Energy Footprints (NEF) Calculation
#'
#' @param data A data.frame or matrix containing nematode genus abundance data. Rows represent samples, and columns represent genera.
#' @param abundance A data.frame containing abundance information for the samples. It must match the row names of the input data.
#' @param AE A named list specifying the assimilation efficiencies for nematode feeding groups.
#' Must contain the following elements:
#' \itemize{
#'   \item Ba - Assimilation efficiency for bacterial feeders (default: 0.6)
#'   \item Fu - Assimilation efficiency for fungal feeders (default: 0.38)
#'   \item Pp - Assimilation efficiency for plant feeders (default: 0.25)
#'   \item Op - Assimilation efficiency for omnivores/predators (default: 0.5)
#' }
#' @param ... Additional arguments (currently unused).
#'
#' @return A list object of class `"NEF"` containing the following components:
#' \describe{
#'   \item{data}{A list with original input data:
#'     \itemize{
#'       \item data - Original genus abundance data.frame or matrix of nematode genera
#'       \item Abundance - Total abundance data used for calculations
#'     }
#'   }
#'   \item{Energy.flux}{A list containing energy flow calculations:
#'     \itemize{
#'       \item Energy.flux: Data frame of energy flows (\eqn{\mu g~C~100g^{-1}~dry~soil}) per feeding group. Columns:
#'         \itemize{
#'           \item Sample.ID - Sample identifier
#'           \item BaEF - Bacterial feeders energy flows
#'           \item FuEF - Fungal feeders energy flows
#'           \item PpEF - Plant feeders energy flows
#'           \item OpEF - Omnivores/Predators energy flows
#'           \item TNEF - Total energy flows of nematodes
#'         }
#'       \item C.flux.node: Data frame of Biomass (\eqn{\mu g~C~100g^{-1}~dry~soil}) per feeding group. Columns:
#'         \itemize{
#'           \item Sample.ID - Sample identifier
#'           \item Ba - Bacterial feeders biomass
#'           \item Fu - Fungal feeders biomass
#'           \item Pp - Plant feeders biomass
#'           \item Op - Omnivores/Predators biomass
#'         }
#'       \item C.flux.path: Data frame of energy flows (\eqn{\mu g~C~100g^{-1}~dry~soil~day^{-1}}). Columns:
#'         \itemize{
#'           \item Sample.ID - Sample identifier
#'           \item R.to.Ba - Carbon flux from Resources to bacterial feeders
#'           \item R.to.Fu - Carbon flux from Resources to fungal feeders
#'           \item R.to.Pp - Carbon flux from Resources to plant feeders
#'           \item Ba.to.Op - Carbon flux from bacterial to omnivorous channels
#'           \item Fu.to.Op - Carbon flux from fungal to omnivorous channels
#'           \item Pp.to.Op - Carbon flux from plant to omnivorous channels
#'         }
#'       \item U: Data frame of ecosystem stability indices. Columns:
#'         \itemize{
#'           \item Sample.ID - Sample identifier
#'           \item U - Energy flow uniformity index
#'         }
#'     }
#'   }
#' }
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
#' result <- NEF(data, abundance)
#' print(result)
#' @export
NEF <- function(data, abundance, AE = list(Ba = 0.6, Fu = 0.38, Pp = 0.25, Op = 0.5), ...) {
  UseMethod("NEF")
}

#' @rdname NEF
#' @method NEF data.frame
#' @exportS3Method Nematode::NEF
NEF.data.frame <- function(data, abundance, AE = list(Ba = 0.6, Fu = 0.38, Pp = 0.25, Op = 0.5), ...) {
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

  if (any(is.na(genus_info$CP_group))) {
    missing_cp <- genus_info$Query.genus[is.na(genus_info$CP_group)]
    stop(
      length(missing_cp), " genus lack CP values: ", paste(head(missing_cp), collapse = ", "),
      ifelse(length(missing_cp) > 6, "...", "")
    )
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
    "Bacterial feeders" = "BaEF",
    "Fungus feeders" = "FuEF",
    "Plant feeders" = "PpEF",
    "Omnivores" = "OpEF",
    "Predators" = "OpEF"
  )
  genus_info$type <- type_map[genus_info$Feeding_habit]
  genus_info$EF <- (0.1 * (genus_info$mass / (genus_info$CP_group * 12))) + (0.0159 * (genus_info$mass^0.75))
  total_ident <- rowSums(data)
  abundance_value <- abundance[match(row.names(data), row.names(abundance)), colnames(abundance)[1]]
  rel_data <- data %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ (.x / total_ident) * abundance_value))
  mul_map_EF <- stats::setNames(genus_info$EF, genus_info$Genus)
  EF_df <- rel_data %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x * mul_map_EF[dplyr::cur_column()]))
  mul_map_biomass <- stats::setNames(genus_info$mass, genus_info$Genus)
  biomass_df <- rel_data %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x * mul_map_biomass[dplyr::cur_column()]))
  conditions <- list(
    BaEF = genus_info$type == "BaEF",
    FuEF = genus_info$type == "FuEF",
    PpEF = genus_info$type == "PpEF",
    OpEF = genus_info$type == "OpEF",
    TNEF = rep(TRUE, nrow(genus_info))
  )

  cal_fun <- function(data, genus_info, condition) {
    sel_genus <- genus_info %>%
      dplyr::filter(!!condition) %>%
      dplyr::pull('Genus')
    data %>%
      dplyr::select(dplyr::all_of(sel_genus)) %>%
      rowSums()
  }
  NEF_result <- lapply(names(conditions), function(names) {
    cal_fun(EF_df, genus_info, conditions[[names]])
  })
  NEF_result_df <- as.data.frame(NEF_result)
  colnames(NEF_result_df) <- names(conditions)
  biomass_result <- lapply(names(conditions)[1:4], function(names) {
    cal_fun(biomass_df, genus_info, conditions[[names]])
  })
  Energy.flux.node <- as.data.frame(biomass_result)
  colnames(Energy.flux.node) <- gsub("EF", "", names(conditions)[1:4])
  Energy.flux.node <- data.frame(
    Sample.ID = row.names(Energy.flux.node),
    Energy.flux.node
  )
  row.names(Energy.flux.node) <- NULL
  prefer <- lapply(names(conditions)[1:3], function(names) {
    cal_fun(data, genus_info, conditions[[names]])
  })
  prefer <- as.data.frame(prefer)
  colnames(prefer) <- gsub("EF", "", names(conditions)[1:3])
  total_prefer <- rowSums(prefer)
  rel_prefer <- prefer %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x / total_prefer))
  loss <- rel_prefer %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x * (NEF_result_df$OpEF / AE$Op)))
  mapping <- c(Ba = "BaEF", Fu = "FuEF", Pp = "PpEF")
  outflow.R <- loss %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ (.x + NEF_result_df[[mapping[dplyr::cur_column()]]]) / AE[[dplyr::cur_column()]]))
  colnames(outflow.R) <- paste0("R.to.", colnames(outflow.R))
  colnames(loss) <- paste0(colnames(loss), ".to.Op")
  Energy.flux.path <- merge(outflow.R, loss, by = "row.names")
  colnames(Energy.flux.path)[1] <- "Sample.ID"
  U <- Energy.flux.path %>%
    dplyr::rowwise() %>%
    dplyr::mutate(U = mean(dplyr::c_across(-dplyr::all_of("Sample.ID"))) / sd(dplyr::c_across(-dplyr::all_of("Sample.ID")))) %>%
    dplyr::select(dplyr::all_of(c("Sample.ID", "U"))) %>%
    as.data.frame()
  result <- list(
    data = list(
      data = data,
      Abundance = abundance
    ),
    Energy.flux = list(
      Energy.flux = NEF_result_df,
      C.flux.node = Energy.flux.node,
      C.flux.path = Energy.flux.path,
      U = U
    )
  )
  class(result) <- "NEF"
  return(result)
}


#' @rdname NEF
#' @method NEF matrix
#' @exportS3Method Nematode::NEF
NEF.matrix <- function(data, abundance, AE = list(Ba = 0.6, Fu = 0.38, Pp = 0.25, Op = 0.5), ...) {
  if (is.null(rownames(data)) || is.null(colnames(data))) {
    stop("Data must have row names and column names")
  }

  if (!all(row.names(data) %in% row.names(abundance))) {
    stop("Please provide the abundance information for all samples!")
  }

  data[is.na(data)] <- 0

  if (!all(is.numeric(data))) {
    stop("The data contains non-numeric characters!")
  }

  if (!all(sapply(abundance, is.numeric))) {
    stop("The abundance contains non-numeric characters!")
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

  if (any(is.na(genus_info$CP_group))) {
    missing_cp <- genus_info$Query.genus[is.na(genus_info$CP_group)]
    stop(
      length(missing_cp), " genus lack CP values: ", paste(head(missing_cp), collapse = ", "),
      ifelse(length(missing_cp) > 6, "...", "")
    )
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
    "Bacterial feeders" = "BaEF",
    "Fungus feeders" = "FuEF",
    "Plant feeders" = "PpEF",
    "Omnivores" = "OpEF",
    "Predators" = "OpEF"
  )
  genus_info$type <- type_map[genus_info$Feeding_habit]
  genus_info$EF <- (0.1 * (genus_info$mass / (genus_info$CP_group * 12))) + (0.0159 * (genus_info$mass^0.75))

  total_ident <- rowSums(data)
  abundance_value <- abundance[match(row.names(data), row.names(abundance)), colnames(abundance)[1]]
  rel_data <- sweep(data, 1, total_ident, "/") * abundance_value

  mul_map_EF <- setNames(genus_info$EF, genus_info$Genus)
  EF_mat <- sweep(rel_data, 2, mul_map_EF[colnames(rel_data)], "*")

  mul_map_biomass <- setNames(genus_info$mass, genus_info$Genus)
  biomass_mat <- sweep(rel_data, 2, mul_map_biomass[colnames(rel_data)], "*")

  conditions <- list(
    BaEF = genus_info$type == "BaEF",
    FuEF = genus_info$type == "FuEF",
    PpEF = genus_info$type == "PpEF",
    OpEF = genus_info$type == "OpEF",
    TNEF = rep(TRUE, nrow(genus_info))
  )

  cal_fun <- function(data_mat, genus_info, condition) {
    sel_genus <- genus_info[["Genus"]][condition]
    rowSums(data_mat[, sel_genus, drop = FALSE])
  }

  NEF_result <- lapply(names(conditions), function(names) {
    cal_fun(EF_mat, genus_info, conditions[[names]])
  })
  NEF_result_mat <- do.call(cbind, NEF_result)
  colnames(NEF_result_mat) <- names(conditions)

  biomass_result <- lapply(names(conditions)[1:4], function(names) {
    cal_fun(biomass_mat, genus_info, conditions[[names]])
  })
  Energy.flux.node <- do.call(cbind, biomass_result)
  colnames(Energy.flux.node) <- gsub("EF", "", names(conditions)[1:4])
  Energy.flux.node <- data.frame(
    Sample.ID = row.names(Energy.flux.node),
    Energy.flux.node,
    stringsAsFactors = FALSE
  )
  row.names(Energy.flux.node) <- NULL

  prefer <- lapply(names(conditions)[1:3], function(names) {
    cal_fun(data, genus_info, conditions[[names]])
  })
  prefer <- do.call(cbind, prefer)
  colnames(prefer) <- gsub("EF", "", names(conditions)[1:3])
  total_prefer <- rowSums(prefer)
  rel_prefer <- sweep(prefer, 1, total_prefer, "/")

  loss <- sweep(rel_prefer, 1, NEF_result_mat[, "OpEF"], "*") / AE$Op
  mapping <- c(Ba = "BaEF", Fu = "FuEF", Pp = "PpEF")

  outflow.R <- vapply(colnames(loss), function(col) {
    (loss[, col] + NEF_result_mat[, mapping[col]]) / AE[[col]]
  }, numeric(nrow(loss)))
  colnames(outflow.R) <- paste0("R.to.", colnames(loss))
  colnames(loss) <- paste0(colnames(loss), ".to.Op")

  Energy.flux.path <- data.frame(
    Sample.ID = row.names(data),
    outflow.R,
    loss,
    stringsAsFactors = FALSE
  )

  U_values <- apply(Energy.flux.path[, -1], 1, function(x) mean(x) / sd(x))
  U <- data.frame(
    Sample.ID = Energy.flux.path$Sample.ID,
    U = U_values,
    stringsAsFactors = FALSE
  )

  row.names(Energy.flux.node) <- NULL
  row.names(Energy.flux.path) <- NULL
  row.names(U) <- NULL
  result <- list(
    data = list(
      data = data,
      Abundance = abundance
    ),
    Energy.flux = list(
      Energy.flux = as.data.frame(NEF_result_mat),
      C.flux.node = Energy.flux.node,
      C.flux.path = Energy.flux.path,
      U = U
    )
  )
  class(result) <- "NEF"
  return(result)
}

#' @rdname NEF
#' @method NEF default
#' @exportS3Method Nematode::NEF
NEF.default <- function(data, abundance, AE = list(Ba = 0.6, Fu = 0.38, Pp = 0.25, Op = 0.5), ...) {
  # Error for unsupported input types
  stop("Unsupported input type: ", paste(class(data), collapse = "/"), call. = FALSE)
}


# summary.NEF <- function(object) {
#   if (!inherits(object, "NEF")) {
#     stop("The object is not of class 'NEF'.")
#   }
#
#   energy_flux <- object$Energy.flux
#
#   cat("Summary of Nematode.Energy.flux:\n")
#   cat("Energy.flux.node:\n")
#   print(energy_flux$Energy.flux.node)
#
#   cat("\nEnergy.flux.path:\n")
#   print(energy_flux$Energy.flux.path)
#
#   cat("\nU values:\n")
#   print(energy_flux$U)
# }
