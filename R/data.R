#' Nematode Taxonomic and Functional Traits
#'
#' A dataset containing taxonomic classification and functional traits of nematode genera,
#' including feeding habits and ecological group (c-p value).
#'
#' @format A data frame with 2524 rows and 4 variables:
#' \describe{
#'   \item{Genus}{Nematode genus name (character), e.g. "Parascaris", "Heterakis"}
#'   \item{Family}{Taxonomic family name (character), e.g. "Ascarididae", "Heterakidae"}
#'   \item{Feeding_habit}{Feeding behavior category (character), e.g. "Bacterial feeders", "Omnivores"}
#'   \item{CP_group}{Colonizer-Persister group (numeric, 1-5)}
#' }
#'
#' @details
#' This dataset is particularly useful for:
#' \itemize{
#'   \item Ecological studies of soil nematode communities
#'   \item Trophic network analysis
#'   \item Calculating nematode maturity indices (e.g. MI, PPI)
#' }
#'
#' @source
#' Nemaplex.UCDavis.edu; Revision Date: 03/18/2026; Accessed 03/18/2026
#' \itemize{
#'   \item Website: \url{http://nemaplex.ucdavis.edu/}
#' }
#'
#' @examples
#' # Load the data
#' data(nematode.info)
#'
#' # Count nematodes by feeding habit
#' table(nematode.info$Feeding_habit)
#'
#' # Find all genera in Ascarididae family
#' subset(nematode.info, Family == "Ascarididae")
"nematode.info"


#' Nematode Genus and Family Average Body Mass
#'
#' A dataset containing the average dry body mass (in micrograms) of nematode genera and families,
#' compiled from morphological measurements and allometric scaling. Essential for metabolic rate calculations
#' and size-spectrum analyses in soil ecology.
#'
#' @format A data frame with 1094 rows and 4 variables:
#' \describe{
#'   \item{Genus}{Nematode genus name (character), taxonomically validated against Nemaplex database}
#'   \item{Family}{Corresponding taxonomic family (character)}
#'   \item{Genus.Average.Mass}{Mean dry mass per genus (numeric, \eqn{\mu}g)}
#'   \item{Family.Average.Mass}{Mean dry mass per family (numeric, \eqn{\mu}g)}
#'   }
#'
#' @source
#' Nemaplex.UCDavis.edu; Revision Date: 02/02/2026; Accessed 03/18/2026
#' \itemize{
#'   \item Website: \url{http://nemaplex.ucdavis.edu/}
#' }
#'
#' @seealso
#' Use \code{\link{nematode.info}} for complementary trait data.
#'
#' @examples
#' # Load data
#' data(nematode.ave.mass)
#'
#' # Find mass range within a family (e.g. Rhabditidae)
#' rhabditidae <- subset(nematode.ave.mass, Family == "Rhabditidae")
#' range(rhabditidae$Genus.Average.Mass, na.rm = TRUE)
#'
#' # Convert to biomass (example: 100 individuals of Acanthopharynx)
#' 100 * subset(nematode.ave.mass, Genus == "Acanthopharynx")$Genus.Average.Mass
"nematode.ave.mass"
