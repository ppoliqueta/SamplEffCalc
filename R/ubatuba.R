#' ubatuba Dataset
#'
#' The ubatuba dataset contains data on polychaete communities collected from sediment samples
#' in shelf environments off the coast of Ubatuba city, São Paulo, Brazil.
#'
#' The dataset includes:
#' Environmental variables (e.g., depth, temperature, sediment characteristics),
#' community parameters (e.g., richness, abundance, diversity indices) and
#' counts of dominant polychaete species.
#'
#' @format A data frame with 54 rows and 28 variables:
#' \describe{
#'   \item{Sample}{Number of the sample}
#'   \item{depth}{local depth in meters}
#'   \item{temp}{temperature in Celsius degrees}
#'   \item{sal}{salinity in PSU}
#'   \item{ox}{oxygen content in mg}
#'   \item{om}{\% of sediment organic matter content}
#'   \item{caco3}{\% of sediment calcium carbonate content}
#'   \item{mud}{\% of mud}
#'   \item{cs}{\% of coarse sand}
#'   \item{ms}{\% of medium sand}
#'   \item{fs}{\% of fine sand}
#'   \item{selec}{selection coefficient in phi}
#'   \item{N}{total abundance}
#'   \item{S}{number of species}
#'   \item{H}{Shannon diversity index}
#'   \item{E}{Pielou's index of evenness}
#'   \item{Hlunl}{abundance of the polychaete *Harmothoe lunulata*}
#'   \item{Ptri}{abundance of the polychaete *Parandalia tricuspis*}
#'   ...
#' }
#'
#' @source Paiva, P.C. 1993. Anelídeos poliquetos da plataforma continental norte do Estado de São Paulo:
#' I. Padrões de densidade e diversidade específica. Boletim do Instituto Oceanográfico, São Paulo, 41(1/2): 69-80.
#'
#' @examples
#' data("ubatuba")
#' head(ubatuba)
"ubatuba"
