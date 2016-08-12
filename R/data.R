#' Occurrence and environmental data.
#'
#' Occurrence and environmental data for modeling the distribution of
#' semi-natural grasslands in Oestfold County, Norway. A dataset containing the
#' presence-only occurrence along with 13 environmental variables.
#'
#' @format A data frame with 17497 rows and 14 variables: \describe{
#'   \item{RV}{response variable, occurrence either presence or uninformed
#'   background} \item{pca1}{explanatory variable 1}
#'   \item{pr_bygall}{explanatory variable 2} \item{pr_tilany}{explanatory
#'   variable 3} \item{teraspif}{explanatory variable 4}
#'   \item{terdem}{explanatory variable 5} \item{terslpdg}{explanatory variable
#'   6} \item{tersolrade}{explanatory variable 7} \item{tertpi09}{explanatory
#'   variable 8} \item{geoberg}{explanatory variable 9}
#'   \item{geolmja1}{explanatory variable 10} \item{lcucor1}{explanatory
#'   variable 11} \item{lcutil_t4}{explanatory variable 12}
#'   \item{terslpps15}{explanatory variable 13} }
#' @source Sabrina Mazzoni and Rune Halvorsen
"grasslandPO"


#' Derived variables and transformation functions.
#'
#' Derived variables and transformation functions for modeling the distribution
#' of semi-natural grasslands in Oestfold County, Norway.
#'
#' @format List with 2 elements: \enumerate{ \item A list of 13 data frames,
#'   with each containing the derived variables produced for a given explanatory
#'   variable. \item A list of all the transformation functions used to produce
#'   the derived variables. }
#' @source produced from \code{\link{grasslandPO}} using
#'   \code{\link{deriveVars}}
"grasslandDVs"