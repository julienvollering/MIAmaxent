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


#' Selected derived variables accompanied by selection trails.
#'
#' Selected derived variables and tables showing forward model selection of
#' derived variables for modeling the distribution of semi-natural grasslands in
#' Oestfold County, Norway.
#'
#' @format List with 2 elements: \enumerate{ \item A list of 10 data frames,
#'   with each containing the derived variables selected for a given explanatory
#'   variable. \item A list of forward model selection trails used to select
#'   derived variables. }
#' @source produced from \code{\link{grasslandDVs}} and \code{\link{grasslandPO}} using
#'   \code{\link{selectDVforEV}}
"grasslandDVselect"


#' Selected explanatory variables accompanied by selection trails.
#'
#' Selected explanatory variables and tables showing forward model selection of
#' explanatory variables for modeling the distribution of semi-natural
#' grasslands in Oestfold County, Norway. Each individual explanatory variable
#' is represented by a group of derived variables.
#'
#' @format List with 2 elements: \enumerate{ \item A list of 13 data frames,
#'   where the first 9 represent selected explanatory variables, and the last 4
#'   represent selected interaction terms between these explanatory variables.
#'   \item A trail of forward model selection used to select explanatory
#'   variables and interaction terms. }
#' @source produced from \code{\link{grasslandDVselect}} and
#'   \code{\link{grasslandPO}} using \code{\link{selectEV}}
"grasslandEVselect"