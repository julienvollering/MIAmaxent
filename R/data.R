#' Simulated occurrence and environmental data.
#'
#' Simulated data set for distribution modeling consisting of occurrence and
#' environmental data. The simulated study area consists of 40 grid cells, with
#' 8 row and 5 columns, in which 10 presences occur.
#'
#' @format A data frame with 40 rows and 5 variables: \describe{
#'   \item{RV}{response variable, occurrence either presence or uninformed
#'   background} \item{EV11}{explanatory variable: northing}
#'   \item{EV12}{explanatory variable: easting} \item{EV13}{explanatory
#'   variable: modified random uniform} \item{EV14}{explanatory variable: random
#'   uniform} }
#' @source Halvorsen, R. (2013) A strict maximum likelihood explanation of
#'   MaxEnt, and some implications for distribution modelling. Sommerfeltia, 36,
#'   1-132.
"sp1po"


#' Derived variables and transformation functions.
#'
#' Derived variables and transformation functions for distribution modeling of a
#' simulated data set.
#'
#' @format List with 2 elements: \enumerate{ \item A list of 4 data frames, with
#'   each containing the derived variables produced for a given explanatory
#'   variable. \item A list of all the transformation functions used to produce
#'   the derived variables. }
#' @source produced from \code{\link{sp1po}} using \code{\link{deriveVars}}
"DVs"


#' Selected derived variables accompanied by selection trails.
#'
#' Selected derived variables and tables showing forward model selection of
#' derived variables for distribution modeling of a simulated data set.
#'
#' @format List with 2 elements: \enumerate{ \item A list of 3 data frames, with
#'   each containing the derived variables selected for a given explanatory
#'   variable. \item A list of forward model selection trails used to select
#'   derived variables. }
#' @source produced from \code{\link{DVs}} and \code{\link{sp1po}} using
#'   \code{\link{selectDVforEV}}
"seldvs"


#' Selected explanatory variables accompanied by selection trails.
#'
#' Selected explanatory variables and tables showing forward model selection of
#' explanatory variables for distribution modeling of a simulated data set. Each
#' individual explanatory variable is represented by a group of derived
#' variables.
#'
#' @format List with 2 elements: \enumerate{ \item A list of 5 data frames,
#'   where the first 3 represent selected explanatory variables, and the last 2
#'   represent selected interaction terms between these explanatory variables.
#'   \item A trail of forward model selection used to select explanatory
#'   variables and interaction terms. }
#' @source produced from \code{\link{seldvs}} and
#'   \code{\link{sp1po}} using \code{\link{selectEV}}
"selevs"