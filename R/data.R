#' Occurrence and environmental toy data.
#'
#' A small, synthetic data set for distribution modeling, consisting of
#' occurrence and environmental data, from Halvorsen (2013). The study area
#' consists of 40 grid cells, with 8 row and 5 columns, in which 10 presences
#' occur.
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
"toydata_sp1po"


#' Derived variables and transformation functions, from toy data.
#'
#' Derived variables and transformation functions for distribution modeling of a
#' small, synthetic data set used in Halvorsen (2013).
#'
#' @format List with 2 elements: \enumerate{ \item A list of 5, with the response
#'   variable followed by data frames each containing the derived variables
#'   produced for a given explanatory variable. \item A list of the response
#'   variable and all the transformation functions used to produce the derived
#'   variables. }
#' @source Produced from \code{\link{toydata_sp1po}} using
#'   \code{\link{deriveVars}}.
#' @references Halvorsen, R. (2013) A strict maximum likelihood explanation of
#'   MaxEnt, and some implications for distribution modelling. Sommerfeltia, 36,
#'   1-132.
"toydata_dvs"


#' Selected derived variables accompanied by selection trails, from toy data.
#'
#' Selected derived variables and tables showing forward model selection of
#' derived variables for distribution modeling of a small, synthetic data set
#' used in Halvorsen (2013).
#'
#' @format List with 2 elements: \enumerate{ \item A list of 3, with the
#'   response variable followed by data frames each containing the derived
#'   variables selected for a given explanatory variable. \item A list of the
#'   response variable and forward model selection trails used to select derived
#'   variables. }
#' @source Produced from \code{\link{toydata_dvs}} using
#'   \code{\link{selectDVforEV}}.
#' @references Halvorsen, R. (2013) A strict maximum likelihood explanation of
#'   MaxEnt, and some implications for distribution modelling. Sommerfeltia, 36,
#'   1-132.
"toydata_seldvs"


#' Selected explanatory variables accompanied by selection trails, from toy
#' data.
#'
#' Selected explanatory variables and tables showing forward model selection of
#' explanatory variables for distribution modeling of a small, synthetic data
#' set used in Halvorsen (2013). Each individual explanatory variable is
#' represented by a group of derived variables.
#'
#' @format List with 3 elements: \enumerate{ \item A list of 3, with the
#'   response variable followed by data frames, represent selected explanatory
#'   variables. \item A trail of forward model selection used to select
#'   explanatory variables and interaction terms. \item The selected model }
#' @source Produced from \code{\link{toydata_seldvs}} using
#'   \code{\link{selectEV}}.
#' @references Halvorsen, R. (2013) A strict maximum likelihood explanation of
#'   MaxEnt, and some implications for distribution modelling. Sommerfeltia, 36,
#'   1-132.
"toydata_selevs"