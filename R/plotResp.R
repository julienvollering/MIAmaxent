#' Plots model response across a given explanatory variable.
#'
#' \code{plotResp} plots  the response of a given Maxent model over any of the
#' included explanatory variables (EVs) in that model. For categorical
#' variables, a box plot is returned rather than a curve. \code{plotResp} also
#' returns a data frame containing the plotted data (for customizable graphics).
#'
#' The response curves generated in this function are single-effect response
#' curves, presenting the response of a model containing the explanatory
#' variable of interest only (cf. marginal-effect response curves).
#'
#' @param rv Response variable vector used to train the model. The RV should
#'   represent presence/background data, coded as: 1/NA.
#' @param ev Explanatory data used to train the model. Named list of data
#'   frames, with each data frame containing 1 or more DVs for a given EV. E.g.
#'   output [[1]] of \code{selectEV}.
#' @param jarpath The pathway to the maxent.jar executable jar file.
#'
#' @return In addition to the graphical output, a data frame containing the
#'   plotted data is returned.
#'
#' @export


plotResp <- function(rv, ev, jarpath = NULL) {

  altrMaxent:::.binaryrvcheck(rv)

  if (file.exists(jarpath) == F) {
    stop("The pathway of the maxent.jar file must be correctly specified by the
      jarpath parameter. \n ")
  }


}
