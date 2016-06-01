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
#' @param evdata Explanatory data used to train the model. Named list of data
#'   frames, with each data frame containing 1 or more DVs for a given EV. E.g.
#'   output [[1]] of \code{selectEV}.
#' @param ev Name or list index of the explanatory variable in \code{evdata} for
#'   which the response curve is to be generated.
#' @param writedir The directory to which Maxent files will be written during
#'   subset selection of DVs. Defaults to the working directory.
#' @param jarpath The pathway to the maxent.jar executable jar file. If
#'   unspecified, the function looks for the file in the writedir.
#'
#' @return In addition to the graphical output, a data frame containing the
#'   plotted data is returned.
#'
#' @export


plotResp <- function(rv, evdata, ev, writedir = NULL, jarpath = NULL) {

  altrMaxent:::.binaryrvcheck(rv)

  if (is.null(writedir)) {
    writedir <- getwd()
  }

  if (is.null(jarpath)) {
    jarpath <- paste(writedir, "\\maxent.jar", sep="")
  }

  if (file.exists(jarpath) == F) {
    stop("maxent.jar file must be present in writedir, or its pathway must be
      specified by the jarpath argument. \n ")
  }

  dir <- paste(writedir, "\\plotResp\\", sep="")
  if (file.exists(dir)) {
    stop("The specified writedir already contains a selection of EVs.
      Please specify a different writedir. \n ")
  } else {
    dir.create(dir)
  }

}
