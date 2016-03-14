#' produces dvs from a given ev by different transformations.
#'
#' @param rv Response variable (vector) representing presence or background,
#'   coded as: 1/NA.
#' @param ev Explanatory variable. May be continuous or categorical.
#' @param writedir Directory to which Maxent runs of spline transformations are
#'   written
#' @param transformtype Set of transformation types to be used.
#' @param allsplines Logical. Keep all spline transformations.
#'
#' @return Dataframe with one column for each DV.


dvfromev <- function(rv, ev, writedir, transformtype, allsplines) {
  if (class(ev) == "numeric" || class(ev) == "integer") {

  }
  if (class(ev) == "factor" || class(ev) == "character") {

  }
}
