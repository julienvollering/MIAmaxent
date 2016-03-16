#' produces dvs from a given ev by different transformations.
#'
#' In MIAT, "D" transformation is only performed if the optimum occurs in the
#' middle 80\% of the EV range. \code{dvfromev} does not currently specify this
#' condition.
#'
#' @param df Dataframe with 2 columns: response variable and explanatory
#'   variable (in that order). Column names are used as identifiers. The
#'   response variable represents presence or background, coded as: 1/NA. The
#'   explanatory variable may be continuous or categorical.
#' @param writedir Directory to which Maxent runs of spline transformations are
#'   written
#' @param transformtype Set of transformation types to be used.
#' @param allsplines Logical. Keep all spline transformations.
#'
#' @return Dataframe with one column for each DV.


dvfromev <- function(df, writedir, transformtype, allsplines) {

  rv <- df[,1]
  ev <- df[,2]
  evname <- colnames(df)[2]
  evdv <- data.frame(df[,2])
  colnames(evdv) <- evname

  if (class(ev) == "numeric" || class(ev) == "integer") {

    if ("L" %in% transformtype) {
      L <- data.frame((ev - range(ev)[1])/diff(range(ev)))
      colnames(L) <- paste(evname, "_L", sep = "")
      evdv <- cbind(evdv, L)
    }

    if ("M" %in% transformtype) {
      L <- (ev - range(ev)[1])/diff(range(ev))
      ZSk <- altrMaxent::.scalex(L, altrMaxent::.minskew(L)$c)
      M <- data.frame((ZSk - range(ZSk)[1])/diff(range(ZSk)))
      colnames(M) <- paste(evname, "_M", sep = "")
      evdv <- cbind(evdv, M)
    }

    if ("D" %in% transformtype) {
      L <- (ev - range(ev)[1])/diff(range(ev))
      opt <- altrMaxent::.fopoptimum(data.frame(rv, L))
      devexp = c(0.5, 1, 2)
      D <- altrMaxent::.devtransf(L, opt, devexp)
      colnames(D) <- paste(evname, "_D", devexp, sep = "")
      evdv <- cbind(evdv, D)
    }
  }

  if (class(ev) == "factor" || class(ev) == "character") {

  }
  return(evdv)
}
