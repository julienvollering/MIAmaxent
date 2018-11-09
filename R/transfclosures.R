#' A lexical closure for defining an "L" transformation
#'
#' @param xnull training data that will be used to parameterize the
#'   transformation
#' @return Function that transforms x into y by the particular L transfomation
#'   defined by the properties of xnull
#' @keywords internal
#' @noRd

.transfL <- function(xnull) {
  function(x) {
    y <- (x - range(xnull)[1])/diff(range(xnull))
    return(y)
  }
}

#' A lexical closure for defining an "M" transformation
#'
#' @param xnull training data that will be used to parameterize the
#'   transformation
#' @return Function that transforms x into y by the particular M transfomation
#'   defined by the properties of xnull
#' @keywords internal
#' @noRd

.transfM <- function(xnull) {
  Lnull <- (xnull - range(xnull)[1])/diff(range(xnull))
  c <- .minskew(Lnull)$c
  ZSknull <- .scalex(Lnull, Lnull, c)

  function(x) {
    L <- (x - range(xnull)[1])/diff(range(xnull))
    ZSk <- .scalex(Lnull, L, c)
    y <- (ZSk - range(ZSknull)[1])/diff(range(ZSknull))
    return(y)
  }
}

#' A lexical closure for defining an "D" transformation
#'
#' @param xnull training data that will be used to parameterize the
#'   transformation
#' @param rv response variable that will be used to parameterize the
#'   transformation
#' @param devexp exponent determining the steepness of the deviation
#'   transformation
#' @return Function that transforms x into y by the particular D transfomation
#'   defined by the properties of xnull and rv, and specified by devexp.
#' @keywords internal
#' @noRd

.transfD <- function(xnull, rv, devexp) {
  optnull <- .fopoptimum(data.frame(rv, xnull))
  uynull <- (abs(xnull - optnull)) ^ devexp

  function(x) {
    uy <- (abs(x - optnull)) ^ devexp
    y <- (uy - range(uynull)[1])/diff(range(uynull))
    return(y)
  }
}

#' A lexical closure for defining a spline ("HF", "HR", "T") transformation
#'
#' @param xnull training data that will be used to parameterize the
#'   transformation
#' @param k Break point between 0 and 1.
#' @param type Specifies type of spline transformation: "HF", "HR", or "T".
#' @return list of functions that transform x into y by the particular spline
#'   transformations defined by the properties of xnull, and specified by k and
#'   type.
#' @keywords internal
#' @noRd

.transfSpline <- function(xnull, k, type) {
  force(k)
  function(x) {
    L <- (x - range(xnull)[1])/diff(range(xnull))
    if (type == "HF") {
      y <- ((L - k) * (1 * (L > k))) / (1 - k)
    }
    if (type == "HR") {
      y <- ((L - k) * (-1 * (L < k))) / k
    }
    if (type == "T") {
      y <- 1 * (L >= k)
    }
    return(y)
  }
}

#' A lexical closure for defining x number of "B" transformations
#'
#' @param xnull categorical training data that will be used to parameterize the
#'   transformation
#' @param lvl string specifying which level of the variable to be made binary
#' @return One function for each level of xnull, which transforms categorical
#'   variable x into binary variable y.
#' @keywords internal
#' @noRd

.transfB <- function(xnull, lvl) {
  force(lvl)
  function(x) {
    y <- unname(sapply(x, function(x) {
      if (is.na(x)) {0} else {
      if (x == lvl) {1} else {0}}}))
    return(y)
  }
}