#' Project model to data.
#'
#' \code{projectModel} Calculates the probability ratio output (PRO) of a given
#' model for any points where values of the explanatory variables in the model
#' are known.The transformations performed on the explanatory variables to build
#' the model must be specified.
#'
#' Missing data (NA) for a continuous variable will result in NA output for that
#' point. Missing data for a categorical variable is treated as belonging to
#' none of the categories.
#'
#' @param data Data frame of all the explanatory variables (EVs) included in the
#'   model, with column names matching EV names. See \code{\link{readData}}.
#' @param transformation Full pathway of the 'transformations.Rdata' file
#'   containing the transformations used to build the model. This file is saved
#'   as a result of the \code{\link{deriveVars}} function. Equivalently, the
#'   second item in the list returned by \code{\link{deriveVars}} can be used
#'   directly.
#' @param model Full pathway of the '.lambdas' file of the model in question.
#'   This file is saved as a result of \code{\link{selectEV}}.
#' @param clamping Logical. Do clamping \emph{sensu} Phillips et al. (2006).
#'   Default is \code{FALSE}.
#'
#' @return List of 2: \enumerate{ \item A data frame with the model output in
#'   column 1 and the corresponding explanatory data in subsequent columns.
#'   \item A data frame showing the range of \code{data} compared to the
#'   training data, on a 0-1 scale.}
#'
#' @references Phillips, S.J., Anderson, R.P. & Schapire, R.E. (2006) Maximum
#'   entropy modeling of species geographic distributions. Ecological Modelling,
#'   190, 231-259.
#'
#' @export


projectModel <- function(data, transformation, model, clamping = FALSE) {

  lambdas <- utils::read.csv(model, header = FALSE)
  dvnames <- as.character(lambdas[1:(nrow(lambdas)-4), 1])
  dvnamesni <- dvnames[-grep(":", dvnames)]
  dvnamesi <- dvnames[grep(":", dvnames)]

  check <- lapply(dvnamesni, function(x) { startsWith(x, colnames(data)) })
  if (any(sapply(check, function(x) { all(x==FALSE) }))) {
    stop("All EVs in the model must be represented in data",
      call. = FALSE)
  }

  alltransf <- .load.transf(transformation)

  evnames <- unique(unname(sapply(dvnamesni, function(x) {
    colnames(data)[startsWith(x, colnames(data))]
  })))
  evnames <- sapply(evnames, .best.match, b = dvnamesni)

  check2 <- lapply(dvnamesni, function(x) { startsWith(names(alltransf), x) })
  if (any(sapply(check2, sum) != 1)) {
    stop("All DVs in the model must be represented in transformation",
      call. = FALSE)
  }

  Ranges <- lapply(evnames, function(x) {
    evdata <- data[, x]
    xnull <- environment(alltransf[startsWith(names(alltransf), x)][[1]])$xnull
    if (class(xnull) %in% c("numeric", "integer")) {
      L <- (evdata - range(xnull)[1]) / diff(range(xnull))
      return(range(L))
    }
    if (class(xnull) %in% c("factor", "character")) {
      if (all(evdata %in% xnull)) {return("inside")} else {return("outside")}
    }
  })
  names(Ranges) <- evnames

  dvdatani <- lapply(dvnamesni, function(x) {
    evname <- evnames[startsWith(x, evnames)]
    evdata <- data[, evname]
    transffunction <- alltransf[startsWith(names(alltransf), x)][[1]]
    y <- transffunction(evdata)
    if (clamping == TRUE) {
      y[y > 1] <- 1
      y[y < 0] <- 0
    }
    return(y)
  })
  names(dvdatani) <- dvnamesni

  prodlist <- strsplit(dvnamesi, ":")
  dvdatai <- lapply(prodlist, function(x) {
    dvdatani[[x[1]]] * dvdatani[[x[2]]]
  })
  names(dvdatai) <- dvnamesi

  dvdata <- data.frame(c(dvdatani, dvdatai), check.names = FALSE)
  modelfunction <- modelfromlambdas(model)
  PRO <- modelfunction(dvdata)[, 1]
  Output <- cbind(PRO, data)

  return(list(output = Output, ranges = Ranges))
}
