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
#' @param data Data frame containing data for all the explanatory variables
#'   (EVs) included in the model, with column names matching EV names.
#' @param transformation Pathway to the .Rdata file containing the named
#'   parameterized transformations used in the model. This file is saved as a
#'   result of the \code{\link{deriveVars}} function. Alternatively, a list
#'   object containing all of the named transformations (e.g. the second item in
#'   the list returned by \code{\link{deriveVars}}).
#' @param model Pathway to the .lambdas file of the model in question. This file
#'   is saved as a result of the \code{selectEV} function.
#' @param clamping logical. Do clamping. Default is \code{FALSE}.
#'
#' @return 1) a data frame with the model output and the corresponding
#'   explanatory data. 2) a data frame showing the range of the data compared to
#'   the training data, on a 0-1 scale.
#'
#' @export


projectModel <- function(data, transformation, model, clamping = FALSE) {

  lambdas <- read.csv(model, header = FALSE)
  dvnames <- as.character(lambdas[1:(nrow(lambdas)-4), 1])
  dvnamesni <- dvnames[-grep(":", dvnames)]
  dvnamesi <- dvnames[grep(":", dvnames)]

  check <- lapply(dvnamesni, function(x) { startsWith(x, colnames(data)) })
  if (any(sapply(check, function(x) { all(x==FALSE) }))) {
    stop("All EVs in the .lambdas file must be represented in data",
      call. = FALSE)
  }

  if (class(transformation) == "character") {
    alltransf <- get(load(transformation))
  } else {
    alltransf <- transformation
  }
  if (!all(sapply(alltransf, class) == "function")) {
    stop("transformation argument should contain functions only", call. = FALSE)
  }

    evnames <- unique(unname(sapply(dvnamesni, function(x) {
    colnames(data)[startsWith(x, colnames(data))]
    })))
    evnames <- sapply(evnames, .best.match, b = dvnamesni)

  Ranges <- lapply(evnames, function(x) {
    evdata <- data[, x]
    xnull <- environment(alltransf[startsWith(names(alltransf), x)][[1]])$xnull
    if (class(xnull) == "numeric" || class(xnull) == "integer") {
      L <- (evdata - range(xnull)[1])/diff(range(xnull))
      return(range(L))
    }
    if (class(xnull) == "factor" || class(xnull) == "character") {
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
