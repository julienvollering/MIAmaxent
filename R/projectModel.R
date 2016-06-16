#' Project model to data.
#'
#' \code{projectModel} Calculates the probability ratio output (PRO) of a given
#' model for any points where values of the explanatory variables in the model
#' are known.The transformations performed on the explanatory variables to build
#' the model must be specified.
#'
#' @param newdata Data frame containing data for all the explanatory variables
#'   included in the model.
#' @param transf Pathway to the .Rdata file containing the named parameterized
#'   transformations used in the model. This file is saved as a result of the
#'   \code{deriveVars} function.
#' @param model Pathway to the .lambdas file of the model in question. This file
#'   is saved as a result of the \code{selectEV} function.
#' @param clamping logical. Do clamping. Default is \code{FALSE}.
#'
#' @return 1) a data frame with the model output and the corresponding
#'   explanatory data. 2) a data frame showing the range of the data compared to
#'   the training data, on a 0-1 scale.
#'
#' @export


projectModel <- function(newdata, transf, model, clamping = FALSE) {

  lambdas <- read.csv(model, header = FALSE)
  dvnames <- as.character(lambdas[1:(nrow(lambdas)-4), 1])
  m <- length(dvnames)
  dvnamesni <- dvnames[-grep(":", dvnames)]
  dvnamesi <- dvnames[grep(":", dvnames)]

  transf <- load(transffile)

  evnames <- unique(unname(sapply(dvnamesni, function(x) {
    colnames(newdata)[startsWith(x, colnames(newdata))]
    })))

  Ranges <- lapply(evnames, function(x) {
    evdata <- newdata[, x]
    xnull <- environment(Storage[startsWith(names(Storage), x)][[1]])$xnull
    if (class(xnull) == "numeric" || class(xnull) == "integer") {
      L <- (evdata - range(xnull)[1])/diff(range(xnull))
      return(range(L))
    }
    if (class(xnull) == "factor" || class(xnull) == "character") {
      if (all(evdata %in% xnull)) {return("inside")} else {return("outside")}
    }
  })
  names(Ranges) <- evnames

  ranges <- list()
  dvdatani <- list()
  lapply(dvnamesni, function(x) {
    evname <- colnames(newdata)[startsWith(x, colnames(newdata))]
    evdata <- newdata[, evname]
    transffunction <- Storage[startsWith(names(Storage), x)][[1]]
    xnull <- environment(transffunction)$xnull
    if (class(xnull) == "numeric" || class(xnull) == "integer") {
      L <- (evdata - range(xnull)[1])/diff(range(xnull))
      ranges[[evname]] <- range(L)
    }
    if (class(xnull) == "factor" || class(xnull) == "character") {
      if (all(evdata %in% xnull)) {
        ranges[[evname]] <- "inside"
      } else {
        ranges[[evname]] <- "outside"
      }
    }
    dvdatani[[x]] <- transffunction(evdata)
    return(list(ranges, dvdatani))
  })
}
