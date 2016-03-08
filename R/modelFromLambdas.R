#' Maxent model from .lambdas file.
#'
#' \code{modelFromLambdas} returns an R function for a given Maxent model,
#' using the .lambdas file produced by Maxent.jar to parameterize the model.
#' The returned function gives predictions of the model in "raw output" format
#' from values of explanatory variables.
#'
#' The \code{modelFromLambdas} function returns a Maxent model, in the form of
#' a function, which calculates Maxent predictions for given values of
#' explanatory variables.
#'
#' Input to the function returned by \code{modelFromLambdas} are values of
#' explanatory variables in the model. The format of this input must be an
#' array (matrix or data frame) with m columns, and column names must match
#' variable names in the .lambdas file used to reproduce the model.
#'
#' @param file pathway to the .lambdas file of a given Maxent model
#' @param name Character(). Specifies the name of the function to be returned.
#' Defaults to the basename of the .lambdas file.
#'
#' @return returns an R function (object).
#'
#' @examples
#' \dontrun {
#' modelFromLambdas("D:\\tutorial-data\\bradypus_variegatus.lambdas", name="mymodel")
#' mymodel(projection.df)
#' }

modelFromLambdas <- function(file, name = NULL) {

  if (is.null(name)) {
    name <- sub(".lambdas", "", basename(file))
  }

  lambdas <- read.csv(file, header = FALSE)
  dvrows <- lambdas[1:(nrow(lambdas)-4),]
  m <- nrow(dvrows)
  thetas <- dvrows[,2]
  xmins <- dvrows[,3]
  xmaxs <- dvrows[,4]
  linPredNorm <- lambdas[nrow(lambdas)-3,2]
  densNorm <- lambdas[nrow(lambdas)-2,2]
  model <- function(X) {
    matchorder <- match(as.character(dvrows[,1]), colnames(X))

    if (ncol(X) != m) {
      stop("Input must have as many columns as there are variables in the model", call. = F)
    }

    if (length((matchorder[!is.na(matchorder)])) < m) {
      stop("Input column names must match the names of the variables in the model", call. = F)
    }

    orderedX <- as.matrix(X[,matchorder])
    thetaX <- matrix(nrow = nrow(orderedX), ncol = m)
    for (j in 1:nrow(orderedX)) {
      for (i in 1:m) {
        thetaX[j,i] <- thetas[i]*((orderedX[j,i]-xmins[i])/(xmaxs[i]-xmins[i]))
      }
    }
    rawoutput <- apply(thetaX, 1, function(x) ((exp(sum(x)-linPredNorm))/densNorm))
    projection <- cbind(orderedX, rawoutput)
    return(projection)
  }
  assign(name, model, pos = ".GlobalEnv")
}
