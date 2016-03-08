#' Maxent model from .lambdas file.
#'
#' \code{modelFromLambdas} returns an R function for a given Maxent model,
#' using the .lambdas file produced by Maxent.jar to parameterize the model.
#' The returned function gives predictions of the model in "raw output" format
#' from values of explanatory variables.
#'
#' @param file pathway to the .lambdas file of a given Maxent model
#'
#' @return returns an R function (object).
#
# This is a function named 'modelFromLambdas' which reproduces a model
# produced by maxent.jar in R, from the model's lambda file.

# Details: the modelFromLambdas function returns a maxent model in the form of a
# function which calculates maxent predictions for the given values of the
# explanatory variables.

# Input to the resulting function must be an array (matrix or df) with m
# columns, where column names match variable names in lambda file.

modelFromLambdas <- function(file) { # file is the pathway to the lambdas file of the desired Maxent model
  lambdas <- read.csv(file, header = FALSE)
  dvrows <- lambdas[1:(nrow(lambdas)-4),]
  m <- nrow(dvrows)
  thetas <- dvrows[,2]
  xmins <- dvrows[,3]
  xmaxs <- dvrows[,4]
  linPredNorm <- lambdas[nrow(lambdas)-3,2]
  densNorm <- lambdas[nrow(lambdas)-2,2]
  function(X) {
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
}
