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
#'
#' @return 1) a data frame with the model output and the corresponding
#'   explanatory data. 2) a data frame showing the range of the data compared to
#'   the training data, on a 0-1 scale.
#'
#' @export


projectModel <- function(data, EV, intervals = 20, smoothwindow = 3,
                    EVranging = FALSE) {

}
