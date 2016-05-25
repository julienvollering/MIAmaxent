#' A lexical closure for defining and storing "L" transformations
#'
#' @param xnull training data that will be used to parameterize the
#'   transformation
#'
#' @return Function that transforms x into y by the particular L transfomation
#'   used on the training data

.Ltransf <- function(xnull) {
  function(x) {
    y <- (x - range(xnull)[1])/diff(range(xnull))
    return(y)
  }
}