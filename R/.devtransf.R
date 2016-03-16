#' calculates derived variables by deviation transformation
#'
#' Note that this function calculates devation AND rescales to [0,1]
#'
#' @param ev Vector of values to be transformed
#' @param optimum Center value from which deviation is calculated
#' @param devexp Vector of exponent values to scale the steepness of the
#'   deviation transformation
#'

devtransf <- function(ev, optimum, devexp = c(0.5, 1, 2)) {
  stopifnot(class(ev) == "numeric" || class(ev) == "integer",
    class(optimum) == "numeric" || class(optimum) == "integer",
    length(optimum) == 1)

  dev <- matrix(nrow = length(ev), ncol = length(devexp))
  for (i in 1:length(devexp)) {
    unscaled <- (abs(ev-optimum))^devexp[i]
    dev[,i] <- (unscaled - range(unscaled)[1])/diff(range(unscaled))
  }
  Dev <- as.data.frame(dev)
  return(Dev)
}