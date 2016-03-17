#' calculates derived variables by hinge transformation
#'
#' Note that this function automatically rescales to [0,1]
#'
#' @param ev Vector of values to be transformed
#' @param knots Number of break points (1 tranformation per break)
#' @param forward Logical value specifying forward or reverse hinge
#'   transformation
#'

.hintransf <- function(ev, knots, forward = TRUE) {
  stopifnot(class(ev) == "numeric" || class(ev) == "integer",
    length(knots) == 1)

  hin <- matrix(nrow = length(ev), ncol = knots)
  for (i in 1:knots) {
    k <-(2 * i - 1) / (2 * knots)
    if (forward == T) {
      hin[,i] <- ((ev - k) * (1 * (ev > k))) / (1 - k)
    } else {
      hin[,i] <- ((ev - k) * (-1 * (ev < k))) / k
    }
  }
  Hin <- as.data.frame(hin)
  return(Hin)
}