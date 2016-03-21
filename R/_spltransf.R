#' calculates derived variables by hinge transformation
#'
#' Note that this function automatically rescales to [0,1]
#'
#' @param ev Vector of values to be transformed
#' @param knots Number of break points (1 tranformation per break)
#' @param type Specifies type of spline transformation: "HF", "HR", or "T"
#'

.spltransf <- function(ev, knots, type) {
  stopifnot(class(ev) == "numeric" || class(ev) == "integer",
    type %in% c("HF", "HR", "T"))

  spl <- matrix(nrow = length(ev), ncol = knots)
  for (i in 1:knots) {
    k <-(2 * i - 1) / (2 * knots)
    if (type == "HF") {
      spl[,i] <- ((ev - k) * (1 * (ev > k))) / (1 - k)
    }
    if (type == "HR") {
      spl[,i] <- ((ev - k) * (-1 * (ev < k))) / k
    }
    if (type == "T") {
      spl[,i] <- 1 * (ev >= k)
    }
  }
  Spl <- as.data.frame(spl)
  return(Spl)
}