#' calculates exponentially weighted moving average
#'
#' @param x numeric. Vector across which the moving average is to be applied.
#' @param n integer. Width of the moving average window. Should be odd,
#'   otherwise the window will be uncentered.
#'
#' \code{DESCRIPTION Imports}: stats
#'

ewma <- function(x, n) {
  if (missing(n)) {
    stop("Specify the width of the moving average window (n)", call. = F)
  }
  if (n < 3) {
    stop("Width of window should be at least 3", call. = F)
  }

  if (n %% 2 != 0) {
    expwindow <- dexp(c(((n-1)/2):0,1:((n-1)/2)))
  } else {
    expwindow <- dexp(c((n/2):0,1:((n-2)/2)))
  }
  weights <- expwindow/sum(expwindow)
  as.numeric(stats::filter(x, weights, sides=2))
}