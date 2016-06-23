#' Best prefix match.
#'
#' Return only the string in 'a' which best matches any of the strings in 'b'
#'
#' @param a Vector of character strings to be matched
#' @param b vector of character strings to match
#'
#' @return single character string

.best.match <- function(a, b) {
  score <- sapply(a, function(a) {
    reg <- regexpr(a, b)
    max(attr(reg, "match.length"))
  })
  a[which.max(score)]
}



#' checks the validity of RV values
#'
#' Presence-only data should be coded as: 1/NA (preferred) or 1/0 (danger of
#' misinterpretation as presence/absence data)
#'
#' @param rv Vector of response variable values

.binaryrvcheck <- function(rv) {
  if (length(levels(as.factor(rv))) > 2) {
    stop("The response variable must contain 2 levels only: presence (1)
      and background (NA/0)", call. = FALSE)
  }
  if (anyNA(rv) && length(levels(as.factor(rv))) > 1) {
    stop("The response variable must contain 2 levels only: presence (1)
      and background (NA/0)", call. = FALSE)
  }
  if (class(rv) != "numeric" && class(rv) != "integer") {
    stop("The response variable must be numeric or integer class: presence (1)
      and background (NA/0)", call. = FALSE)
  }
}



#' calculates exponentially weighted moving average
#'
#' @param x numeric. Vector across which the moving average is to be applied.
#' @param n integer. Width of the moving average window. Should be odd,
#'   otherwise the window will be uncentered.
#'
#' \code{DESCRIPTION Imports}: stats
#' @return vector of moving average values

.ewma <- function(x, n) {
  if (missing(n)) {
    stop("Specify the width of the moving average window (n)", call. = FALSE)
  }
  if (n < 3) {
    stop("Width of window should be at least 3", call. = FALSE)
  }

  if (n %% 2 != 0) {
    expwindow <- dexp(c(((n-1)/2):0,1:((n-1)/2)))
  } else {
    expwindow <- dexp(c((n/2):0,1:((n-2)/2)))
  }
  weights <- expwindow/sum(expwindow)
  as.numeric(stats::filter(x, weights, sides=2))
}


#' Make regular intervals
#'
#' @param a Numeric vector
#' @param b number of intervals
#'
#' @return factor variable with 1 level for each interval

.reg.interval <- function(a, b) {
  intwidth <- (max(a) - min(a)) / b
  cutpts <- seq(min(a), max(a), by = intwidth)
  Hmisc::cut2(a, cuts = cutpts, oneval = FALSE)
}
