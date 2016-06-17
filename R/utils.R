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
