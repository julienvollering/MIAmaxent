#' Best prefix match.
#'
#' Return only the string in 'a' which best matches any of the strings in 'b'
#'
#' @param a Vector of character strings to be matched
#' @param b vector of character strings to match
#'
#' @return single character string


best.match <- function(a, b) {
  score <- sapply(a, function(a) {
    reg <- regexpr(a, b)
    max(attr(reg, "match.length"))
  })
  a[which.max(score)]
}
