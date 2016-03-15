#' skewness transformation using constant c
#'
#' \code{DESCRIPTION Imports}: e1071
#'
#' @param x Vector of data.
#' @param c Constant


scalex <- function(x, c) {
  if(e1071::skewness(x,na.rm=TRUE,type=2) < 0) # specify the type 2
    return(exp(c*x))
  else
    return(log(x+c))
}
