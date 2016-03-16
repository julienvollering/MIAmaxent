#' calculates skewness of a vector
#'
#' Also calculates the constant 'c' needed for zero-skewness transformation in
#' \code{scalex}
#'
#' \code{DESCRIPTION Imports}: e1071
#'
#' @param x Vector of data. Must have scale [0,1]!


minskew <- function(x) {
  cmin <- min(x)-10*(max(x)-min(x))
  cmax <- max(x)+10*(max(x)-min(x))
  if(e1071::skewness(x, na.rm=TRUE, type=2) >= 0 && cmin < -min(x)) {
    cmin <- -min(x)
  }
  cmid <- (cmin+cmax)/2;
  skew <- e1071::skewness(altrMaxent::.scalex(x, cmid), na.rm=TRUE);
  while (abs(skew) > 1*10^-05 && min(abs(c(cmax, cmin)-cmid)) > 10^-10) {
    #_OptCode print(c(cmin,cmid,cmax,skew));
    sleft <- e1071::skewness(altrMaxent::.scalex(x, (cmid+cmin)/2),
      na.rm=TRUE, type=2)
    sright <- e1071::skewness(altrMaxent::.scalex(x, (cmid+cmax)/2),
      na.rm=TRUE, type=2)
    if (abs(sleft) < abs(skew) && abs(sleft) < abs(sright)) {
      cmax <- cmid
      skew <- sleft
    }
    else if (abs(sright) < abs(skew)) {
      cmin <- cmid
      skew <- sright
    }
    else {
      cmin <- (cmin+cmid)/2;
      cmax <- (cmax+cmid)/2;
    }
    cmid <- (cmin+cmax)/2;
  }
  return(list(c=cmid, skew=skew));
}
