#' calculate optimum ev value based on frequency of presence
#'
#' The optimum that is retuned is based on the smoothed data, unless a maximum
#' exists at the extremes of the EV (outside the 5-interval smoothing window).
#'
#' \code{DESCRIPTION Imports}: dplyr, Hmisc, scales
#'
#' @param data Dataframe containing the response variable in the first column and
#'   explanatory variables in the second column. The response variable should
#'   represent presence or background, coded as: 1/NA.
#' @param intervals Number of intervals into which the continuous EV is divided.
#'   Irrelevant for categorical EVs.
#' @param smoothwindow Width of the smoothing window (in an exponentially
#'   weighted moving average). Irrelevant for categorical EVs.
#' @param EVranging if \code{TRUE}, will range the EV scale to [0,1]. Irrelevant
#'   for categorical EVs.
#'
#' @return the EV value at which FOP is highest (\code{EVoptimum})


.fopoptimum <- function(data, intervals = 20, smoothwindow = 5,
  EVranging = F) {

  df <- data.frame(RV = data[,1], EV = data[,2])
  altrMaxent:::.binaryrvcheck(df[,1])
  df[,1][is.na(df[,1])] <- 0

  if (class(df[,2]) == "numeric" || class(df[,2]) == "integer") {
    if (EVranging == T) {
      df[,2] <- scales::rescale(df[,2], to = c(0,1))
    }
    intwidth <- (range(df[,2])[2]-range(df[,2])[1])/intervals
    cutpts <- seq(range(df[,2])[1], range(df[,2])[2], by = intwidth)
    df$int <- Hmisc::cut2(df[,2], cuts = cutpts)

    grouped <- dplyr::group_by(df, int)
    FOPdf <- dplyr::summarise(grouped,
      intEV = mean(EV),
      intRV = mean(RV, na.rm=F)
      )

    FOPdf$smoothRV <- altrMaxent:::.ewma(FOPdf$intRV, smoothwindow)
    maxRV <- FOPdf$smoothRV
    maxRV[is.na(maxRV)] <- FOPdf$intRV[is.na(maxRV)]
    EVoptimum <- FOPdf$intEV[which(maxRV == max(maxRV))]

  } else {
    stop("EVoptimum is calculated for numeric or integer class EVs only",
      call. = F)
  }

  return(EVoptimum)
}
