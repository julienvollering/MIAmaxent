#' Plot Frequency of Observed Presence (FOP)
#'
#' \code{plotFOP} produces a Frequency of Observed Presence (FOP) plot for a
#' given explanatory variable. For continuous variables, the exponetially
#' weighted moving average of the FOP values is added. \code{plotFOP} also
#' returns a list containing the optimum EV value, and a data frame with the
#' plotted data (for customizable plotting).
#'
#' The EVoptimum that is retuned is based on the smoothed data, unless a maximum
#' exists at the extremes of the EV (outside the smoothing window). Note that if
#' the response variable represents presence/absence data, the result is an
#' empirical frequency of presence curve, rather than a observed frequency of
#' presence curve (see Stoea et al. [in press], Sommerfeltia).
#'
#' \code{DESCRIPTION Imports}: dplyr, Hmisc, scales
#'
#' @param data Dataframe containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent presence or background, coded as: 1/NA. See Details for
#'   information regarding presence/absence data.
#' @param EV Name or column number of the explanatory variable for which to
#'   calculate FOP.
#' @param intervals Number of intervals into which the continuous EV is divided.
#'   Irrelevant for categorical EVs.
#' @param smoothwindow Width of the smoothing window. Represents the number of
#'   intervals included in an exponentially weighted moving average. Should be
#'   odd, otherwise the window will be uncentered. Irrelevant for categorical
#'   EVs.
#' @param EVranging if \code{TRUE}, will range the EV scale to [0,1]. This is
#'   equivalent to plotting FOP over the linear transformation produced by
#'   deriveVars. Irrelevant for categorical EVs.
#'
#' @return In addition to the plotted output, a list is returned containing 1)
#'   the EV value at which FOP is highest (\code{EVoptimum}) and 2) a data frame
#'   with the plotted data (\code{FOPdata}).
#'
#' @export


plotFOP <- function(data, EV, intervals = 20, smoothwindow = 3,
                    EVranging = FALSE) {

  df <- data.frame(RV = data[,1], EV = data[,EV])
  RVname <- colnames(data)[1]
  if (class(EV)=="character") {
    EVname <- EV
  } else {
    EVname <- colnames(data)[EV]
  }

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
      n = n(),
      intEV = mean(EV),
      intRV = mean(RV, na.rm=F)
      )

    FOPdf$smoothRV <- altrMaxent:::.ewma(FOPdf$intRV, smoothwindow)

    plot(FOPdf$intRV ~ FOPdf$intEV,
      xlab = EVname, ylab = RVname,
      main = "Frequency of Observed Presence")
    lines(FOPdf$intEV,FOPdf$smoothRV, col="grey")

    maxRV <- FOPdf$smoothRV
    maxRV[is.na(maxRV)] <- FOPdf$intRV[is.na(maxRV)]

    FOP <- list()
    FOP[[1]] <- FOPdf$intEV[which(maxRV == max(maxRV))]
    FOP[[2]] <- data.frame("n"=FOPdf$n, "meanEV"=FOPdf$intEV,
      "freqRV"=FOPdf$intRV, "smooth_freqRV"=FOPdf$smoothRV)
    names(FOP) <- c("EVoptimum", "FOPdata")

  }

  if (class(df[,2]) == "factor" || class(df[,2]) == "character") {
    grouped <- dplyr::group_by(df, EV)
    FOPdf <- dplyr::summarise(grouped,
      n = n(),
      intRV = mean(RV, na.rm=F)
    )

    barplot(FOPdf$intRV, names.arg = FOPdf$EV,
      xlab = EVname, ylab = RVname,
      main = "Frequency of Observed Presence")

    FOP <- list()
    FOP[[1]] <- FOPdf$EV[which(FOPdf$intRV == max(FOPdf$intRV))]
    FOP[[2]] <- data.frame("n"=FOPdf$n, "level"=FOPdf$EV,
      "freqRV"=FOPdf$intRV)
    names(FOP) <- c("EVoptimum", "FOPdata")
  }

  return(FOP)
}
