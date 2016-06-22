#' Plot Frequency of Observed Presence (FOP)
#'
#' \code{plotFOP} produces a Frequency of Observed Presence (FOP) plot for a
#' given explanatory variable. For continuous variables, the exponetially
#' weighted moving average of the FOP values is added. \code{plotFOP} also
#' returns a list containing the optimum EV value, and a data frame containing
#' the plotted data (for customizable graphics).
#'
#' The EVoptimum that is retuned is based on the smoothed data, unless a maximum
#' exists at the extremes of the EV (outside the smoothing window). Note that if
#' the response variable represents presence/absence data, the result is an
#' empirical frequency of presence curve, rather than a observed frequency of
#' presence curve (see Stoea et al. [in press], Sommerfeltia).
#'
#' \code{DESCRIPTION Imports}: dplyr, Hmisc
#'
#' @param data Dataframe containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent presence or background, coded as: 1/NA. See Details for
#'   information regarding presence/absence data.
#' @param EV Name or column index of the explanatory variable for which to
#'   calculate FOP.
#' @param smoothwindow Width of the smoothing window. Represents the number of
#'   intervals included in an exponentially weighted moving average. Should be
#'   odd, otherwise the window will be uncentered. Irrelevant for categorical
#'   EVs.
#' @param EVranging if \code{TRUE}, will range the EV scale to [0,1]. This is
#'   equivalent to plotting FOP over the linear transformation produced by
#'   deriveVars. Irrelevant for categorical EVs.
#' @param intervals Number of intervals into which the continuous EV is divided.
#'   Defaults to the minimum of N/50 and 100. Irrelevant for categorical EVs.
#'
#' @return In addition to the graphical output, a list is returned with 1) the
#'   EV value at which FOP is highest (\code{EVoptimum}) and 2) a data frame
#'   containing the plotted data (\code{FOPdata}).
#'
#' @export


plotFOP <- function(data, EV, smoothwindow = 5, EVranging = FALSE,
                    intervals = NULL) {

  df <- data.frame(RV = data[, 1], EV = data[, EV])
  EVname <- colnames(data[, EV, drop = FALSE])

  .binaryrvcheck(df[, 1])
  df[, 1][is.na(df[, 1])] <- 0

  if (class(df[, 2]) %in% c("numeric", "integer")) {
    if (EVranging == T) {
      df[, 2] <- (df[, 2] - min(df[, 2])) / diff(range(df[, 2]))
    }
    if (is.null(intervals)) {intervals <- min(c(ceiling(nrow(df) / 50), 100))}
    df$int <- .reg.interval(df[, 2], intervals)

    grouped <- dplyr::group_by(df, int)
    FOPdf <- as.data.frame(dplyr::summarise(grouped, n = n(), intEV = mean(EV),
      intRV = mean(RV, na.rm=F)))

    FOPdf$smoothRV <- .ewma(FOPdf$intRV, smoothwindow)

    plot(FOPdf$intRV ~ FOPdf$intEV, main = paste0("FOP plot: ", EVname),
      xlab = EVname, ylab = "Frequency of Observed Presence (FOP)")
    lines(FOPdf$intEV, FOPdf$smoothRV, col="grey")

    maxRV <- FOPdf$smoothRV
    maxRV[is.na(maxRV)] <- FOPdf$intRV[is.na(maxRV)]
    EVoptimum = FOPdf$intEV[which(maxRV == max(maxRV))]

    while (length(EVoptimum) > 1) {
      intervals <- intervals - 1
      df$int <- .reg.interval(df[, 2], intervals)
      grouped <- dplyr::group_by(df, int)
      FOPdf <- dplyr::summarise(grouped, intEV = mean(EV),
        intRV = mean(RV, na.rm=F))
      FOPdf$smoothRV <- .ewma(FOPdf$intRV, smoothwindow)
      maxRV <- FOPdf$smoothRV
      maxRV[is.na(maxRV)] <- FOPdf$intRV[is.na(maxRV)]
      EVoptimum <- FOPdf$intEV[which(maxRV == max(maxRV))]
    }

    FOP <- list(EVoptimum = EVoptimum, FOPdata = FOPdf)
  }

  if (class(df[, 2]) %in% c("factor", "character")) {
    grouped <- dplyr::group_by(df, EV)
    FOPdf <- dplyr::summarise(grouped, n = n(), intRV = mean(RV, na.rm=F))

    barplot(FOPdf$intRV, names.arg = FOPdf$EV,
      main = paste0("FOP plot: ", EVname), xlab = EVname,
      ylab = "Frequency of Observed Presence (FOP)")

    FOP <- list(EVoptimum = FOPdf$EV[which(FOPdf$intRV == max(FOPdf$intRV))],
      FOPdata = data.frame(n = FOPdf$n, level = FOPdf$EV,
        freqRV = FOPdf$intRV))
  }

  return(FOP)
}
