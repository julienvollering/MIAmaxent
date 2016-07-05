#' Plot Frequency of Observed Presence (FOP).
#'
#' \code{plotFOP} produces a Frequency of Observed Presence (FOP) plot for a
#' given explanatory variable. For continuous variables, the exponetially
#' weighted moving average of the FOP values is added. \code{plotFOP} also
#' returns a list containing the optimum EV value, and a data frame containing
#' the plotted data (for customizable graphics).
#'
#' If the response variable in \code{data} represents presence/absence data, the
#' result is an empirical frequency of presence curve, rather than a observed
#' frequency of presence curve (see Stoea et al. [in press], Sommerfeltia).
#'
#' The returned value of 'EVoptimum' is based on the smoothed FOP values, such
#' that an outlying maximum in FOP may, in some cases, not be considered the
#' optimal value of EV.
#'
#' \code{DESCRIPTION Imports}: dplyr, Hmisc
#'
#' @param data Data frame containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent presence or background, coded as: 1/NA. See Details for
#'   information regarding presence/absence data.
#' @param EV Name or column index of the explanatory variable in \code{data} for
#'   which to calculate FOP.
#' @param smoothwindow Width of the smoothing window. Represents the number of
#'   intervals included in an exponentially weighted moving average. Should be
#'   odd, otherwise the window will be uncentered. Irrelevant for categorical
#'   EVs.
#' @param EVranging Logical. If \code{TRUE}, will range the EV scale to [0,1].
#'   This is equivalent to plotting FOP over the linear transformation produced
#'   by deriveVars. Irrelevant for categorical EVs.
#' @param intervals Number of intervals into which the continuous EV is divided.
#'   Defaults to the minimum of N/50 and 100. Irrelevant for categorical EVs.
#' @param ... Arguments to be passed to \code{plot} to control the appearance of
#'   the plot. For example: \itemize{ \item \code{cex} for size of points \item
#'   \code{col} for color \item \code{xlim} for range of the x-axis }
#'
#' @return In addition to the graphical output, a list of 2: \enumerate{ \item
#'   The EV value at which FOP is highest (\code{EVoptimum}) \item a data frame
#'   containing the plotted data (\code{FOPdata}).}
#'
#' @export


plotFOP <- function(data, EV, smoothwindow = 5, EVranging = FALSE,
                    intervals = NULL, ...) {

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
    FOPdf <- as.data.frame(dplyr::summarise(grouped, n = n(),
      intEV = mean(EV), intRV = mean(RV, na.rm=F)))

    FOPdf$smoothRV <- .ewma(FOPdf$intRV, smoothwindow)

    graphics::plot(FOPdf$intRV ~ FOPdf$intEV,
      main = paste0("FOP plot: ", EVname),
      xlab = EVname, ylab = "Frequency of Observed Presence (FOP)", ...)
    graphics::lines(FOPdf$intEV, FOPdf$smoothRV, col="grey")

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

    graphics::barplot(FOPdf$intRV, names.arg = FOPdf$EV,
      main = paste0("FOP plot: ", EVname), xlab = EVname,
      ylab = "Frequency of Observed Presence (FOP)", ...)

    FOP <- list(EVoptimum = FOPdf$EV[which(FOPdf$intRV == max(FOPdf$intRV))],
      FOPdata = data.frame(n = FOPdf$n, level = FOPdf$EV,
        freqRV = FOPdf$intRV))
  }

  return(FOP)
}
