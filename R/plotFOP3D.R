#' Plot Frequency of Observed Presence (FOP) in 3D
#'
#' \code{plotFOP3D} does ...
#'
#' \code{DESCRIPTION Imports}: dplyr, Hmisc, scatterplot3d
#'
#' @param data Data frame containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent presence or background, coded as: 1/NA. See Details for
#'   information regarding presence/absence data.
#' @param EV Name or column index of the explanatory variable in \code{data} for
#'   which to calculate FOP.
#' @param EVranging Logical. If \code{TRUE}, will range the EV scale to [0,1].
#'   This is equivalent to plotting FOP over the linear transformation produced
#'   by deriveVars. Irrelevant for categorical EVs.
#' @param intervals Number of intervals into which the continuous EV is divided.
#'   Defaults to the minimum of N/50 and 100. Irrelevant for categorical EVs.
#'
#'
#'


plotFOP <- function(data, EV1, EV2, smoothwindow = 5, EVranging = FALSE) {

  .binaryrvcheck(data[, 1])
  data[, 1][is.na(data[, 1])] <- 0

  df <- data.frame(RV = data[, 1], EV1 = data[, EV1], EV2 = data[, EV2])
  EVnames <- colnames(data[, c(EV1, EV2), drop = FALSE])

  if (EVranging == TRUE) {
    for (i in c(2:3)) {
      if (class(df[, 1]) %in% c("numeric", "integer")) {
        df[, i] <- (df[, i] - min(df[, i])) / diff(range(df[, i]))
      }
    }
  }

  nintervals <- min(c(ceiling(nrow(df) / 50), 100))
  df$int1 <- .reg.interval(df[, 2], nintervals)
  df$int2 <- .reg.interval(df[, 3], nintervals)

  gdf <- dplyr::group_by(df, int1, int2)
  FOPdf <- as.data.frame(dplyr::summarise(gdf, n = n(), intEV1 = mean(EV1),
    intEV2 = mean(EV2), intRV = mean(RV, na.rm=F)))
  par(mfrow=c(2,2))
  for (i in c(36, 144, 72, 108)) {
    scatterplot3d::scatterplot3d(FOPdf$intEV1, FOPdf$intEV2, FOPdf$intRV,
      highlight.3d = TRUE, col.axis = "blue", col.grid = "lightblue",
      pch = 20, type = "h", xlab = EVnames[1], ylab = EVnames[2], zlab = "FOP",
      angle = i)
  }
}
