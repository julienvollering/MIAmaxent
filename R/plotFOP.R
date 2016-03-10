#' plot Frequency of Observed Presence
#'
#' \code{plotFOP} produces a Frequency of Observed Presence plot for a given
#' explanatory variable, and provides a data frame with the plotted data (for
#' customizable plotting).
#'
#' @param data dataframe. Containing the response variable in the first column
#' and explanatory variables in subsequenct columns.
#'
#' Imports: dplyr, Hmisc


plotFOP <- function(data, EV, intervals = 20) {
  df <- data.frame(RV = data[,1], EV = data[,EV])
  df[,1][is.na(df[,1])] <- 0

  intwidth <- (range(df[,2])[2]-range(df[,2])[1])/intervals
  cutpts <- seq(range(df[,2])[1], range(df[,2])[2], by = intwidth)
  df$int <- Hmisc::cut2(df[,2], cuts = cutpts)

  binned <- dplyr::group_by(df, int)
  FOP <- dplyr::summarise(binned,
    count = n(),
    fopEV = mean(EV),
    fopRV = mean(RV, na.rm=F)
    )

  plot(FOP$fopRV ~ FOP$fopEV)
}