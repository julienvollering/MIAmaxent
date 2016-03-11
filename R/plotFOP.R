#' plot Frequency of Observed Presence
#'
#' \code{plotFOP} produces a Frequency of Observed Presence plot for a given
#' explanatory variable, and provides a data frame with the plotted data (for
#' customizable plotting).
#'
#' @param data dataframe. Containing the response variable in the first column
#' and explanatory variables in subsequenct columns.
#'
#' Imports: dplyr, Hmisc, stats


plotFOP <- function(data, EV, intervals = 20, smoothwindow = 5) {
  df <- data.frame(RV = data[,1], EV = data[,EV])

  if (anyNA(df[,1]) && length(levels(as.factor(df[,1]))) > 1) {
    stop("The response variable must contain 2 levels only: presence (1)
      and background (NA/0)", call. = FALSE)
  }

  df[,1][is.na(df[,1])] <- 0

  if (class(df[,2]) == "numeric" | class(df[,2]) == "integer") {
  intwidth <- (range(df[,2])[2]-range(df[,2])[1])/intervals
  cutpts <- seq(range(df[,2])[1], range(df[,2])[2], by = intwidth)
  df$int <- Hmisc::cut2(df[,2], cuts = cutpts)

  grouped <- dplyr::group_by(df, int)
  FOPdf <- dplyr::summarise(grouped,
    n = n(),
    intEV = mean(EV),
    intRV = mean(RV, na.rm=F)
    )

  plot(FOPdf$intRV ~ FOPdf$intEV)

  s <- smoothwindow
  if (s %% 2 != 0) {
    expwindow <- dexp(c(((s-1)/2):0,1:((s-1)/2)))
  } else {
    expwindow <- dexp(c((s/2):0,1:((s-2)/2)))
  }
  expweights <- expwindow/sum(expwindow)
  FOPdf$smoothRV <- filter(FOPdf$intRV,expweights, sides=2)
  lines(FOPdf$intEV,FOPdf$smoothRV)
  }

}
