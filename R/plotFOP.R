#' plot Frequency of Observed Presence
#'
#' \code{plotFOP} produces a Frequency of Observed Presence plot for a given
#' explanatory variable, and provides a data frame with the plotted data (for
#' customizable plotting).The EVoptimum that is retuned is based on the smoothed
#' data, unless a maximum exists at the extrames of the EV (outside the smoothing
#' window).
#'
#' @param data dataframe. Containing the response variable in the first column
#' and explanatory variables in subsequenct columns.
#'
#' Imports: dplyr, Hmisc


plotFOP <- function(data, EV, intervals = 20, smoothwindow = 3) {
  df <- data.frame(RV = data[,1], EV = data[,EV])

  if (anyNA(df[,1]) && length(levels(as.factor(df[,1]))) > 1) {
    stop("The response variable must contain 2 levels only: presence (1)
      and background (NA/0)", call. = FALSE)
  }

  if (class(df[,1]) != "numeric" && class(df[,1]) != "integer") {
    stop("The response variable must be numeric or integer class: presence (1)
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

  s <- smoothwindow
  FOPdf$smoothRV <- ewma(FOPdf$intRV, s)

  RVname <- colnames(data)[1]
  if (class(EV)=="character") {
    EVname <- EV
  } else {
    EVname <- colnames(data)[EV]
  }
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

  # need to add support for categorical EVs

  return(FOP)
}
