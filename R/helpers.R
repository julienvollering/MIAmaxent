#' Executes basic maxent.jar run from R
#'
#' @param rv Vector of response variable values
#' @param ev Data frame of explanatory variables
#' @param maxbkg Maximum number of uninformed background points to use for
#'   training
#' @param dir Directory to which Maxent files will be written
#' @param jarpath Pathway to the maxent.jar executable jar file

.runjar <- function(rv, ev, maxbkg = 10000, dir, jarpath) {
  df <- data.frame("RV" = rv, "X" = -9999, "Y" = -9999, ev, check.names = FALSE)
  samplesdf <- na.omit(df)
  environlayersdf <- df
  csvfiles <- file.path(dir, c("samples.csv", "environlayers.csv"))
  write.csv(samplesdf, csvfiles[1], row.names = F)
  write.csv(environlayersdf, csvfiles[2], row.names = F)

  jarflags1 <- " removeduplicates=FALSE addsamplestobackground=FALSE"
  jarflags2 <- " autofeature=FALSE betamultiplier=0"
  jarflags3 <- " quadratic=FALSE product=FALSE hinge=FALSE threshold=FALSE"
  jarflags4 <- " outputformat=raw writebackgroundpredictions=TRUE"
  jarflags5 <- " outputgrids=FALSE pictures=FALSE"
  jarflags6 <- " extrapolate=FALSE writemess=FALSE plots=FALSE"
  jarflags7 <- " doclamp=FALSE writeclampgrid=FALSE"
  jarflags8 <- " autorun=TRUE threads=8 visible=FALSE warnings=FALSE"
  jarflags <- paste0(jarflags1, jarflags2, jarflags3, jarflags4, jarflags5,
    jarflags6, jarflags7, jarflags8)

  command <- paste0("java -mx512m -jar ",
    "\"", jarpath, "\"",
    jarflags, " maximumbackground=", maxbkg,
    " samplesfile=","\"", csvfiles[1], "\"",
    " environmentallayers=", "\"", csvfiles[2], "\"",
    " outputdirectory=", "\"", dir, "\\", "\"")
  javacommand <- gsub("\\\\","/", command)
  system(paste(javacommand), wait = TRUE)
}


#' calculates optimum EV value based on FOP
#'
#' The optimum that is returned is based on the smoothed data, unless a maximum
#' exists at the extremes of the EV (outside the 5-interval smoothing window).
#'
#' @param data Dataframe containing the response variable in the first column and
#'   explanatory variables in the second column. The response variable should
#'   represent presence or background, coded as: 1/NA.
#' @param smoothwindow Width of the smoothing window (in an exponentially
#'   weighted moving average). Irrelevant for categorical EVs.
#' @param intervals Number of intervals into which the continuous EV is divided.
#'   Defaults to the minimum of N/50 and 100. Irrelevant for categorical EVs.
#'
#' @return the EV value at which FOP is highest (\code{EVoptimum})

.fopoptimum <- function(data, smoothwindow = 5, intervals = NULL) {

  df <- data.frame(RV = data[, 1], EV = data[, 2])
  .binaryrvcheck(df[, 1])
  df[, 1][is.na(df[, 1])] <- 0

  if (!class(df[, 2]) %in% c("numeric", "integer")) {
    stop("EVoptimum is calculated for numeric or integer class EVs only",
      call. = F)
  }

  if (is.null(intervals)) {intervals <- min(c(ceiling(nrow(df) / 50), 100))}
  df$int <- .reg.interval(df[, 2], intervals)

  grouped <- dplyr::group_by(df, int)
  FOPdf <- dplyr::summarise(grouped, intEV = mean(EV), intRV = mean(RV, na.rm=F))

  FOPdf$smoothRV <- .ewma(FOPdf$intRV, smoothwindow)
  maxRV <- FOPdf$smoothRV
  maxRV[is.na(maxRV)] <- FOPdf$intRV[is.na(maxRV)]
  EVoptimum <- FOPdf$intEV[which(maxRV == max(maxRV))]

  while (length(EVoptimum) > 1) {
    intervals <- intervals - 1
    df$int <- .reg.interval(df[, 2], intervals)
    grouped <- dplyr::group_by(df, int)
    FOPdf <- dplyr::summarise(grouped, intEV = mean(EV), intRV = mean(RV, na.rm=F))
    FOPdf$smoothRV <- .ewma(FOPdf$intRV, smoothwindow)
    maxRV <- FOPdf$smoothRV
    maxRV[is.na(maxRV)] <- FOPdf$intRV[is.na(maxRV)]
    EVoptimum <- FOPdf$intEV[which(maxRV == max(maxRV))]
  }

  return(EVoptimum)
}