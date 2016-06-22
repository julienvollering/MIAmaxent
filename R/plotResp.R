#' Plot model response across a given explanatory variable.
#'
#' \code{plotResp} plots  the response of a given Maxent model over any of the
#' included explanatory variables (EVs) in that model. For categorical
#' variables, a box plot is returned rather than a curve. \code{plotResp} also
#' returns a data frame containing the plotted data (for customizable graphics).
#'
#' The response curves generated in this function are single-effect response
#' curves, presenting the response of a model containing the explanatory
#' variable of interest only (cf. marginal-effect response curves).
#'
#' The \code{ev} specified in \code{dvdata} must not be an interaction term.
#'
#' @param data Data frame of training data, with response variable (1/NA) in the
#'   first column and explanatory variables in subsequenct columns.
#' @param ev Name or list index of the explanatory variable in \code{dvdata} for
#'   which the response curve is to be generated. Interaction terms not allowed.
#' @param dvdata List of derived variables used to train the model, with each
#'   list item a data frame containing 1 or more DVs for a given EV. E.g. output
#'   [[1]] of \code{selectEV}.
#' @param dir The directory to which Maxent files will be written. Defaults to
#'   the working directory.
#' @param jarpath The pathway to the maxent.jar executable jar file. If
#'   unspecified, the function looks for the file in \code{dir}.
#'
#' @return In addition to the graphical output, a data frame containing the
#'   plotted data is returned.
#'
#' @export


plotResp <- function(data, ev, dvdata, dir = NULL, jarpath = NULL,
                     logscale = FALSE) {

  .binaryrvcheck(data[, 1])

  if (is.null(dir)) {dir <- getwd()}
  if (is.null(jarpath)) {jarpath <- file.path(dir, "maxent.jar")}
  if (file.exists(jarpath) == F) {
    stop("maxent.jar file must be present in dir, or its pathway must be
       specified by the jarpath argument. \n ", call. = FALSE)
  }

  dir <- file.path(dir, "plotResp")
  dir.create(dir, showWarnings = FALSE)
  evname <- names(dvdata[ev])
  if (!(evname %in% colnames(data))) {
    stop("The specified EV must be present in the untransformed data.
       E.g. interaction terms between multiple EVs are not supported. \n ",
      call. = FALSE)
  }
  modeldir <- file.path(dir, paste0("response", evname))
  if (file.exists(modeldir)) {
    stop("The response to this EV has already been evaluated in the given dir.
       Please select a different EV or specify a different dir. \n ",
      call. = FALSE)
  } else {
    dir.create(modeldir)
  }

  df <- data.frame("RV" = data[, 1], "X" = -9999, "Y" = -9999, dvdata[[ev]])
  samplesdf <- na.omit(df)
  environlayersdf <- df
  csvfiles <- file.path(modeldir, c("samples.csv", "environlayers.csv"))
  write.csv(samplesdf, csvfiles[1], row.names = F)
  write.csv(environlayersdf, csvfiles[2], row.names = F)

  jarflags1 <- " removeduplicates=FALSE addsamplestobackground=FALSE"
  jarflags2 <- " maximumbackground=100000 autofeature=FALSE betamultiplier=0"
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
    jarflags,
    " samplesfile=","\"", csvfiles[1], "\"",
    " environmentallayers=", "\"", csvfiles[2], "\"",
    " outputdirectory=", "\"", modeldir, "\\", "\"")
  javacommand <- gsub("\\\\","/", command)
  system(paste(javacommand), wait = TRUE)

  output <- read.csv(file.path(modeldir, "1_backgroundPredictions.csv"))
  respPts <- data.frame(data[,evname], output[,3]*length(output[,3]))
  colnames(respPts) <- c("EV", "PRO")
  if (logscale == TRUE) {respPts$PRO <- log10(respPts$PRO)}

  if (class(respPts[, 1]) %in% c("numeric", "integer")) {
    plot(respPts[, 2] ~ respPts[, 1], pch = 20, col="grey",
      main = paste0("Single-effect model response to ", evname), xlab = evname,
      ylab = ifelse(logscale == TRUE, "log Probability Ratio Output (logPRO)",
        "Probability Ratio Output (PRO)"))

    intervals <- min(c(ceiling(nrow(respPts) / 50), 100))
    intwidth <- (max(respPts[, 1]) - min(respPts[, 1])) / intervals
    cutpts <- seq(min(respPts[, 1]), max(respPts[, 1]), by = intwidth)
    respPts$int <- Hmisc::cut2(respPts[, 1], cuts = cutpts, oneval = FALSE)
    grouped <- dplyr::group_by(respPts, int)
    respLine <- as.data.frame(dplyr::summarise(grouped, n = n(),
      intEV = mean(EV, na.rm = TRUE),
      intPRO = mean(PRO, na.rm = TRUE)))
    respLine$smoothPRO <- .ewma(respLine$intPRO, 3)
    lines(respLine$smoothPRO ~ respLine$intEV, col="red", lwd = 2)

    if (logscale == TRUE) {abline(h = 0, lty = 3)} else {abline(h = 1, lty = 3)}

    result <- list("respPts" = respPts, "respLine" = respLine)
  }

  if (class(respPts[, 1]) %in% c("factor", "character")) {
    plot(respPts[, 2] ~ respPts[, 1],
      main = paste0("Single-effect model response to ", evname), xlab = evname,
      ylab = ifelse(logscale == TRUE, "log Probability Ratio Output (logPRO)",
        "Probability Ratio Output (PRO)"))
    if (logscale == TRUE) {abline(h = 0, lty = 3)} else {abline(h = 1, lty = 3)}
    result <- respPts
  }

  return(result)
}
