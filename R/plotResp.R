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
#' @param rv Response variable vector used to train the model. The RV should
#'   represent presence/background data, coded as: 1/NA.
#' @param ev Name or list index of the explanatory variable in \code{dvdata} for
#'   which the response curve is to be generated.
#' @param evdata Explanatory variables used to train the model, in their
#'   original (untransformed) form. Data frame.
#' @param dvdata Derived variables used to train the model. Named list of data
#'   frames, with each data frame containing 1 or more DVs for a given EV. E.g.
#'   output [[1]] of \code{selectEV}.
#' @param writedir The directory to which Maxent files will be written during
#'   subset selection of DVs. Defaults to the working directory.
#' @param jarpath The pathway to the maxent.jar executable jar file. If
#'   unspecified, the function looks for the file in the writedir.
#'
#' @return In addition to the graphical output, a data frame containing the
#'   plotted data is returned.
#'
#' @export


plotResp <- function(rv, ev, evdata, dvdata, writedir = NULL, jarpath = NULL) {

  altrMaxent:::.binaryrvcheck(rv)

  if (is.null(writedir)) {
    writedir <- getwd()
  }

  if (is.null(jarpath)) {
    jarpath <- paste(writedir, "\\maxent.jar", sep="")
  }

  if (file.exists(jarpath) == F) {
    stop("maxent.jar file must be present in writedir, or its pathway must be
       specified by the jarpath argument. \n ", call. = FALSE)
  }

  dir <- paste(writedir, "\\plotResp", sep="")
  dir.create(dir, showWarnings = FALSE)
  evname <- names(dvdata[ev])
  modeldir <- paste(dir, "\\response", evname, sep="")
  if (file.exists(modeldir)) {
    stop("The response to this EV has already been evaluated in the given writedir.
       Please select a different EV or specify a different writedir. \n ",
      call. = FALSE)
  } else {
    dir.create(modeldir)
  }

  df <- data.frame("RV" = rv, "X" = -9999, "Y" = -9999, dvdata[[ev]])
  samplesdf <- na.omit(df)
  environlayersdf <- df
  csvfiles <- paste0(modeldir, c("\\samples.csv", "\\environlayers.csv"))
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
  system(paste(javacommand), wait=T)

  output <- read.csv(paste(modeldir, "\\1_backgroundPredictions.csv", sep=""))
  plotdf <- data.frame(evdata[,evname], output[,3]*length(output[,3]))
  colnames(plotdf) <- c(evname, "PRO")
  plotdf <- plotdf[order(plotdf[,evname]),]

  plot(plotdf[,"PRO"] ~ plotdf[,evname], pch = 20,
    xlab = evname, ylab = "Probability Ratio Output (PRO)",
    main = paste0("Single-effect response to ", evname))
  if (class(plotdf[,evname]) %in% c("numeric", "integer")) {
    lines(plotdf[,"PRO"] ~ plotdf[,evname], col="red", lwd = 2)
  }

  return(plotdf)

}
