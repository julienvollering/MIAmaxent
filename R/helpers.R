#' Executes basic maxent.jar run from R
#'
#' @param rv Vector of response variable values
#' @param ev Data frame of explanatory variables
#' @param maxbkg Maximum number of uninformed background points to use for
#'   training
#' @param dir Directory to which Maxent files will be written
#' @param jarpath Pathway to the maxent.jar executable jar file

.runjar <- function(rv, ev, maxbkg = 10000, dir, jarpath) {
  df <- data.frame("RV" = rv, "X" = -9999, "Y" = -9999, ev)
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
