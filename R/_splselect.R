#' selects a subset of spline dvs based on Maxent FTVE
#'
#' @param rv Vector of response variable values
#' @param dv Dataframe of spline dvs to be selected from (HF, HR, or Th)
#' @param dir Directory to which Maxent runs are
#'   written
#' @param jarpath Pathway to maxent.jar
#'

.splselect <- function(rv, dv, dir, jarpath) {

  comparison <- data.frame(DV=character(), n=integer(), N=integer(),
                Entropy=numeric(), trainingAUC=numeric(), FTA=numeric(),
                Directory=character())

  for (i in 1:ncol(dv))
    dvname <- colnames(dv)[i]
    df <- data.frame("RV" = rv, "X" = -9999, "Y" = -9999, dv[,i])
    colnames(df)[4] <- dvname

    dvdir <- paste(dir, "\\", dvname, sep = "")
    dir.create(dvdir)

    samplesdf <- na.omit(df)
    environlayersdf <- df

    csvfiles <- paste(dvdir, c("\\samples.csv", "\\environlayers.csv"), sep="")
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
    jarflags <- paste(jarflags1, jarflags2, jarflags3, jarflags4, jarflags5,
      jarflags6, jarflags7, jarflags8, sep="")

    command <- paste("java -mx512m -jar ",
                     "\"", jarpath, "\"",
                     jarflags,
                     " samplesfile=","\"", csvfiles[1], "\"",
                     " environmentallayers=", "\"", csvfiles[2], "\"",
                     " outputdirectory=", "\"", dvdir, "\\", "\"",
                     sep="")
    javacommand <- gsub("\\\\","/", command)
    system(paste(javacommand), wait=T)

    maxRes <- read.csv(paste(dvdir, "\\maxentResults.csv", sep=""))
    comparison$DV[i] <- dvname
    comparison$n[i] <- maxRes$X.Training.samples
    comparison$N[i] <- maxRes$X.Background.points
    comparison$Entropy[i] <- maxRes$Entropy
    comparison$trainingAUC[i] <- maxRes$Training.AUC
    comparison$FTA[i] <- (log(comparison$N[i]) - comparison$Entropy[i]) /
                         (log(comparison$N[i]) - log(comparison$n[i]))
    comparison
}