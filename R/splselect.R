#' selects a subset of spline dvs based on Maxent FTVE
#'
#' @param rv Vector of response variable values
#' @param dv List of spline dvs to be selected from (HF, HR, or Th)
#' @param dir Directory to which Maxent runs are written
#' @param jarpath Pathway to maxent.jar
#'

.splselect <- function(rv, dv, dir, jarpath) {

  n <- length(dv)
  comparison <- data.frame(DV=character(n), KnotPosition=numeric(n),
    n=integer(n), N=integer(n), Entropy=numeric(n), trainingAUC=numeric(n),
    FVA=numeric(n), df=integer(n), Fstatistic=numeric(n), Pvalue=numeric(n),
    Directory=character(n), stringsAsFactors = F)

  pb <- txtProgressBar(min = 0, max = n, style = 3)

  for (i in 1:n) {
    dvname <- names(dv)[[i]]
    df <- data.frame("RV" = rv, "X" = -9999, "Y" = -9999, dv[[i]])
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
    comparison$KnotPosition[i] <- (2 * i - 1) / (2 * n)
    comparison$n[i] <- maxRes$X.Training.samples
    comparison$N[i] <- maxRes$X.Background.points
    comparison$Entropy[i] <- maxRes$Entropy
    comparison$trainingAUC[i] <- maxRes$Training.AUC
    comparison$FVA[i] <- (log(comparison$N[i]) - comparison$Entropy[i]) /
                         (log(comparison$N[i]) - log(comparison$n[i]))
    comparison$df[i] <- comparison$N[i] - comparison$n[i] - 3
    comparison$Fstatistic[i] <- (comparison$FVA[i] * comparison$df[i]) /
                                ((1-comparison$FVA[i]) * 1)
    comparison$Pvalue[i] <- 1 - pf(comparison$Fstatistic[i], 1, comparison$df[i])
    comparison$Directory[i] <- dvdir

    setTxtProgressBar(pb, i)
  }

  write.csv(comparison, paste(dir, "\\splineselection.csv", sep=""), row.names = F)

  selected <- character()
  for (i in 3:(nrow(comparison)-2)) {
    if (comparison$FVA[i] >= comparison$FVA[i-2] &&
        comparison$FVA[i] >= comparison$FVA[i-1] &&
        comparison$FVA[i] >= comparison$FVA[i+1] &&
        comparison$FVA[i] >= comparison$FVA[i+2] &&
        comparison$Pvalue[i] < 0.05) {
      selected <- append(selected, comparison$DV[i])
    }
  }

  ptsx <- comparison$KnotPosition[names(dv) %in% selected]
  ptsy <- comparison$FVA[names(dv) %in% selected]

  png(paste(dir, "\\Vknotplot.png", sep=""))
  plot(comparison$KnotPosition, comparison$FVA, lty = "solid",
    main = "V-knot plot",
    xlab = "Position of knot",
    ylab = "Fraction of variation accounted for (FVA)")
  if (length(selected) > 0) {
    points(ptsx, ptsy, col="red", pch=16)
    text(ptsx, ptsy, labels=selected, cex= 0.9, pos=1)
  }

  dev.off()

  return(selected)
}