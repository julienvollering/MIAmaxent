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
    dvdir <- .dirpath.create(dir, dvname)
    df <- data.frame(dv[[i]])
    colnames(df) <- dvname
    .runjar(rv, df, maxbkg = length(rv) + 1, dvdir, jarpath)

    maxRes <- read.csv(file.path(dvdir, "maxentResults.csv"))
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