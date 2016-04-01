#' Forward selection of DVs
#'
#' @param rv Vector of response variable values.
#' @param dv Dataframe of DVs to be selected from.
#' @param alpha Alpha level for F-test.
#' @param dir Directory to which Maxent runs are written
#' @param jarpath Pathway to maxent.jar
#'

.parsdvs <- function(rv, dv, alpha, dir, jarpath) {

  selectedset <- character(length=0)
  remainingset <- colnames(dv)
  cyclenumber <- 0
  bestFVA <- 0

  iterationexit <- FALSE
  while (iterationexit == FALSE) {

    cyclenumber <- cyclenumber + 1
    cycledir <- paste(dir, paste0("cycle", cyclenumber), sep="\\")
    dir.create(cycledir)
    cyclemodels <- lapply(remainingset, function(x) c(selectedset, x))

    nrows <- length(cyclemodels)
    ctable <- data.frame(cycle=integer(nrows), model=integer(nrows),
      DV=character(nrows), m=integer(nrows), trainAUC=numeric(nrows),
      Entropy=numeric(nrows), FVA=numeric(nrows), addedFVA=numeric(nrows),
      Fstatistic=numeric(nrows), dfe=integer(nrows), dfa=integer(nrows),
      Pvalue=numeric(nrows), Directory=character(nrows),
      stringsAsFactors = F)

    for (i in 1:length(cyclemodels)) { #JV: consider replacing with lapply + function
      dvnames <- cyclemodels[[i]]
      df <- data.frame("RV" = rv, "X" = -9999, "Y" = -9999, dv[,dvnames])
      colnames(df)[4:ncol(df)] <- dvnames

      modeldir <- paste(cycledir, paste0("model", i), sep = "\\")
      dir.create(modeldir)

      samplesdf <- na.omit(df)
      environlayersdf <- df
      csvfiles <- paste(modeldir, c("\\samples.csv", "\\environlayers.csv"), sep="")
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

      maxRes <- read.csv(paste(modeldir, "\\maxentResults.csv", sep=""))
      ctable$cycle[i] <- cyclenumber
      ctable$model[i] <- i
      ctable$DV[i] <- paste(dvnames, collapse = " ")
      ctable$m[i] <- length(dvnames)
      ctable$trainAUC[i] <- maxRes$Training.AUC
      ctable$Entropy[i] <- maxRes$Entropy
      n <- maxRes$X.Training.samples
      N <- maxRes$X.Background.points
      ctable$FVA[i] <- (log(N) - ctable$Entropy[i]) / (log(N) - log(n))
      ctable$addedFVA[i] <- ctable$FVA[i] - bestFVA
      df <- N - n - 3
      ctable$Fstatistic[i] <- (table$FVA[i] * df) / ((1 - table$FVA[i]) * 1)
      ctable$dfe[i] <- # check Sommerfeltia
      ctable$dfa[i] <-
      ctable$Pvalue[i] <-
      ctable$Directory[i] <-
    }

    bestFVA <-
    if (length(addeddv) < 1) {
      iterationexit <- TRUE
    }
  }

#######################################################################


  pb <- txtProgressBar(min = 0, max = ncol(dv), style = 3)

  for (i in 1:ncol(dv)) {
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
    table$DV[i] <- dvname
    table$KnotPosition[i] <- (2 * i - 1) / (2 * ncol(dv))
    table$n[i] <- maxRes$X.Training.samples
    table$N[i] <- maxRes$X.Background.points
    table$Entropy[i] <- maxRes$Entropy
    table$trainingAUC[i] <- maxRes$Training.AUC
    table$FVA[i] <- (log(table$N[i]) - table$Entropy[i]) /
                         (log(table$N[i]) - log(table$n[i]))
    table$df[i] <- table$N[i] - table$n[i] - 3
    table$Fstatistic[i] <- (table$FVA[i] * table$df[i]) /
                                ((1-table$FVA[i]) * 1)
    table$Pvalue[i] <- 1 - pf(table$Fstatistic[i], 1, table$df[i])
    table$Directory[i] <- dvdir

    setTxtProgressBar(pb, i)
  }

  write.csv(table, paste(dir, "\\splineselection.csv", sep=""), row.names = F)

  selected <- character()
  for (i in 3:(nrow(table)-2)) {
    if (table$FVA[i] >= table$FVA[i-2] &&
        table$FVA[i] >= table$FVA[i-1] &&
        table$FVA[i] >= table$FVA[i+1] &&
        table$FVA[i] >= table$FVA[i+2] &&
        table$Pvalue[i] < 0.05) {
      selected <- append(selected, table$DV[i])
    }
  }
  selecteddf <- as.data.frame(dv[, colnames(dv) %in% selected, drop=F])

  return(selecteddf)
}