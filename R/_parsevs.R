#' Forward selection of EVs
#'
#' @param rv Vector of response variable values.
#' @param ev List of EVs to be selected from.
#' @param alpha Alpha level for F-test.
#' @param interaction Logical. Allows interaction terms.
#' @param dir Directory to which Maxent runs are written.
#' @param jarpath Pathway to maxent.jar
#'

.parsevs <- function(rv, ev, alpha, interaction, dir, jarpath) {

  selectedset <- character(length=0)
  remainingset <- names(ev)
  modeltable <- data.frame()
  cyclenumber <- 0
  mnull <- 0
  bestFVA <- 0

  iterationexit <- FALSE
  while (iterationexit == FALSE) {

    cyclenumber <- cyclenumber + 1
    cycledir <- paste(dir, paste0("cycle", cyclenumber), sep="\\")
    dir.create(cycledir)
    cyclemodels <- lapply(remainingset, function(x) c(selectedset, x))

    nrows <- length(cyclemodels)
    ctable <- data.frame(cycle=integer(nrows), model=integer(nrows),
      EV=character(nrows), m=integer(nrows), trainAUC=numeric(nrows),
      Entropy=numeric(nrows), FVA=numeric(nrows), addedFVA=numeric(nrows),
      Fstatistic=numeric(nrows), dfe=integer(nrows), dfu=integer(nrows),
      Pvalue=numeric(nrows), Directory=character(nrows),
      stringsAsFactors = F)

    for (i in 1:length(cyclemodels)) { #JV: consider replacing with lapply + function
      evnames <- cyclemodels[[i]]
      dvnames <- unlist(lapply(ev[evnames], names), use.names = FALSE)
      df <- data.frame("RV" = rv, "X" = -9999, "Y" = -9999, ev[evnames])
      colnames(df)[4:ncol(df)] <- dvnames

      modeldir <- paste(cycledir, paste0("model", i), sep = "\\")
      dir.create(modeldir)

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

      maxRes <- read.csv(paste(modeldir, "\\maxentResults.csv", sep=""))
      ctable$cycle[i] <- cyclenumber
      ctable$model[i] <- i
      ctable$EV[i] <- paste(evnames, collapse = " ")
      ctable$m[i] <- length(dvnames)
      ctable$trainAUC[i] <- maxRes$Training.AUC
      ctable$Entropy[i] <- maxRes$Entropy
      n <- maxRes$X.Training.samples
      N <- maxRes$X.Background.points
      ctable$FVA[i] <- (log(N) - ctable$Entropy[i]) / (log(N) - log(n))
      ctable$addedFVA[i] <- ctable$FVA[i] - bestFVA
      dfe <- ctable$m[i] - mnull
      dfu <- (N - n) - (ctable$m[i] + 1) - 1
      ctable$Fstatistic[i] <- (ctable$addedFVA[i] * dfu) /
        ((1 - ctable$FVA[i]) * dfe)
      ctable$dfe[i] <- dfe
      ctable$dfu[i] <- dfu
      ctable$Pvalue[i] <- 1 - pf(ctable$Fstatistic[i], dfe, dfu)
      ctable$Directory[i] <- modeldir
    }

    ctable <- ctable[order(-ctable$Fstatistic), ]
    modeltable <- rbind(modeltable, ctable, make.row.names = FALSE)

    if (ctable$Pvalue[1] < alpha) {
      selectedset <- unlist(strsplit(ctable$EV[1], split=" "))
      mnull <- ctable$m[1]
      bestFVA <- ctable$FVA[1]
      addedEV <- unlist(lapply(strsplit(ctable$EV[2:nrow(ctable)], split=" "),
        function(x) {x[cyclenumber]}))
      remainingset <- addedEV[ctable$Pvalue[2:nrow(ctable)] < alpha]
    }

    if (nrow(ctable) == 1 || ctable$Pvalue[1] > alpha ||
        (ctable$Pvalue[1] < alpha &&
            all(ctable$Pvalue[2:nrow(ctable)] >= alpha))) {
      iterationexit <- TRUE
    }
  }

  return(list(ev[selectedset], modeltable))
}