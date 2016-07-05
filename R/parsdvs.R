#' Forward selection of DVs
#'
#' @param rv Vector of response variable values.
#' @param dv Dataframe of DVs to be selected from.
#' @param alpha Alpha level for F-test.
#' @param dir Directory to which Maxent runs are written
#'

.parsdvs <- function(rv, dv, alpha, dir) {

  selectedset <- character(length=0)
  remainingset <- colnames(dv)
  evtable <- data.frame()
  cyclenumber <- 0
  bestFVA <- 0

  iterationexit <- FALSE
  while (iterationexit == FALSE) {

    cyclenumber <- cyclenumber + 1
    cycledir <- .dirpath.create(dir, paste0("cycle", cyclenumber))
    cyclemodels <- lapply(remainingset, function(x) c(selectedset, x))

    nrows <- length(cyclemodels)
    ctable <- data.frame(cycle=integer(nrows), model=integer(nrows),
      DV=character(nrows), m=integer(nrows), trainAUC=numeric(nrows),
      Entropy=numeric(nrows), FVA=numeric(nrows), addedFVA=numeric(nrows),
      Fstatistic=numeric(nrows), dfe=integer(nrows), dfu=integer(nrows),
      Pvalue=numeric(nrows), Directory=character(nrows),
      stringsAsFactors = F)

    for (i in 1:length(cyclemodels)) {
      dvnames <- cyclemodels[[i]]
      modeldir <- .dirpath.create(cycledir, paste0("model", i))
      df <- data.frame(dv[, dvnames])
      colnames(df) <- dvnames
      .runjar(rv, df, maxbkg = length(rv) + 1, modeldir)

      maxRes <- read.csv(file.path(modeldir, "maxentResults.csv"))
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
      dfe <- 1
      dfu <- (N - n) - (ctable$m[i] + 1) - 1
      ctable$Fstatistic[i] <- (ctable$addedFVA[i] * dfu) /
        ((1 - ctable$FVA[i]) * dfe)
      ctable$dfe[i] <- dfe
      ctable$dfu[i] <- dfu
      ctable$Pvalue[i] <- 1 - stats::pf(ctable$Fstatistic[i], dfe, dfu)
      ctable$Directory[i] <- modeldir
    }

    ctable <- ctable[order(-ctable$Fstatistic), ]
    evtable <- rbind(evtable, ctable, make.row.names = FALSE)

    if (ctable$Pvalue[1] < alpha) {
      selectedset <- unlist(strsplit(ctable$DV[1], split=" "))
      bestFVA <- ctable$FVA[1]
      addedDV <- sapply(strsplit(ctable$DV[seq(nrow(ctable))[-1]], split=" "),
        function(x) {x[cyclenumber]})
      remainingset <- addedDV[ctable$Pvalue[seq(nrow(ctable))[-1]] < alpha]
    }

    if (nrow(ctable) == 1 || ctable$Pvalue[1] > alpha ||
        (ctable$Pvalue[1] < alpha &&
            all(ctable$Pvalue[seq(nrow(ctable))[-1]] >= alpha))) {
      iterationexit <- TRUE
    }
  }

  return(list(dv[,selectedset, drop = FALSE], evtable))
}