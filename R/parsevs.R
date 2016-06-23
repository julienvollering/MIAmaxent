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

    for (i in 1:length(cyclemodels)) {
      evnames <- cyclemodels[[i]]
      dvnames <- unlist(lapply(ev[evnames], names), use.names = FALSE)
      modeldir <- file.path(cycledir, paste0("model", i))
      dir.create(modeldir)
      df <- data.frame(ev[evnames])
      colnames(df) <- dvnames
      .runjar(rv, df, maxbkg = length(rv) + 1, modeldir, jarpath)

      maxRes <- read.csv(file.path(modeldir, "maxentResults.csv"))
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
      addedEV <- sapply(strsplit(ctable$EV[seq(nrow(ctable))[-1]], split=" "),
        function(x) {x[cyclenumber]})
      remainingset <- addedEV[ctable$Pvalue[seq(nrow(ctable))[-1]] < alpha]
    }

    message(paste0("Cycle ", cyclenumber, " complete."))

    if (nrow(ctable) == 1 || ctable$Pvalue[1] > alpha ||
        (ctable$Pvalue[1] < alpha &&
            all(ctable$Pvalue[seq(nrow(ctable))[-1]] >= alpha))) {
      iterationexit <- TRUE
    }
  }

  if (interaction == FALSE) {
    return(list(ev[selectedset], modeltable))
  }

  message(paste0("Forward selection of interaction terms between ",
    length(selectedset), " EVs"))

  combos <- t(combn(selectedset, 2))
  products <- vector("list", nrow(combos))
  names(products) <- apply(combos, 1, function(x) {paste(x, collapse=":")})
  for (i in 1:nrow(combos)) {
    ev1 <- ev[[combos[i,1]]]
    ev2 <- ev[[combos[i,2]]]
    dvcombos <- expand.grid(names(ev1), names(ev2))
    productdvs <- matrix(nrow = nrow(ev[[1]]), ncol = nrow(dvcombos))
    for (j in 1:nrow(dvcombos)) {
      productdvs[, j] <- ev1[[dvcombos[j,1]]] * ev2[[dvcombos[j,2]]]
    }
    colnames(productdvs) <- apply(dvcombos, 1, function(x) {
      paste(x, collapse=":")})
    products[[i]] <- as.data.frame(productdvs)
  }

  ev <- append(ev, products)
  remainingset <- names(products)

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

    for (i in 1:length(cyclemodels)) {
      evnames <- cyclemodels[[i]]
      dvnames <- unlist(lapply(ev[evnames], names), use.names = FALSE)
      modeldir <- file.path(cycledir, paste0("model", i))
      dir.create(modeldir)
      df <- data.frame(ev[evnames])
      colnames(df) <- dvnames
      .runjar(rv, df, maxbkg = length(rv) + 1, modeldir, jarpath)

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
      addedEV <- sapply(strsplit(ctable$EV[seq(nrow(ctable))[-1]], split=" "),
        function(x) {x[cyclenumber]})
      remainingset <- addedEV[ctable$Pvalue[seq(nrow(ctable))[-1]] < alpha]
    }

    message(paste0("Cycle ", cyclenumber, " complete."))

    if (nrow(ctable) == 1 || ctable$Pvalue[1] > alpha ||
        (ctable$Pvalue[1] < alpha &&
            all(ctable$Pvalue[seq(nrow(ctable))[-1]] >= alpha))) {
      iterationexit <- TRUE
    }
  }

  return(list(ev[selectedset], modeltable))

}