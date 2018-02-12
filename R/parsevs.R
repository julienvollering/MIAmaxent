#' Forward selection of EVs
#'
#' @param rv Vector of response variable values.
#' @param ev List of EVs to be selected from.
#' @param alpha Alpha level for F-test.
#' @param test Character string matching either "Chisq" or "F".
#' @param interaction Logical. Allows interaction terms.
#' @param formula Model formula specifying a starting point for model selection.
#'

.parsevs <- function(rv, ev, alpha, interaction, dir, formula) {

  if (!is.null(formula) && length(labels(stats::terms(formula))) != 0) {
    roundnumber <- 0
    mnull <- 0
    bestFVA <- 0
    rounddir <- .dirpath.create(dir, "round0")
    roundmodel <- labels(stats::terms(formula))
    ctable <- data.frame(round=integer(1), model=integer(1), EV=character(1),
      m=integer(1), trainAUC=numeric(1), Entropy=numeric(1), FVA=numeric(1),
      addedFVA=numeric(1), dfe=integer(1), dfu=integer(1),
      Fstatistic=numeric(1), Pvalue=numeric(1), Directory=character(1),
      stringsAsFactors = F)
    evnames <- roundmodel
    dvnames <- unlist(lapply(ev[evnames], names), use.names = FALSE)
    modeldir <- .dirpath.create(rounddir, "model1")
    df <- data.frame(ev[evnames])
    colnames(df) <- dvnames
    .runjar(rv, df, maxbkg = length(rv) + 1, modeldir)
    maxRes <- utils::read.csv(file.path(modeldir, "maxentResults.csv"))
    ctable$round[1] <- 0
    ctable$model[1] <- 1
    ctable$EV[1] <- paste(evnames, collapse = " ")
    ctable$m[1] <- length(dvnames)
    ctable$trainAUC[1] <- maxRes$Training.AUC
    ctable$Entropy[1] <- maxRes$Entropy
    n <- maxRes$X.Training.samples
    N <- maxRes$X.Background.points
    ctable$FVA[1] <- (log(N) - ctable$Entropy[1]) / (log(N) - log(n))
    ctable$addedFVA[1] <- ctable$FVA[1] - bestFVA
    ctable$dfe[1] <- ctable$m[1] - mnull
    ctable$dfu[1] <- (N - n) - (ctable$m[1] + 1) - 1
    ctable$Fstatistic[1] <- (ctable$addedFVA[1] * ctable$dfu[1]) /
      ((1 - ctable$FVA[1]) * ctable$dfe[1])
    ctable$Pvalue[1] <- 1 - stats::pf(ctable$Fstatistic[1], ctable$dfe[1],
      ctable$dfu[1])
    ctable$Directory[1] <- modeldir
    modeltable <- ctable
    selectedset <- roundmodel
    mnull <- ctable$m[1]
    bestFVA <- ctable$FVA[1]
    remainingset <- names(ev)[!(names(ev) %in% roundmodel)]
    message("Round 0 complete.")
  } else {
    selectedset <- character(length=0)
    remainingset <- names(ev)
    modeltable <- data.frame()
    roundnumber <- 0
    mnull <- 0
    bestFVA <- 0
  }

  iterationexit <- FALSE
  while (iterationexit == FALSE) {

    roundnumber <- roundnumber + 1
    rounddir <- .dirpath.create(dir, paste0("round", roundnumber))
    roundmodels <- lapply(remainingset, function(x) c(selectedset, x))

    nrows <- length(roundmodels)
    ctable <- data.frame(round=integer(nrows), model=integer(nrows),
      EV=character(nrows), m=integer(nrows), trainAUC=numeric(nrows),
      Entropy=numeric(nrows), FVA=numeric(nrows), addedFVA=numeric(nrows),
      dfe=integer(nrows), dfu=integer(nrows), Fstatistic=numeric(nrows),
      Pvalue=numeric(nrows), Directory=character(nrows),
      stringsAsFactors = F)

    for (i in 1:length(roundmodels)) {
      evnames <- roundmodels[[i]]
      dvnames <- unlist(lapply(ev[evnames], names), use.names = FALSE)
      modeldir <- .dirpath.create(rounddir, paste0("model", i))
      df <- data.frame(ev[evnames])
      colnames(df) <- dvnames
      .runjar(rv, df, maxbkg = length(rv) + 1, modeldir)

      maxRes <- utils::read.csv(file.path(modeldir, "maxentResults.csv"))
      ctable$round[i] <- roundnumber
      ctable$model[i] <- i
      ctable$EV[i] <- paste(evnames, collapse = " ")
      ctable$m[i] <- length(dvnames)
      ctable$trainAUC[i] <- maxRes$Training.AUC
      ctable$Entropy[i] <- maxRes$Entropy
      n <- maxRes$X.Training.samples
      N <- maxRes$X.Background.points
      ctable$FVA[i] <- (log(N) - ctable$Entropy[i]) / (log(N) - log(n))
      ctable$addedFVA[i] <- ctable$FVA[i] - bestFVA
      ctable$dfe[i] <- ctable$m[i] - mnull
      ctable$dfu[i] <- (N - n) - (ctable$m[i] + 1) - 1
      ctable$Fstatistic[i] <- (ctable$addedFVA[i] * ctable$dfu[i]) /
        ((1 - ctable$FVA[i]) * ctable$dfe[i])
      ctable$Pvalue[i] <- 1 - stats::pf(ctable$Fstatistic[i], ctable$dfe[i],
        ctable$dfu[i])
      ctable$Directory[i] <- modeldir
    }

    ctable <- ctable[order(ctable$Pvalue, -ctable$Fstatistic), ]
    modeltable <- rbind(modeltable, ctable, make.row.names = FALSE)

    if (ctable$Pvalue[1] < alpha) {
      selectedset <- unlist(strsplit(ctable$EV[1], split=" "))
      mnull <- ctable$m[1]
      bestFVA <- ctable$FVA[1]
      addedEV <- sapply(strsplit(ctable$EV[seq(nrow(ctable))[-1]], split=" "),
        function(x) {x[length(selectedset)]})
      remainingset <- addedEV[ctable$Pvalue[seq(nrow(ctable))[-1]] < alpha]
    }

    message(paste0("Round ", roundnumber, " complete."))

    if (nrow(ctable) == 1 || ctable$Pvalue[1] > alpha ||
        (ctable$Pvalue[1] < alpha &&
            all(ctable$Pvalue[seq(nrow(ctable))[-1]] >= alpha))) {
      iterationexit <- TRUE
    }
  }

  if (interaction == FALSE || length(selectedset) < 2) {
    return(list(ev[selectedset], modeltable))
  }

  message(paste0("Forward selection of interaction terms between ",
    length(selectedset), " EVs"))

  combos <- t(utils::combn(selectedset, 2))
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

    roundnumber <- roundnumber + 1
    rounddir <- .dirpath.create(dir, paste0("round", roundnumber))
    roundmodels <- lapply(remainingset, function(x) c(selectedset, x))

    nrows <- length(roundmodels)
    ctable <- data.frame(round=integer(nrows), model=integer(nrows),
      EV=character(nrows), m=integer(nrows), trainAUC=numeric(nrows),
      Entropy=numeric(nrows), FVA=numeric(nrows), addedFVA=numeric(nrows),
      dfe=integer(nrows), dfu=integer(nrows), Fstatistic=numeric(nrows),
      Pvalue=numeric(nrows), Directory=character(nrows),
      stringsAsFactors = F)

    for (i in 1:length(roundmodels)) {
      evnames <- roundmodels[[i]]
      dvnames <- unlist(lapply(ev[evnames], names), use.names = FALSE)
      modeldir <- .dirpath.create(rounddir, paste0("model", i))
      df <- data.frame(ev[evnames])
      colnames(df) <- dvnames
      .runjar(rv, df, maxbkg = length(rv) + 1, modeldir)

      maxRes <- utils::read.csv(file.path(modeldir, "maxentResults.csv"))
      ctable$round[i] <- roundnumber
      ctable$model[i] <- i
      ctable$EV[i] <- paste(evnames, collapse = " ")
      ctable$m[i] <- length(dvnames)
      ctable$trainAUC[i] <- maxRes$Training.AUC
      ctable$Entropy[i] <- maxRes$Entropy
      n <- maxRes$X.Training.samples
      N <- maxRes$X.Background.points
      ctable$FVA[i] <- (log(N) - ctable$Entropy[i]) / (log(N) - log(n))
      ctable$addedFVA[i] <- ctable$FVA[i] - bestFVA
      ctable$dfe[i] <- ctable$m[i] - mnull
      ctable$dfu[i] <- (N - n) - (ctable$m[i] + 1) - 1
      ctable$Fstatistic[i] <- (ctable$addedFVA[i] * ctable$dfu[i]) /
        ((1 - ctable$FVA[i]) * ctable$dfe[i])
      ctable$Pvalue[i] <- 1 - stats::pf(ctable$Fstatistic[i], ctable$dfe[i],
        ctable$dfu[i])
      ctable$Directory[i] <- modeldir
    }

    ctable <- ctable[order(ctable$Pvalue, -ctable$Fstatistic), ]
    modeltable <- rbind(modeltable, ctable, make.row.names = FALSE)

    if (ctable$Pvalue[1] < alpha) {
      selectedset <- unlist(strsplit(ctable$EV[1], split=" "))
      mnull <- ctable$m[1]
      bestFVA <- ctable$FVA[1]
      addedEV <- sapply(strsplit(ctable$EV[seq(nrow(ctable))[-1]], split=" "),
        function(x) {x[length(selectedset)]})
      remainingset <- addedEV[ctable$Pvalue[seq(nrow(ctable))[-1]] < alpha]
    }

    message(paste0("Round ", roundnumber, " complete."))

    if (nrow(ctable) == 1 || ctable$Pvalue[1] > alpha ||
        (ctable$Pvalue[1] < alpha &&
            all(ctable$Pvalue[seq(nrow(ctable))[-1]] >= alpha))) {
      iterationexit <- TRUE
    }
  }

  return(list(ev[selectedset], modeltable))

}