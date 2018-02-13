#' Forward selection of EVs
#'
#' @param data List of data frames with response variable (first element) and EVs to
#'   be selected from (remaining elements).
#' @param alpha Alpha level for F-test.
#' @param test Character string matching either "Chisq" or "F".
#' @param interaction Logical. Allows interaction terms.
#' @param formula Model formula specifying a starting point for model selection.
#'

.parsevs <- function(data, alpha, interaction, dir, formula) {

  dfnames <-  unlist(lapply(data, names))
  df <- data.frame(do.call(cbind, data))
  names(df) <- dfnames
  roundnumber <- 0

  if (!is.null(formula) && length(labels(stats::terms(formula))) != 0) {
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
    remainingset <- names(data)[-1]
    modeltable <- data.frame()
    refformula <- stats::formula(paste(names(df)[1], "~ 1"))
  }

  iterationexit <- FALSE
  while (iterationexit == FALSE) {

    roundnumber <- roundnumber + 1
    remainingsetdvs <- lapply(data[remainingset], names)
    formulas <- lapply(remainingsetdvs, function(x) {
      stats::update.formula(refformula,
                            paste("~ . +", paste(x, collapse = " + ")))})
    ctable <- .compare(formulas, refformula, df, test=test)
    if (length(selectedset)==0) {
      variables <- remainingset
    } else { variables <- paste(paste(selectedset, collapse=" + "),
                                remainingset, sep=" + ")}
    ctable$variables <- variables
    ctable <- ctable[order(ctable$P,
                           -ctable[, match(test, names(ctable))]), ]
    modeltable <- rbind(modeltable, data.frame("round"=roundnumber, ctable),
                     make.row.names = FALSE)

    if (!is.na(ctable$P[1]) && ctable$P[1] < alpha) {
      selectedset <- unlist(strsplit(ctable$variables[1], split=" + ",
                                     fixed=TRUE))
      selectedsetdvs <- unlist(lapply(data[selectedset], names))
      refformula <- stats::formula(paste(names(df)[1], "~",
                                         paste(selectedsetdvs, collapse=" + ")))
      testedEVs <- sapply(strsplit(ctable$variables[seq(nrow(ctable))[-1]],
                                   split=" + ", fixed=TRUE),
                          function(x) {x[length(selectedset)]})
      remainingset <- testedEVs[ctable$P[seq(nrow(ctable))[-1]] < alpha]
      remainingset <- remainingset[!is.na(remainingset)]
    }

    message(paste("Round", roundnumber, "complete."))

    if (nrow(ctable) == 1 ||
        ctable$P[1] > alpha ||
        (all(ctable$P[seq(nrow(ctable))[-1]] >= alpha |
             is.na(ctable$P[seq(nrow(ctable))[-1]])))) {
      iterationexit <- TRUE
    }
  }

  if (interaction == FALSE || length(selectedset) < 2) {
    return(list(data[selectedset], modeltable))
  }

  message(paste("Forward selection of interaction terms between",
    length(selectedset), "EVs"))

  # combos <- t(utils::combn(selectedset, 2))
  # products <- vector("list", nrow(combos))
  # names(products) <- apply(combos, 1, function(x) {paste(x, collapse=":")})
  # for (i in 1:nrow(combos)) {
  #   ev1 <- data[[combos[i,1]]]
  #   ev2 <- data[[combos[i,2]]]
  #   dvcombos <- expand.grid(names(ev1), names(ev2))
  #   productdvs <- matrix(nrow = nrow(data[[1]]), ncol = nrow(dvcombos))
  #   for (j in 1:nrow(dvcombos)) {
  #     productdvs[, j] <- ev1[[dvcombos[j,1]]] * ev2[[dvcombos[j,2]]]
  #   }
  #   colnames(productdvs) <- apply(dvcombos, 1, function(x) {
  #     paste(x, collapse=":")})
  #   products[[i]] <- as.data.frame(productdvs)
  # }
  #
  # data <- c(data, products)
  # dfnames <-  unlist(lapply(data, names))
  # df <- data.frame(do.call(cbind, data))
  # names(df) <- dfnames
  # remainingset <- names(products)


  remainingset <- apply(utils::combn(selectedset, 2), 2,
                        function(x) {paste(x, collapse=":")})

  iterationexit <- FALSE
  while (iterationexit == FALSE) {

    roundnumber <- roundnumber + 1
    remainingsetdvs <- lapply(strsplit(remainingset, ":", fixed=TRUE),
                              function(x, data) {
                                dvs1 <- names(data[[x[1]]])
                                dvs2 <- names(data[[x[2]]])
                                apply(expand.grid(dvs1, dvs2), 1, function(y) {
                                  paste(y, collapse=":")})
                              }, data=data)
    formulas <- lapply(remainingsetdvs, function(x) {
      stats::update.formula(refformula,
                            paste("~ . +", paste(x, collapse = " + ")))})
    ctable <- .compare(formulas, refformula, df, test=test)
    if (length(selectedset)==0) {
      variables <- remainingset
    } else { variables <- paste(paste(selectedset, collapse=" + "),
                                remainingset, sep=" + ")}
    ctable$variables <- variables
    ctable <- ctable[order(ctable$P,
                           -ctable[, match(test, names(ctable))]), ]
    modeltable <- rbind(modeltable, data.frame("round"=roundnumber, ctable),
                        make.row.names = FALSE)

    if (!is.na(ctable$P[1]) && ctable$P[1] < alpha) {
      selectedset <- unlist(strsplit(ctable$variables[1], split=" + ",
                                     fixed=TRUE))

      #Bookmark
      selectedsetdvs <- unlist(lapply(data[selectedset], names))
      refformula <- stats::formula(paste(names(df)[1], "~",
                                         paste(selectedsetdvs, collapse=" + ")))
      testedEVs <- sapply(strsplit(ctable$variables[seq(nrow(ctable))[-1]],
                                   split=" + ", fixed=TRUE),
                          function(x) {x[length(selectedset)]})
      remainingset <- testedEVs[ctable$P[seq(nrow(ctable))[-1]] < alpha]
      remainingset <- remainingset[!is.na(remainingset)]
    }

    message(paste("Round", roundnumber, "complete."))

    if (nrow(ctable) == 1 ||
        ctable$P[1] > alpha ||
        (all(ctable$P[seq(nrow(ctable))[-1]] >= alpha |
             is.na(ctable$P[seq(nrow(ctable))[-1]])))) {
      iterationexit <- TRUE
    }
  }

  return(list(ev[selectedset], modeltable))

}