#' Forward selection of DVs
#'
#' @param df Data frame with response variable in first column and DVs to be
#'   selected from in subsequent columns.
#' @param alpha Alpha level for F-test.
#'

.parsdvs <- function(df, alpha) {

  selectedset <- character(length=0)
  remainingset <- names(df)[-1]
  evtable <- data.frame()
  roundnumber <- 0
  bestFVA <- 0

  iterationexit <- FALSE
  while (iterationexit == FALSE) {

    roundnumber <- roundnumber + 1
    roundmodels <- lapply(remainingset, function(x) c(selectedset, x))

    nrows <- length(roundmodels)
    ctable <- data.frame(round=integer(nrows), DV=character(nrows),
                         m=integer(nrows), Entropy=numeric(nrows),
                         FVA=numeric(nrows), addedFVA=numeric(nrows),
                         dfe=integer(nrows), dfu=integer(nrows),
                         Fstatistic=numeric(nrows), Pvalue=numeric(nrows),
                         devianceF=numeric(nrows), devianceP=numeric(nrows),
                         stringsAsFactors = F)

    for (i in 1:length(roundmodels)) {
      dvnames <- roundmodels[[i]]
      formula <- stats::formula(paste(paste(colnames(df)[1], "~"),
                                      paste0("`", dvnames, "`",
                                             collapse = " + ")))
      iwlr <- .runIWLR(formula, df)
      ctable$round[i] <- roundnumber
      ctable$DV[i] <- paste(dvnames, collapse = " ")
      ctable$m[i] <- length(dvnames)
      ctable$Entropy[i] <- iwlr$entropy
      n <- sum(df[, 1]==1, na.rm=TRUE)
      N <- nrow(df)
      ctable$FVA[i] <- (log(N) - ctable$Entropy[i]) / (log(N) - log(n))
      ctable$addedFVA[i] <- ctable$FVA[i] - bestFVA
      ctable$dfe[i] <- 1
      ctable$dfu[i] <- (N - n) - (ctable$m[i] + 1) - 1
      ctable$Fstatistic[i] <- (ctable$addedFVA[i] * ctable$dfu[i]) /
        ((1 - ctable$FVA[i]) * ctable$dfe[i])
      ctable$Pvalue[i] <- 1 - stats::pf(ctable$Fstatistic[i], ctable$dfe[i],
        ctable$dfu[i])
      ctable$devianceF[i] <- anova(iwlr, test="F")$F[2]
      ctable$devianceP[i] <- anova(iwlr, test="F")$Pr[2]
    }

    ctable <- ctable[order(ctable$Pvalue, -ctable$Fstatistic), ]
    evtable <- rbind(evtable, ctable, make.row.names = FALSE)

    if (ctable$Pvalue[1] < alpha) {
      selectedset <- unlist(strsplit(ctable$DV[1], split=" "))
      bestFVA <- ctable$FVA[1]
      addedDV <- sapply(strsplit(ctable$DV[seq(nrow(ctable))[-1]], split=" "),
        function(x) {x[length(selectedset)]})
      remainingset <- addedDV[ctable$Pvalue[seq(nrow(ctable))[-1]] < alpha]
    }

    if (nrow(ctable) == 1 || ctable$Pvalue[1] > alpha ||
        (ctable$Pvalue[1] < alpha &&
            all(ctable$Pvalue[seq(nrow(ctable))[-1]] >= alpha))) {
      iterationexit <- TRUE
    }
  }

  return(list(df[,selectedset, drop = FALSE], evtable))
}