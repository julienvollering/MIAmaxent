#' selects a subset of spline dvs based on Maxent FTVE
#'
#' @param rv Vector of response variable values
#' @param dvs List of spline dvs to be selected from (HF, HR, or Th)
#'

.splselect <- function(rv, dvs) {

  n <- length(dvs)
  ctable <- data.frame(DV=character(n), KnotPosition=numeric(n), n=integer(n),
                       N=integer(n), Entropy=numeric(n), FVA=numeric(n),
                       dfu=integer(n), Fstatistic=numeric(n), Pvalue=numeric(n),
                       stringsAsFactors = F)

  pb <- utils::txtProgressBar(min = 0, max = n, style = 3)

  for (i in 1:n) {
    dvname <- names(dvs)[[i]]
    df <- data.frame(rv, dvs[[i]])
    colnames(df) <- c("RV", dvname)
    formula <- stats::formula(paste("RV ~", dvname))
    iwlr <- .runIWLR(formula, df)
    ctable$DV[i] <- dvname
    ctable$KnotPosition[i] <- (2 * i - 1) / (2 * n)
    ctable$n[i] <- sum(df[,"RV"]==1, na.rm=TRUE)
    ctable$N[i] <- nrow(df)
    ctable$Entropy[i] <- iwlr$entropy
    ctable$FVA[i] <- (log(ctable$N[i]) - ctable$Entropy[i]) /
                         (log(ctable$N[i]) - log(ctable$n[i]))
    ctable$dfu[i] <- (ctable$N[i] - ctable$n[i]) - 2 - 1
    ctable$Fstatistic[i] <- (ctable$FVA[i] * ctable$dfu[i]) /
                                ((1-ctable$FVA[i]) * 1)
    ctable$Pvalue[i] <- 1 - stats::pf(ctable$Fstatistic[i], 1, ctable$dfu[i])

    utils::setTxtProgressBar(pb, i)
  }


  selected <- character()
  for (i in 3:(nrow(ctable)-2)) {
    if (ctable$FVA[i] >= ctable$FVA[i-2] &&
        ctable$FVA[i] >= ctable$FVA[i-1] &&
        ctable$FVA[i] >= ctable$FVA[i+1] &&
        ctable$FVA[i] >= ctable$FVA[i+2] &&
        ctable$Pvalue[i] < 0.05) {
      selected <- append(selected, ctable$DV[i])
    }
  }

  return(selected)
}