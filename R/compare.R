#' Compares a list of models to a reference model and returns a comparison table
#'
#' @param formulas List of formulas for models to be compared to the null model.
#' @param refformula Formula for the reference model.
#' @param data Data frame containing the variables referenced in all formulas.
#' @param test Character string matching either "Chisq" or "F". F-test is sensu
#'   Halvorsen (2013, 2015).

.compare <- function(formulas, refformula, data, test="Chisq") {
  test <- match.arg(test, choices = c("Chisq", "F"))
  n <- length(formulas)

  if (test == "Chisq") {
    ctable <- data.frame(variables=character(n), m=integer(n),
                         DSquared=numeric(n), ChiSquared=numeric(n),
                         df=integer(n), P=numeric(n), stringsAsFactors = F)

    refiwlr <- .runIWLR(refformula, data)

    for (i in 1:length(formulas)) {
      iwlr <- .runIWLR(formulas[[i]], data)
      ctable$variables[i] <- paste(labels(terms(formulas[[i]])),
                                   collapse = " + ")
      ctable$m[i] <- length(iwlr$coefficients)-1
      a1 <- stats::anova(iwlr, test="Chisq")
      ctable$DSquared[i] <- round(sum(a1$Deviance, na.rm = TRUE) /
                                    a1$`Resid. Dev`[1], digits = 3)
      a2 <- stats::anova(refiwlr, iwlr, test="Chisq")
      ctable$ChiSquared[i] <- round(a2$Deviance[2], digits = 3)
      ctable$df[i] <- a2$Df[2]
      ctable$P[i] <- signif(a2$`Pr(>Chi)`[2], digits = 3)
    }
  }

  if (test == "F") {
    ctable <- data.frame(variables=character(n), m=integer(n),
                         DSquared=numeric(n), F=numeric(n), dfe=integer(n),
                         dfu=integer(n), P=numeric(n), stringsAsFactors = F)

    refiwlr <- .runIWLR(refformula, data)

    for (i in 1:length(formulas)) {
      iwlr <- .runIWLR(formulas[[i]], data)
      ctable$variables[i] <- paste(labels(terms(formulas[[i]])),
                                   collapse = " + ")
      ctable$m[i] <- length(iwlr$coefficients)-1
      a1 <- stats::anova(iwlr)
      DSquared <- sum(a1$Deviance, na.rm = TRUE) / a1$`Resid. Dev`[1]
      ctable$DSquared[i] <- round(DSquared, digits = 3)
      N <- nrow(data)
      n <- sum(data[,1]==1, na.rm = TRUE)
      a2 <- stats::anova(refiwlr, iwlr)
      addedDSquared <- a2$Deviance[2] / a2$`Resid. Dev`[1]
      ctable$dfe[i] <- a2$Df[2]
      ctable$dfu[i] <- (N - n) - (ctable$m[i] + 1) - 1
      ctable$F[i] <- (addedDSquared * ctable$dfu[i]) /
        ((1 - DSquared) * ctable$dfe[i])
      ctable$P[i] <- 1 - stats::pf(ctable$F[i], ctable$dfe[i], ctable$dfu[i])
    }
  }

  return(ctable)
}