#' Compares a list of models to a reference model and returns a comparison table
#'
#' @param formulas List of formulas for models to be compared to the null model.
#' @param refformula Formula for the reference model.
#' @param data Data frame containing the variables referenced in all formulas.
#' @param test Character string matching either "Chisq" or "F". F-test is sensu
#'   Halvorsen (2013, 2015).
#' @param algorithm Character string matching either "maxent" or "LR".
#'
#' @keywords internal
#' @noRd

.compare <- function(formulas, refformula, data, test="Chisq", algorithm="maxent") {

  test <- match.arg(test, choices = c("Chisq", "F"))
  algorithm <- match.arg(algorithm, choices = c("maxent", "LR"))
  nformulas <- length(formulas)

  if (algorithm == "maxent") {
    refmod <- .runIWLR(refformula, data)
  } else {
    refmod <- .runLR(refformula, data)
  }

  if (test == "Chisq") {
    ctable <- data.frame(variables=character(nformulas), m=integer(nformulas),
                         Dsq=numeric(nformulas), Chisq=numeric(nformulas), df=integer(nformulas),
                         P=numeric(nformulas), stringsAsFactors = F)

    for (i in 1:length(formulas)) {
      if (algorithm == "maxent") {
        mod <- .runIWLR(formulas[[i]], data)
      } else {
        mod <- .runLR(formulas[[i]], data)
      }
      ctable$variables[i] <- paste(labels(stats::terms(formulas[[i]])),
                                   collapse = " + ")
      ctable$m[i] <- length(mod$coefficients)-1
      ctable$Dsq[i] <- round(1 - (mod$deviance/mod$null.deviance), digits = 3)
      a2 <- stats::anova(refmod, mod, test="Chisq")
      ctable$Chisq[i] <- round(a2$Deviance[2], digits = 3)
      ctable$df[i] <- a2$Df[2]
      ctable$P[i] <- signif(a2$`Pr(>Chi)`[2], digits = 3)
    }
  }

  if (test == "F") {
    ctable <- data.frame(variables=character(nformulas), m=integer(nformulas),
                         Dsq=numeric(nformulas), F=numeric(nformulas), dfe=integer(nformulas),
                         dfu=integer(nformulas), P=numeric(nformulas), stringsAsFactors = F)

    for (i in 1:length(formulas)) {
      if (algorithm == "maxent") {
        mod <- .runIWLR(formulas[[i]], data)
      } else {
        mod <- .runLR(formulas[[i]], data)
      }
      ctable$variables[i] <- paste(labels(stats::terms(formulas[[i]])),
                                   collapse = " + ")
      ctable$m[i] <- length(mod$coefficients)-1
      Dsq <- 1 - (mod$deviance/mod$null.deviance)
      ctable$Dsq[i] <- round(Dsq, digits = 3)
      N <- nrow(data)
      n <- sum(data[,1]==1, na.rm = TRUE)
      a2 <- stats::anova(refmod, mod)
      ctable$dfe[i] <- a2$Df[2]
      ctable$dfu[i] <- (N - n) - (ctable$m[i] + 1) - 1
      ctable$F[i] <- round((a2$Deviance[2] * ctable$dfu[i]) /
        (a2$`Resid. Dev`[2] * ctable$dfe[i]), digits = 3)
      ctable$P[i] <- signif(1 - stats::pf(ctable$F[i], ctable$dfe[i], ctable$dfu[i]),
                            digits = 3)
    }
  }

  return(ctable)
}