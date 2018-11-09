#' Creates (standard) logistic regression model without regularization
#'
#' @param formula Object of class "formula": a symbolic description of the model
#'   to be fitted.
#' @param data Data frame containing the variables in the model. Response
#'   variable as binary (1/NA or 1/0).
#'
#' @keywords internal
#' @noRd

.runLR <- function(formula, data) {
  RV <- all.vars(formula)[1]
  data[,RV][is.na(data[,RV])] <- 0

  withCallingHandlers({
    model <- stats::glm(formula=formula, family="binomial", data=data)
  }, warning = function(w) {
    if(grepl("fitted probabilities numerically 0 or 1", conditionMessage(w)) ||
       grepl("glm.fit: algorithm did not converge", conditionMessage(w))){
      invokeRestart("muffleWarning")
    }
  })

  if (any(is.na(model$coefficients))) {
    nacoef <- names(model$coefficients)[is.na(model$coefficients)]
    model$formula.narm <- stats::update(model$formula,
                                        paste("~ . -", paste(nacoef, collapse = " - ")))
    model$betas <- model$coefficients[-1][!is.na(model$coefficients[-1])]
  } else {
    model$formula.narm <- model$formula
    model$betas <- model$coefficients[-1]
  }
  class(model) <- c("lr", class(model))
  return(model)
}