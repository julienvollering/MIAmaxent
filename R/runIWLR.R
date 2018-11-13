#' Creates infinitely weighted logistic regression model (equivalent to maxent)
#' without regularization
#'
#' @param formula Object of class "formula": a symbolic description of the model
#'   to be fitted. Do not use '.' term, as weights are concatenated to the data
#'   object.
#' @param data Data frame containing the variables in the model. Response
#'   variable as (1/NA).
#'
#' @keywords internal
#' @noRd

.runIWLR <- function(formula, data) {
  RV <- all.vars(formula)[1]
  data[, RV][is.na(data[, RV])] <- 0
  padd <- data[data[, RV]==1, ]
  padd[, RV] <- 0
  padddata <- rbind(data, padd)
  # Code below this line was modified from the MIT-licensed 'maxnet' library
  # Copyright (c) 2016, Steven Phillips
  wgts <- padddata[, RV] + (1 - padddata[, RV])*100
  glmdata <- cbind(padddata, wgts)

  withCallingHandlers({
    model <- stats::glm(formula=formula, family="binomial", data=glmdata,
                        weights=wgts)
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
  bkg <- stats::model.matrix(model$formula.narm,
                             padddata[padddata[, RV]==0, ])[, -1, drop=FALSE]
  model$alpha <- 0
  link <- (bkg %*% model$betas) + model$alpha
  rr <- exp(link)
  raw <- rr / sum(rr)
  model$entropy <- -sum(raw * log(raw), na.rm = TRUE)
  model$alpha <- -log(sum(rr))
  class(model) <- c("iwlr", class(model))
  return(model)
  # Code above this line was modified from the MIT-licensed 'maxnet' library
  # Copyright (c) 2016, Steven Phillips
}



#' Predict method for infinitely-weighted logistic regression
#'
#' Returns model predictions for new data in "PRO" or "raw" format.
#'
#' @param object Model of class "iwlr"
#' @param newdata Data frame containing variables with which to predict
#' @param type Type of model output: "PRO" or "raw"
#'
#' @keywords internal
#'
#' @export
#'
#' @method predict iwlr

predict.iwlr <- function(object, newdata, type="PRO", ...) {
  mmformula <- stats::update.formula(object$formula.narm, NULL ~ . - 1)
  newdata <- stats::model.matrix(mmformula, newdata)
  raw <- exp((newdata %*% object$betas) + object$alpha)
  RV <- all.vars(object$formula)[1]
  N <- sum(object$data[,RV] == 0)
  PRO <- raw * N
  if (type == "PRO") {return(PRO)} else {return(raw)}
}