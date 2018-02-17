#' Creates infinitely weighted logistic regression model (equivalent to maxent)
#' without regularization
#'
#' @param formula Object of class "formula": a symbolic description of the model
#'   to be fitted. Do not use '.' term, as weights are concatenated to the data
#'   object.
#' @param data Data frame containing the variables in the model. Response
#'   variable as (1/NA).

.runIWLR <- function(formula, data) {
  RV <- all.vars(formula)[1]
  data[,RV][is.na(data[,RV])] <- 0
  padd <- data[data[,RV]==1, ]
  padd[, RV] <- 0
  padddata <- rbind(data, padd)
  # Code below this line was modified from the MIT-licensed 'maxnet' library
  wgts <- padddata[,1]+(1-padddata[,1])*100
  glmdata <- cbind(padddata, wgts)

  withCallingHandlers({
    model <- stats::glm(formula=formula, family=binomial, data=glmdata,
                        weights=wgts)
  }, warning = function(w) {
    if(grepl("fitted probabilities numerically 0 or 1", conditionMessage(w))){
      invokeRestart("muffleWarning")
    }
  })

  if (any(is.na(model$coefficients))) {
    nacoef <- names(model$coefficients)[is.na(model$coefficients)]
    model$betas <- model$coefficients[-1][!is.na(model$coefficients[-1])]
    formula <- stats::update(formula, paste("~ . -", nacoef))
  } else {
    model$betas <- model$coefficients[-1]
  }
  bkg <- model.matrix(formula, padddata[padddata[, RV]==0, ])[,-1, drop=FALSE]
  model$alpha <- 0
  link <- (bkg %*% model$betas) + model$alpha
  rr <- exp(link)
  raw <- rr / sum(rr)
  model$entropy <- -sum(raw * log(raw), na.rm = TRUE)
  model$alpha <- -log(sum(rr))
  return(model)
  # Code above this line was modified from the MIT-licensed 'maxnet' library
}