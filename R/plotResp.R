#' Plot single-effect model response.
#'
#' \code{plotResp} plots the single-effect response of a given model over any of
#' the explanatory variables (EVs) included in that model. For categorical
#' variables, a bar plot is returned rather than a scatter plot. Single-effect
#' response curves present the response of a model containing the explanatory
#' variable of interest only (cf. marginal-effect response curves;
#' \code{\link{plotResp2}}).
#'
#' @param data Data frame containing the observations used to train the
#'   \code{model}, with the response variable in the first column and
#'   explanatory variables in subsequent columns. The response variable should
#'   represent presence/background data, coded as: 1/NA. See
#'   \code{\link{readData}}.
#' @param transformations Transformation functions used to create the derived
#'   variables in the model. I.e. the 'transformations' returned by
#'   \code{\link{deriveVars}}. Equivalently, the full file pathway of the
#'   'transformations.Rdata' file saved as a result of \code{\link{deriveVars}}.
#' @param model The model for which the response is to be plotted, represented
#'   by an object of class 'glm'. This may be the object returned by
#'   \code{\link{chooseModel}}, or the 'selectedmodel' returned by
#'   \code{\link{selectEV}}.
#' @param EV Name of the explanatory variable or its column index in \code{data}
#'   for which the response curve is to be plotted. Interaction terms not
#'   allowed.
#' @param logscale Logical. Plot the common logarithm of PRO rather than PRO
#'   itself.
#' @param ... Arguments to be passed to \code{plot} or \code{barplot} to control
#'   the appearance of the plot. For example: \itemize{ \item \code{lwd} for
#'   line width \item \code{cex.main} for size of plot title \item \code{space}
#'   for space between bars }
#'
#' @examples
#'
#' @export


plotResp <- function(data, transformations, model, EV, logscale = FALSE, ...) {

  .binaryrvcheck(data[, 1])

  if (EV==1) {
    stop("'EV' cannot be the first column of 'data', which must be the response variable")
  }

  evname <- names(data[EV])
  betas <- model$betas[grep(paste0(evname, "_"), names(model$betas))]
  if (length(betas)==0) {
    stop("The 'EV' specified is not in the model")
  }

  alltransf <- .load.transf(transformations)
  evtransfs <- alltransf[match(paste0(names(betas), "_transf"),
                               names(alltransf), nomatch = 0)]
  # if (!(length(evtransfs)==length(betas))) {
  #   stop("The transformation function for at least one DV in the model is missing")
  # }

  dvdata <- lapply(evtransfs, function(f, x) { f(x) }, x=data[,evname])
  names(dvdata) <- names(betas)
  traindata <- data.frame("RV"=data[,1], dvdata)
  formula <- stats::formula(paste("RV ~", paste(names(betas), collapse = " + ")))
  smodel <- .runIWLR(formula, traindata)

  if (class(data[, evname]) %in% c("numeric", "integer")) {
    seq <- seq(min(data[,evname]), max(data[,evname]), length.out = 100)
  }
  if (class(data[, evname]) %in% c("factor", "character")) {
    seq <- levels(as.factor(data[, evname]))
  }
  newdata <- do.call(cbind, lapply(evtransfs, function(f, x) { f(x) }, x=seq))
  colnames(newdata) <- names(betas)
  raw <- exp((newdata %*% smodel$betas) + smodel$alpha)

  resp <- data.frame(EV = seq, PRO = raw*nrow(data))
  if (logscale == TRUE) {resp$PRO <- log10(resp$PRO)}

  if (class(resp[, 1]) %in% c("numeric", "integer")) {
    graphics::plot(resp[, 2] ~ resp[, 1], type="l", col="red", ...,
      main = paste0("Single-effect response plot: ", evname), xlab = evname,
      ylab = ifelse(logscale == TRUE, "log Probability Ratio Output (logPRO)",
        "Probability Ratio Output (PRO)"))
  }

  if (class(resp[, 1]) %in% c("factor", "character")) {
    graphics::barplot(resp[, 2], names.arg = resp[, 1], col="red", ...,
      main = paste0("Single-effect response plot: ", evname), xlab = evname,
      ylab = ifelse(logscale == TRUE, "log Probability Ratio Output (logPRO)",
        "Probability Ratio Output (PRO)"))
  }

  if (logscale == TRUE) { graphics::abline(h = 0, lty = 3)
  } else { graphics::abline(h = 1, lty = 3) }

}
