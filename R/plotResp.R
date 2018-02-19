#' Plot single-effect model response.
#'
#' \code{plotResp} plots the single-effect response of a given model over any of
#' the explanatory variables (EVs) included in that model. For categorical
#' variables, a bar plot is returned rather than a scatter plot. Single-effect
#' response curves present the response of a model containing the explanatory
#' variable of interest only (cf. marginal-effect response curves;
#' \code{\link{plotResp2}}).
#'
#' @param model The model for which the response is to be plotted, represented
#'   by an object of class 'glm'. This may be the object returned by
#'   \code{\link{chooseModel}}, or the 'selectedmodel' returned by
#'   \code{\link{selectEV}}.
#' @param transformations Transformation functions used to create the derived
#'   variables in the model. I.e. the 'transformations' returned by
#'   \code{\link{deriveVars}}. Equivalently, the full file pathway of the
#'   'transformations.Rdata' file saved as a result of \code{\link{deriveVars}}.
#' @param EV Name of the explanatory variable for which the response curve is to
#'   be plotted. Interaction terms not allowed.
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


plotResp <- function(model, transformations, EV, logscale = FALSE, ...) {

  evbetas <- model$betas[grep(paste0(EV, "_"), names(model$betas))]
  evbetasni <- evbetas[!grepl(":", names(evbetas), fixed=TRUE)]
  if (length(evbetasni)==0) {
    stop("The 'EV' specified cannot be found in the model")
  }

  alltransf <- .load.transf(transformations)
  evtransfs <- alltransf[match(paste0(names(evbetasni), "_transf"),
                               names(alltransf), nomatch = 0)]
  if (!(length(evtransfs)==length(evbetasni))) {
    stop("The transformation function for at least one DV in the model is missing")
  }

  dvdata <- lapply(evtransfs, function(f) {
    x <- environment(f)$xnull
    f(x) })
  names(dvdata) <- names(evbetasni)
  traindata <- data.frame("RV"=alltransf[[1]], dvdata)
  formula <- stats::formula(paste("RV ~", paste(names(evbetasni), collapse = " + ")))
  smodel <- .runIWLR(formula, traindata)

  evnull <- environment(evtransfs[[1]])$xnull
  if (class(evnull) %in% c("numeric", "integer")) {
    seq <- seq(min(evnull), max(evnull), length.out = 100)
  }
  if (class(evnull) %in% c("factor", "character")) {
    seq <- levels(as.factor(evnull))
  }
  newdata <- as.data.frame(do.call(cbind,
                                   lapply(evtransfs, function(f, x) {
                                     f(x) }, x=seq)))
  names(newdata) <- names(evbetasni)
  mmformula <- stats::update.formula(smodel$formula.narm, NULL ~ . - 1)
  newdata <- model.matrix(mmformula, newdata)

  raw <- exp((newdata %*% smodel$betas) + smodel$alpha)

  resp <- data.frame(EV = seq, PRO = raw*nrow(data))
  if (logscale == TRUE) {resp$PRO <- log10(resp$PRO)}

  if (class(resp[, 1]) %in% c("numeric", "integer")) {
    graphics::plot(resp[, 2] ~ resp[, 1], type="l", col="red", ...,
      main = paste0("Single-effect response plot: ", EV), xlab = EV,
      ylab = ifelse(logscale == TRUE, "log Probability Ratio Output (logPRO)",
        "Probability Ratio Output (PRO)"))
  }

  if (class(resp[, 1]) %in% c("factor", "character")) {
    graphics::barplot(resp[, 2], names.arg = resp[, 1], col="red", ...,
      main = paste0("Single-effect response plot: ", EV), xlab = EV,
      ylab = ifelse(logscale == TRUE, "log Probability Ratio Output (logPRO)",
        "Probability Ratio Output (PRO)"))
  }

  if (logscale == TRUE) { graphics::abline(h = 0, lty = 3)
  } else { graphics::abline(h = 1, lty = 3) }

}
