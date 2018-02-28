#' Plot model response.
#'
#' Plots the response of a given model over any of the explanatory variables
#' (EVs) included in that model. For categorical variables, a bar plot is
#' returned rather than a line plot. Single-effect response curves show the
#' response of a model containing the explanatory variable of interest only,
#' while marginal effect response curves show the response of the model when all
#' other explanatory variables are held constant at their mean values (cf.
#' \code{plotResp}, \code{plotResp2}).
#'
#' @param model The model for which the response is to be plotted, represented
#'   by an object of class 'glm'. This may be the object returned by
#'   \code{\link{chooseModel}}, or the 'selectedmodel' returned by
#'   \code{\link{selectEV}}.
#' @param transformations Transformation functions used to create the derived
#'   variables in the model. I.e. the 'transformations' returned by
#'   \code{\link{deriveVars}}. Equivalently, the full file pathway of the
#'   'transformations.Rdata' file saved as a result of \code{\link{deriveVars}}.
#' @param EV Character. Name of the explanatory variable for which the response
#'   curve is to be plotted. Interaction terms not allowed.
#' @param logscale Logical. Plot the common logarithm of PRO rather than PRO
#'   itself.
#' @param ... Arguments to be passed to \code{plot} or \code{barplot} to control
#'   the appearance of the plot. For example: \itemize{ \item \code{lwd} for
#'   line width \item \code{cex.main} for size of plot title \item \code{space}
#'   for space between bars }
#'
#' @describeIn plotResp Plot single-effect model response.
#'
#' @examples
#' \dontrun{
#' # From vignette:
#' plotResp(grasslandmodel, grasslandDVs$transformations, "pr.bygall")
#' plotResp(grasslandmodel, grasslandDVs$transformations, "geolmja1")
#'
#' plotResp2(grasslandmodel, grasslandDVs$transformations, "pr.bygall")
#' }
#'
#' @export


plotResp <- function(model, transformations, EV, logscale = FALSE, ...) {

  if (!(class(model)[1] %in% c("iwlr", "lr"))) {
    stop("'model' should be of the class produced by 'selectEV' or 'chooseModel'", call. = FALSE)
  }
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

  if (class(model)[1] == "iwlr") {
    smodel <- .runIWLR(formula, traindata)
  } else if (class(model)[1] == "lr") {
    smodel <- .runLR(formula, traindata)
  }

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
  type <- ifelse(class(model)[1] == "iwlr", "PRO", "response")
  preds <- stats::predict(smodel, newdata, type)
  resp <- data.frame(EV = seq, preds = preds)

  if (logscale == TRUE) { resp$preds <- log10(resp$preds) }
  ylab <- ifelse(type == "PRO", "Probability Ratio Output (PRO)", "Predicted probability")
  if (logscale == TRUE) { ylab <- paste("log", ylab) }
  args1 <- list(main = paste0("Single-effect response plot: ", EV), xlab = EV,
                ylab = ylab, col="red")
  inargs <- list(...)
  args1[names(inargs)] <- inargs

  if (class(resp[, 1]) %in% c("numeric", "integer")) {
    do.call(graphics::plot, c(list(x=resp[, 1], y=resp[, 2], type="l"), args1))
  }

  if (class(resp[, 1]) %in% c("factor", "character")) {
    do.call(graphics::barplot, c(list(height=resp[, 2], names.arg=resp[, 1]), args1))
  }

  if (type == "PRO") {
    if (logscale == TRUE) { graphics::abline(h = 0, lty = 3)
    } else { graphics::abline(h = 1, lty = 3) }
  }

}
