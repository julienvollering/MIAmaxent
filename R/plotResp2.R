#' Plot marginal-effect model response.
#'
#' \code{plotResp2} plots the marginal-effect response of a given model over any
#' of the explanatory variables (EVs) included in that model. For categorical
#' variables, a bar plot is returned rather than a scatter plot. Marginal-effect
#' response curves present the response of the model when all other explanatory
#' variables are held constant at their mean values (cf. single-effect response
#' curves; \code{\link{plotResp}}).
#'
#' @param model The model for which the response is to be plotted, represented
#'   by an object of class 'glm'. This may be the object returned by
#'   \code{\link{chooseModel}}, or the 'selectedmodel' returned by
#'   \code{\link{selectEV}}.
#' @param transformations Transformation functions used to create the derived
#'   variables in the model. I.e. the 'transformations' returned by
#'   \code{\link{deriveVars}}. Equivalently, the full file pathway of the
#'   'transformations.Rdata' file saved as a result of \code{\link{deriveVars}}.
#' @param EV Character. Name of the explanatory variable for which the response curve is to
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


plotResp2 <- function(model, transformations, EV, logscale = FALSE, ...) {

  evbetas <- model$betas[grep(paste0(EV, "_"), names(model$betas))]
  evbetasni <- evbetas[!grepl(":", names(evbetas), fixed=TRUE)]
  if (length(evbetasni)==0) {
    stop("The 'EV' specified cannot be found in the model")
  }

  betasni <- model$betas[!grepl(":", names(model$betas), fixed=TRUE)]
  alltransf <- .load.transf(transformations)
  modtransfs <- alltransf[match(paste0(names(betasni), "_transf"),
                               names(alltransf), nomatch = 0)]
  if (!(length(modtransfs)==length(betasni))) {
    stop("The transformation function for at least one DV in the model is missing")
  }

  anevtransf <- modtransfs[[paste0(names(evbetasni)[1], "_transf")]]
  evnull <- environment(anevtransf)$xnull
  if (class(evnull) %in% c("numeric", "integer")) {
    seq <- seq(min(evnull), max(evnull), length.out = 100)
  }
  if (class(evnull) %in% c("factor", "character")) {
    seq <- levels(as.factor(evnull))
  }

  marginal <- grepl(paste0(EV, "_"), names(betasni))

  newdata <- mapply(function(f, marginal, seq) {
    if (marginal) {
      dv <- f(seq)
    } else {
      x <- environment(f)$xnull
      if (class(x) %in% c("numeric", "integer")) {
        dvmean <- f(mean(x, na.rm = TRUE))
        dv <- rep(dvmean, length(seq))
      }
      if (class(x) %in% c("factor", "character")) {
        ux <- unique(x)
        mode <- ux[which.max(tabulate(match(x, ux)))]
        dvmode <- f(mode)
        dv <- rep(dvmode, length(seq))
      }
    }
   return(dv)
  }, modtransfs, marginal, MoreArgs = list("seq"=seq))
  colnames(newdata) <- names(betasni)
  newdata <- as.data.frame(newdata)
  mmformula <- stats::update.formula(model$formula.narm, NULL ~ . - 1)
  newdata <- model.matrix(mmformula, newdata)

  raw <- exp((newdata %*% model$betas) + model$alpha)
  resp <- data.frame(EV = seq, PRO = raw*length(transformations[[1]]))
  if (logscale == TRUE) {resp$PRO <- log10(resp$PRO)}

  args1 <- list(main = paste0("Single-effect response plot: ", EV), xlab = EV,
                ylab = ifelse(logscale == TRUE,
                              "log Probability Ratio Output (logPRO)",
                              "Probability Ratio Output (PRO)"), col="red")
  inargs <- list(...)
  args1[names(inargs)] <- inargs

  if (class(resp[, 1]) %in% c("numeric", "integer")) {
    do.call(graphics::plot, c(list(x=resp[, 1], y=resp[, 2], type="l"), args1))
  }

  if (class(resp[, 1]) %in% c("factor", "character")) {
    do.call(graphics::barplot, c(list(height=resp[, 2], names.arg=resp[, 1]), args1))
  }

  if (logscale == TRUE) { graphics::abline(h = 0, lty = 3)
  } else { graphics::abline(h = 1, lty = 3) }

}
