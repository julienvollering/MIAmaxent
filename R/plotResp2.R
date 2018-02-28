#' @describeIn plotResp Plot marginal-effect model response.
#'
#' @examples
#'
#' @export


plotResp2 <- function(model, transformations, EV, logscale = FALSE, ...) {

  if (!(class(model)[1] %in% c("iwlr", "lr"))) {
    stop("'model' should be of the class produced by 'selectEV' or 'chooseModel'", call. = FALSE)
  }
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
  type <- ifelse(class(model)[1] == "iwlr", "PRO", "response")
  preds <- stats::predict(model, newdata, type)
  resp <- data.frame(EV = seq, preds = preds)

  if (logscale == TRUE) { resp$preds <- log10(resp$preds) }
  ylab <- ifelse(type == "PRO", "Probability Ratio Output (PRO)", "Predicted probability")
  if (logscale == TRUE) { ylab <- paste("log", ylab) }
  args1 <- list(main = paste0("Marginal-effect response plot: ", EV), xlab = EV,
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
