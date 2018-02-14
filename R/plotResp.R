#' Plot single-effect model response.
#'
#' \code{plotResp} plots the single-effect response of a given Maxent model over
#' any of the included explanatory variables (EVs) in that model. For
#' categorical variables, a bar plot is returned rather than a scatter plot.
#' \code{plotResp} also returns a data frame containing the plotted data (for
#' customizable graphics). Single-effect response curves present the response of
#' a model containing the explanatory variable of interest only (cf.
#' marginal-effect response curves; \code{\link{plotResp2}}).
#'
#' The plot contains points, representing the model response across individual
#' data points, as well as a line, representing an exponentially weighted moving
#' average of the model response over intervals of the EV.
#'
#' Note that the rows in \code{data} and \code{dvdata} must correspond to the
#' same observations, so that the transformed data relate correctly to their
#' untransformed values. This will automatically be the case if
#' \code{\link{deriveVars}} and \code{\link{selectDVforEV}} are used to create
#' \code{dvdata} from \code{data}.
#'
#' @param data Data frame containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent presence/background data, coded as: 1/NA. See
#'   \code{\link{readData}}.
#' @param dvdata List of explanatory variables used to train the model, where
#'   each list item is a data frame containing \emph{selected} DVs for a
#'   \emph{selected} EV (e.g. the first item in the list returned by
#'   \code{\link{selectEV}}).
#' @param EV Name of the explanatory variable or its list index in \code{dvdata}
#'   for which the response curve is to be generated. Interaction terms not
#'   allowed.
#' @param logscale Logical. Plot the common logarithm of PRO rather than PRO
#'   itself.
#' @param ... Arguments to be passed to \code{plot} to control the appearance of
#'   the plot. For example: \itemize{ \item \code{cex} for size of points \item
#'   \code{col} for color \item \code{pch} for type }
#'
#' @return In addition to the graphical output, the plotted data is returned. In
#'   the case of a continuous EV, the plotted data consists of: EV values ("EV")
#'   and Probability Ratio Output of the model ("PRO").
#'
#'   In the case of a categorical EV, the plotted data is a data frame
#'   containing the number of points in the level ("n"), the level name
#'   ("level"), and the mean Probability Ratio Output of the level ("levelPRO").
#'
#' @examples
#'
#' @export


plotResp <- function(data, dvdata, EV, logscale = FALSE, ...) {

  .binaryrvcheck(data[, 1])

  evname <- names(dvdata[EV])
  if (!(evname %in% colnames(data) && evname %in% names(dvdata))) {
    stop("The specified EV must be present in the both untransformed 'data' and the transformed 'dvdata'.
       E.g. interaction terms between multiple EVs are not supported. \n ",
      call. = FALSE)
  }

  df <- data.frame("RV"=data[,1], dvdata[[evname]])
  formula <- stats::formula(paste("RV ~", paste(names(df)[-1], collapse=" + ")))
  model <- .runIWLR(formula, df)
  df[,1][is.na(df[,1])] <- 0
  newdata <- model.matrix(stats::update(formula, ~. -1), df)
  raw <- exp((newdata %*% model$betas) + model$alpha)

  respPts <- data.frame(EV = data[, evname], PRO = raw*nrow(df))
  if (logscale == TRUE) {respPts$PRO <- log10(respPts$PRO)}

  if (class(respPts[, 1]) %in% c("numeric", "integer")) {
    graphics::plot(respPts[, 2] ~ respPts[, 1], ...,
      main = paste0("Single-effect response plot: ", evname), xlab = evname,
      ylab = ifelse(logscale == TRUE, "log Probability Ratio Output (logPRO)",
        "Probability Ratio Output (PRO)"))

    if (logscale == TRUE) { graphics::abline(h = 0, lty = 3) }
    else { graphics::abline(h = 1, lty = 3) }

    result <- respPts
  }

  if (class(respPts[, 1]) %in% c("factor", "character")) {
    respBar <- as.data.frame(dplyr::summarise(dplyr::group_by(respPts, EV),
      n = n(), levelPRO = mean(PRO, na.rm = TRUE)))
    graphics::barplot(respBar[, 3], names.arg = respBar[, 1], ...,
      main = paste0("Single-effect response plot: ", evname), xlab = evname,
      ylab = ifelse(logscale == TRUE, "log Probability Ratio Output (logPRO)",
        "Probability Ratio Output (PRO)"))
    if (logscale == TRUE) {
      graphics::abline(h = 0, lty = 3)
    } else {
      graphics::abline(h = 1, lty = 3)
    }
    result <- data.frame(n = respBar[, 2], level = respBar[, 1],
      levelPRO = respBar[, 3])
  }

  return(result)
}
