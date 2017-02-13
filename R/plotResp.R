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
#' The \code{EV} specified in \code{dvdata} must not be an interaction term.
#'
#' @param data Data frame containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent presence/background data, coded as: 1/NA. See
#'   \code{\link{readData}}.
#' @param dvdata List of explanatory variables used to train the model, where
#'   each list item is a data frame containing \emph{selected} DVs for a
#'   \emph{selected} EV (e.g. the first item in the list returned by
#'   \code{\link{selectEV}}).
#' @param EV Name or list index of the explanatory variable in \code{dvdata} for
#'   which the response curve is to be generated. Interaction terms not allowed.
#' @param dir Directory to which files will be written. Defaults to the working
#'   directory.
#' @param logscale Logical. Plot the common logarithm of PRO rather than PRO
#'   itself.
#' @param ... Arguments to be passed to \code{plot} to control the appearance of
#'   the plot. For example: \itemize{ \item \code{cex} for size of points \item
#'   \code{col} for color \item \code{pch} for type }
#'
#' @return In addition to the graphical output, the plotted data is returned. In
#'   the case of a continuous EV, the plotted data is a list of 2: \enumerate{
#'   \item \code{respPts}. Model response across individual data points. Columns
#'   in this data frame represent the following: EV value ("EV"), Probability
#'   Ratio Output of the model ("PRO"), and corresponding EV interval ("int").
#'   \item \code{respLine}. Model response across intervals of the EV. Columns
#'   in this data frame represent the following: EV interval ("int"), number of
#'   points in the interval ("n"), mean EV value of the points in the interval
#'   ("intEV"), mean Probability Ratio Output of the points in the interval
#'   ("intPRO"), and exponentially weighted moving average of intPRO
#'   ("smoothPRO").}
#'
#'   In the case of a categorical EV, the plotted data is a data frame
#'   containing the number of points in the level ("n"), the level name
#'   ("level"), and the mean Probability Ratio Output of the level ("levelRV").
#'
#' @examples
#' \dontrun{
#' responseEV1 <- plotResp(dat, deriveddat, "EV1", dir = "D:/path/to/modeling/directory")
#'
#' # From vignette
#' pr_bygallResp <- plotResp(grasslandPO, grasslandEVselect[[1]], "pr_bygall")
#' par(mfrow = c(1,2))
#' pr_bygallFOP <- plotFOP(grasslandPO, "pr_bygall", intervals=50)
#' pr_bygallResp <- plotResp(grasslandPO, grasslandEVselect[[1]], "pr_bygall")
#' par(mfrow = c(1,1))
#' }
#'
#' @export


plotResp <- function(data, dvdata, EV, dir = NULL, logscale = FALSE, ...) {

  .binaryrvcheck(data[, 1])

  if (is.null(dir)) { dir <- getwd()}

  fdir <- file.path(dir, "plotResp")
  dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  evname <- names(dvdata[EV])
  if (!(evname %in% colnames(data))) {
    stop("The specified EV must be present in the untransformed data.
       E.g. interaction terms between multiple EVs are not supported. \n ",
      call. = FALSE)
  }
  modeldir <- file.path(fdir, paste0("response", evname))
  dir.create(modeldir, showWarnings = FALSE)

  .runjar(data[, 1], dvdata[[EV]], maxbkg = nrow(data) + 1, modeldir)

  output <- utils::read.csv(file.path(modeldir, "1_backgroundPredictions.csv"))
  respPts <- data.frame(EV = data[, evname], PRO = output[,3]*length(output[,3]))
  if (logscale == TRUE) {respPts$PRO <- log10(respPts$PRO)}

  if (class(respPts[, 1]) %in% c("numeric", "integer")) {
    graphics::plot(respPts[, 2] ~ respPts[, 1], ...,
      main = paste0("Single-effect response plot: ", evname), xlab = evname,
      ylab = ifelse(logscale == TRUE, "log Probability Ratio Output (logPRO)",
        "Probability Ratio Output (PRO)"))

    if (nrow(respPts) >= 250) {
      intervals <- min(c(ceiling(nrow(respPts) / 50), 100))
      respPts$int <- cut(respPts[, 1], max(2, intervals))
      grouped <- dplyr::group_by(respPts, int)
      respLine <- as.data.frame(dplyr::summarise(grouped, n = n(),
        intEV = mean(EV, na.rm = TRUE),
        intPRO = mean(PRO, na.rm = TRUE)))
      respLine$smoothPRO <- .ewma(respLine$intPRO, 5)
      graphics::lines(respLine$smoothPRO ~ respLine$intEV, col="red", lwd = 2)
      result <- list("respPts" = respPts, "respLine" = respLine)
    } else {
      result <- respPts
      warning("Exponentially weighted moving average line not drawn due to few intervals",
        call. = FALSE)
    }

    if (logscale == TRUE) {
      graphics::abline(h = 0, lty = 3)
    } else {
      graphics::abline(h = 1, lty = 3)
    }
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
