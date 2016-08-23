#' Plot marginal-effect model response.
#'
#' \code{plotResp2} plots the marginal-effect response of a given Maxent model
#' over any of the included explanatory variables (EVs) in that model. For
#' categorical variables, a bar plot is returned rather than a scatter plot.
#' \code{plotResp2} also returns a data frame containing the plotted data (for
#' customizable graphics). Marginal-effect response curves present the response
#' of the model when all other explanatory variables are held constant at their
#' mean values (cf. single-effect response curves; \code{\link{plotResp}}).
#'
#' The plot contains points, representing the model response across individual
#' data points, as well as a line, representing an exponentially weighted moving
#' average of the model response over intervals of the EV.
#'
#' Model response is commonly plotted across EV values of the training data, but
#' it is possible to plot the model response over any EV values supplied in
#' \code{data}.
#'
#' The \code{EV} specified in \code{data} must not be an interaction term.
#'
#' @param data Data frame of explanatory variables (EVs) included in the model,
#'   with column names matching EV names. See \code{\link{readData}}.
#' @param EV Name or column index of the explanatory variable in \code{data} for
#'   which the response curve is to be generated.
#' @param transformation  Full pathway of the 'transformations.Rdata' file
#'   containing the transformations used to build the model. This file is saved
#'   as a result of the \code{\link{deriveVars}} function. Equivalently, the
#'   second item in the list returned by \code{\link{deriveVars}} can be used
#'   directly.
#' @param model Full pathway of the '.lambdas' file of the model in question.
#'   This file is saved as a result of \code{\link{selectEV}}.
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
#' responseEV1 <- plotResp2(dat, "EV1",
#'    transformation = "D:/path/to/modeling/directory/deriveVars/transformations.Rdata",
#'    model = "D:/path/to/modeling/directory/selectEV/round/model/1.lambdas")
#'
#' # From vignette:
#' pr_bygallResp2 <- plotResp2(grasslandPO, "pr_bygall",
#'    transformation = system.file("extdata", "transformations.Rdata", package = "maxentmodelselectr"),
#'    model = system.file("extdata", "1.lambdas", package = "maxentmodelselectr"))
#' }
#'
#' @export


plotResp2 <- function(data, EV, transformation, model, logscale = FALSE,
                         ...) {

  lambdas <- utils::read.csv(model, header = FALSE)
  dvnames <- as.character(lambdas[1:(nrow(lambdas)-4), 1])
  dvnamesni <- dvnames[grep(":", dvnames, invert = TRUE)]
  dvnamesi <- dvnames[grep(":", dvnames)]

  check <- lapply(dvnamesni, function(x) { startsWith(x, colnames(data)) })
  if (any(sapply(check, function(x) { all(x==FALSE) }))) {
    stop("All EVs in the model must be represented in data",
      call. = FALSE)
  }

  alltransf <- .load.transf(transformation)

  evnames <- unique(unname(sapply(dvnamesni, function(x) {
    colnames(data)[startsWith(x, colnames(data))]
  })))
  evnames <- sapply(evnames, .best.match, b = dvnamesni)

  check2 <- lapply(dvnamesni, function(x) { startsWith(names(alltransf), x) })
  if (any(sapply(check2, sum) < 1)) {
    stop("All DVs in the model must be represented in transformation",
      call. = FALSE)
  }

  EV <- colnames(data[, EV, drop = FALSE])
  evdata <- data[, EV]
  xnull <- environment(alltransf[startsWith(names(alltransf), EV)][[1]])$xnull
  if (class(xnull) %in% c("numeric", "integer")) {
    L <- (evdata - range(xnull)[1]) / diff(range(xnull))
    Range <- (range(L))
  }
  if (class(xnull) %in% c("factor", "character")) {
    if (all(evdata %in% xnull)) {Range <- "inside"} else {Range <- "outside"}
  }

  reg <- regexpr(EV, dvnamesni)
  margdvnamesni <- dvnamesni[which(startsWith(dvnamesni, EV) &
      attr(reg, "match.length") == max(attr(reg, "match.length")))]
  cnstdvnamesni <-  dvnamesni[!(dvnamesni %in% margdvnamesni)]

  margdvdatani <- lapply(margdvnamesni, function(x) {
    evname <- evnames[startsWith(x, evnames)]
    evdata <- data[, evname]
    transffunction <- alltransf[.best.match.ind(names(alltransf), x)][[1]]
    y <- transffunction(evdata)
    return(y)
  })
  names(margdvdatani) <- margdvnamesni

  cnstdvdatani <- lapply(cnstdvnamesni, function(x) {
    xnull <- environment(alltransf[startsWith(names(alltransf), x)][[1]])$xnull
    if (class(xnull) %in% c("numeric", "integer")) {
      evdata <- rep(mean(xnull), nrow(data))
    }
    if (class(xnull) %in% c("factor", "character")) {
      evdata <- rep(names(which.max(table(xnull))), nrow(data))
    }
    transffunction <- alltransf[.best.match.ind(names(alltransf), x)][[1]]
    y <- transffunction(evdata)
    return(y)
  })
  names(cnstdvdatani) <- cnstdvnamesni

  dvdatani <- c(margdvdatani, cnstdvdatani)

  prodlist <- strsplit(dvnamesi, ":")
  dvdatai <- lapply(prodlist, function(x) {
    dvdatani[[x[1]]] * dvdatani[[x[2]]]
  })
  names(dvdatai) <- dvnamesi

  dvdf <- data.frame(c(dvdatani, dvdatai), check.names = FALSE)
  modelfunction <- modelfromlambdas(model)
  respPts <- data.frame(EV = evdata, PRO = modelfunction(dvdf)[, 1])
  if (logscale == TRUE) {respPts$PRO <- log10(respPts$PRO)}

  if (class(respPts[, 1]) %in% c("numeric", "integer")) {
    graphics::plot(respPts[, 2] ~ respPts[, 1], ...,
      main = paste0("Marginal-effect response plot: ", EV), xlab = EV,
      ylab = ifelse(logscale == TRUE, "log Probability Ratio Output (logPRO)",
        "Probability Ratio Output (PRO)"))

    intervals <- min(c(ceiling(nrow(respPts) / 50), 100))
    respPts$int <- .reg.interval(respPts[, 1], intervals)
    grouped <- dplyr::group_by(respPts, int)
    respLine <- as.data.frame(dplyr::summarise(grouped, n = n(),
      intEV = mean(EV, na.rm = TRUE),
      intPRO = mean(PRO, na.rm = TRUE)))
    respLine$smoothPRO <- .ewma(respLine$intPRO, 5)
    graphics::lines(respLine$smoothPRO ~ respLine$intEV, col="red", lwd = 2)

    if (logscale == TRUE) {
      graphics::abline(h = 0, lty = 3)
    } else {
      graphics::abline(h = 1, lty = 3)
      }

    result <- list("respPts" = respPts, "respLine" = respLine)
  }

  if (class(respPts[, 1]) %in% c("factor", "character")) {
    respBar <- as.data.frame(dplyr::summarise(dplyr::group_by(respPts, EV),
      n = n(), intPRO = mean(PRO, na.rm = TRUE)))
    graphics::barplot(respBar[, 3], names.arg = respBar[, 1], ...,
      main = paste0("Marginal-effect response plot: ", EV), xlab = EV,
      ylab = ifelse(logscale == TRUE, "log Probability Ratio Output (logPRO)",
        "Probability Ratio Output (PRO)"))
    if (logscale == TRUE) {
      graphics::abline(h = 0, lty = 3)
    } else {
      graphics::abline(h = 1, lty = 3)
    }
    result <- respBar
  }

  return(result)
}
