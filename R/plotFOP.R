#' Plot Frequency of Observed Presence (FOP).
#'
#' \code{plotFOP} produces a Frequency of Observed Presence (FOP) plot for a
#' given explanatory variable. An FOP plot shows the rate of occurrence of the
#' response variable across intervals or levels of the explanatory variable. For
#' continuous variables, a local regression ("loess") of the FOP values is added
#' to the plot as a line. Data density is plotted in the background (grey) to
#' help visualize where FOP values are more or less certain.
#'
#' A list of the optimum EV value and a data frame containing the plotted data
#' is returned invisibly. Store invisibly returned output by assigning it to an
#' object.
#'
#' In the local regression ("loess"), the plotted FOP values are regressed
#' against their EV values. The points are weighted by the number of
#' observations they represent, such that an FOP value from an interval with
#' many observations is given more weight.
#'
#' For continuous variables, the returned value of 'EVoptimum' is based on the
#' loess-smoothed FOP values, such that a point maximum in FOP may not always be
#' considered the optimal value of EV.
#'
#' If the response variable in \code{data} represents presence/absence data, the
#' result is an empirical frequency of presence curve, rather than a observed
#' frequency of presence curve (see Støa et al. [2018], Sommerfeltia).
#'
#' @param data Data frame containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent either presence and background (coded as 1/NA) or presence
#'   and absence (coded as 1/0). See Details for information regarding
#'   implications of occurrence data type. See also \code{\link{readData}}.
#' @param EV Name or column index of the explanatory variable in \code{data} for
#'   which to calculate FOP.
#' @param span The proportion of FOP points included in the local regression
#'   neighborhood. Should be between 0 and 1. Irrelevant for categorical EVs.
#' @param intervals Number of intervals into which the continuous EV is divided.
#'   Defaults to the minimum of N/10 and 100. Irrelevant for categorical EVs.
#' @param ranging Logical. If \code{TRUE}, will range the EV scale to [0,1].
#'   This is equivalent to plotting FOP over the linear transformation produced
#'   by deriveVars. Irrelevant for categorical EVs.
#' @param densitythreshold Numeric. Intervals containing fewer than this number
#'   of observations will be represented with an open symbol in the plot.
#'   Irrelevant for categorical EVs.
#' @param ... Arguments to be passed to \code{plot} or \code{barplot} to control
#'   the appearance of the plot. For example: \itemize{ \item \code{lwd} for
#'   line width \item \code{cex.main} for size of plot title \item \code{space}
#'   for space between bars }
#'
#' @return In addition to the graphical output, a list of 2: \enumerate{ \item
#'   \code{EVoptimum}. The EV value (or level, for categorical EVs) at which FOP
#'   is highest \item  \code{FOPdata}. A data frame containing the plotted data.
#'   Columns in this data frame represent the following: EV interval ("int"),
#'   number of observations in the interval ("n"), mean EV value of the
#'   observations in the interval ("intEV"), mean RV value of the observations
#'   in the interval ("intRV"), and local regression predicted intRV ("loess").
#'   For categorical variables, only the level name ("level"), the number of
#'   observations in the level ("n"), and the mean RV value of the level
#'   ("levelRV") are used.}
#'
#' @references Støa, B., R. Halvorsen, S. Mazzoni, and V. I. Gusarov. (2018).
#'   Sampling bias in presence-only data used for species distribution
#'   modelling: theory and methods for detecting sample bias and its effects on
#'   models. Sommerfeltia 38:1–53.


#'
#' @examples
#' FOPev11 <- plotFOP(toydata_sp1po, 2)
#' FOPev12 <- plotFOP(toydata_sp1po, "EV12", intervals = 8)
#' FOPev12$EVoptimum
#' FOPev12$FOPdata
#'
#' \dontrun{
#' # From vignette:
#' teraspifFOP <- plotFOP(grasslandPO, "teraspif")
#' terslpdgFOP <- plotFOP(grasslandPO, "terslpdg")
#' terslpdgFOP <- plotFOP(grasslandPO, "terslpdg", span = 0.75, intervals = 20)
#' terslpdgFOP
#' geobergFOP <- plotFOP(grasslandPO, 10)
#' geobergFOP
#' }
#'
#' @export


plotFOP <- function(data, EV, span = 0.5, intervals = NULL, ranging = FALSE,
                    densitythreshold = NULL, ...) {

  if (EV==1) {
    stop("'EV' cannot be the first column of 'data', which must be the response variable")
  }
  df <- data.frame(RV = data[, 1], EV = data[, EV])
  evname <- names(data[, EV, drop = FALSE])

  .binaryrvcheck(df[, 1])
  df[, 1][is.na(df[, 1])] <- 0

  if (class(df[, 2]) %in% c("numeric", "integer")) {
    if (ranging == T) {
      df[, 2] <- (df[, 2] - min(df[, 2])) / diff(range(df[, 2]))
    }
    if (is.null(intervals)) {intervals <- min(c(ceiling(nrow(df)/10), 100))}
    df$int <- cut(df[, 2], breaks=max(2, intervals))

    grouped <- dplyr::group_by(df, int)
    FOPdf <- as.data.frame(dplyr::summarise(grouped, n = n(),
      intEV = mean(EV), intRV = mean(RV, na.rm=FALSE)))
    FOPdf$loess <- stats::predict(stats::loess(intRV~intEV, FOPdf,
                                               weights=FOPdf$n, span=span))

    if (any(is.na(FOPdf$loess))) {
      evoptimum <- FOPdf$intEV[which.max(FOPdf$intRV)]
    } else { evoptimum <- FOPdf$intEV[which.max(FOPdf$loess)]  }

    FOP <- list(EVoptimum = evoptimum,
                FOPdata = FOPdf)

    op <- graphics::par(mar=(c(5, 4, 4, 4) + 0.3))
    on.exit(graphics::par(op))
    dens <- stats::density(df[, 2])
    graphics::plot(range(dens$x), range(dens$y), type="n", axes=FALSE, ann=FALSE)
    graphics::polygon(x=c(min(dens$x), dens$x, max(dens$x)), y=c(0, dens$y, 0),
                      border=NA, col="grey90")
    graphics::axis(side=4, col="grey60", col.axis="grey60", las=1)
    graphics::mtext("Kernel estimated data density", side=4, line=3, col="grey60")
    graphics::par(new=TRUE)
    if (is.null(densitythreshold)) {densitythreshold <- 0}
    closedsymbols <- as.numeric(FOPdf$n >= densitythreshold)+1
    args1 <- list(bty="n", main = paste0("FOP plot: ", evname), xlab = evname,
                  ylab = "Frequency of Observed Presence (FOP)", las=1,
                  pch=c(1,20)[closedsymbols])
    inargs <- list(...)
    args1[names(inargs)] <- inargs
    do.call(graphics::plot, c(list(x=FOPdf$intEV, y=FOPdf$intRV), args1))
    graphics::points(FOPdf$loess ~ FOPdf$intEV, type="l", lwd=2, col="red")
  }

  if (class(df[, 2]) %in% c("factor", "character")) {
    grouped <- dplyr::group_by(df, EV)
    FOPdf <- as.data.frame(dplyr::summarise(grouped, n = n(),
                                            lvlRV = mean(RV, na.rm=FALSE)))

    FOP <- list(EVoptimum = FOPdf$EV[which.max(FOPdf$lvlRV)],
                FOPdata = data.frame(level=FOPdf$EV, n=FOPdf$n,
                                     levelRV=FOPdf$lvlRV))

    op <- graphics::par(mar=(c(5, 4, 4, 4) + 0.3))
    on.exit(graphics::par(op))
    graphics::barplot(FOPdf$n, axes=FALSE, ann=FALSE, col="grey90", border=NA)
    graphics::axis(side=4, col="grey60", col.axis="grey60", las=1)
    graphics::mtext("Number of observations in data", side=4, line=3, col="grey60")
    graphics::par(new=TRUE)
    args1 <- list(names.arg = FOPdf$EV, main = paste0("FOP plot: ", evname),
                  xlab = evname, ylab = "Frequency of Observed Presence (FOP)",
                  density=rep(20, nrow(FOPdf)), col="black", las = 2)
    inargs <- list(...)
    args1[names(inargs)] <- inargs
    do.call(graphics::barplot, c(list(FOPdf$lvlRV), args1))
  }

  invisible(FOP)
}
