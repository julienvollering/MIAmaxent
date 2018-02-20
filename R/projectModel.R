#' Project model to data.
#'
#' \code{projectModel} Calculates model predictions for any points where values
#' of the explanatory variables in the model are known. \code{projectModel} can
#' be used to get model predictions for the training data, or to project the
#' model to a new space or time.
#'
#' Missing data (NA) for a continuous variable will result in NA output for that
#' point. Missing data for a categorical variable is treated as belonging to
#' none of the categories.
#'
#' When \code{rescale = FALSE} the scale of the model output (PRO or raw)
#' returned by this function is dependent on the data used to train the model.
#' For example, a location with PRO = 2 can be interpreted as having a
#' probability of presence twice as high as an average site in the
#' \emph{training} data (Halvorsen, 2013, Halvorsen et al., 2015). When
#' \code{rescale = TRUE}, the output is linearly rescaled with respect to the
#' data onto which the model is projected. In this case, a location with PRO = 2
#' can be interpreted as having a probability of presence twice as high as an
#' average site in the \emph{projection} data. Similarly, raw values are on a
#' scale which is dependent on the size of either the training data extent
#' (\code{rescale = FALSE}) or projection data extent (\code{rescale = TRUE}).
#'
#' @param model The model to be projected, represented by an object of class
#'   'glm'. This may be the object returned by \code{\link{chooseModel}}, or the
#'   'selectedmodel' returned by \code{\link{selectEV}}.
#' @param transformations Transformation functions used to create the derived
#'   variables in the model. I.e. the 'transformations' returned by
#'   \code{\link{deriveVars}}. Equivalently, the full file pathway of the
#'   'transformations.Rdata' file saved as a result of \code{\link{deriveVars}}.
#' @param data Data frame of all the explanatory variables (EVs) included in the
#'   model (see \code{\link{readData}}). Alternatively, an object of class
#'   'RasterStack' containing rasters for all EVs included in the model. Column
#'   or raster names must match EV names.
#' @param clamping Logical. Do clamping \emph{sensu} Phillips et al. (2006).
#'   Default is \code{FALSE}.
#' @param rescale Logical. Linearly rescale model output (PRO or raw) with
#'   respect to the projection \code{data}? This has implications for the
#'   interpretation of output values with respect to reference values (e.g. PRO
#'   = 1). See details.
#' @param raw Logical. Return raw Maxent output instead of probability ratio
#'   output (PRO)? Default is FALSE.
#'
#' @return List of 2: \enumerate{ \item A data frame with the model output in
#'   column 1 and the corresponding explanatory data in subsequent columns.
#'   \item A list showing the range of \code{data} compared to the training
#'   data, on a 0-1 scale.} If \code{data} is a RasterStack, the model output is
#'   plotted as a raster, and only the list of ranges is returned.
#'
#'
#' @references Halvorsen, R. (2013) A strict maximum likelihood explanation of
#'   MaxEnt, and some implications for distribution modelling. Sommerfeltia, 36,
#'   1-132.
#' @references Halvorsen, R., Mazzoni, S., Bryn, A. & Bakkestuen, V. (2015)
#'   Opportunities for improved distribution modelling practice via a strict
#'   maximum likelihood interpretation of MaxEnt. Ecography, 38, 172-183.
#' @references Phillips, S.J., Anderson, R.P. & Schapire, R.E. (2006) Maximum
#'   entropy modeling of species geographic distributions. Ecological Modelling,
#'   190, 231-259.
#'
#' @export


projectModel <- function(model, transformations, data, clamping = FALSE,
                         rescale = FALSE, raw = FALSE) {

  dvnamesni <- names(model$betas)[grep(":", names(model$betas), invert = TRUE)]
  evnames <- unique(sub("_.*", "", dvnamesni))

  map <- FALSE
  if (class(data) == "RasterStack") {
    map <- TRUE
    evstack <- data[[evnames]]
    data <- raster::as.data.frame(evstack, na.rm = TRUE)
    cells <- as.numeric(row.names(data))
  }

  for (i in evnames) {
    if (sum(colnames(data) == i) != 1) {
      stop(paste(a, "must be represented in 'data' (exactly once)"),
           call. = FALSE) }
  }

  alltransf <- .load.transf(transformations)
  .check.dvs.in.transf(dvnamesni, alltransf)

  Ranges <- lapply(evnames, function(x) {
    evdata <- data[, x]
    anevtransf <- alltransf[grepl(paste0(x, "_"), names(alltransf))][[1]]
    xnull <- environment(anevtransf)$xnull
    if (class(xnull) %in% c("numeric", "integer")) {
      L <- (evdata - range(xnull)[1]) / diff(range(xnull))
      return(range(L, na.rm = TRUE))
    }
    if (class(xnull) %in% c("factor", "character")) {
      if (all(evdata %in% xnull)) {return("inside")} else {return("outside")}
    }
  })
  names(Ranges) <- evnames

  dvdatani <- lapply(dvnamesni, function(x) {
    evdata <- data[, sub("_.*", "", x)]
    y <- alltransf[[paste0(x, "_transf")]](evdata)
    if (clamping == TRUE) {
      y[y > 1] <- 1
      y[y < 0] <- 0
    }
    return(y)
  })
  dvdatani <- as.data.frame(do.call(cbind, dvdatani))
  names(dvdatani) <- dvnamesni

  mmformula <- stats::update.formula(model$formula.narm, NULL ~ . - 1)
  newdata <- model.matrix(mmformula, dvdatani)
  newdata <- newdata[match(rownames(dvdatani), rownames(newdata)), ]
  rownames(newdata) <- rownames(dvdatani)
  if (raw == TRUE) {
    modeloutput <- exp((newdata %*% model$betas) + model$alpha)
  } else {
    Ntrain <- length(alltransf[[1]])
    modeloutput <- exp((newdata %*% model$betas) + model$alpha) * Ntrain
  }

  if (rescale == TRUE) {
    if (raw == TRUE) {
      modeloutput <- modeloutput/sum(modeloutput, na.rm = TRUE)
    } else {
      modeloutput <- (modeloutput/sum(modeloutput, na.rm = TRUE)) *
        sum(!is.na(modeloutput))
    }
  }

  if (raw == TRUE) {
    Output <- data.frame("raw" = modeloutput, data)
  } else {
    Output <- data.frame("PRO" = modeloutput, data)
  }

  if (map == TRUE) {
    values <- rep(NA, raster::ncell(evstack))
    values[cells] <- Output[,1]
    outraster <- evstack[[1]]
    outraster <- raster::setValues(outraster, values)
    raster::plot(outraster)
    return(Ranges)
  } else {
    return(list(output = Output, ranges = Ranges))
  }

}
