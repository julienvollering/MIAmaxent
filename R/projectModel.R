#' Project model to data.
#'
#' \code{projectModel} Calculates the probability ratio output (PRO) of a given
#' model for any points where values of the explanatory variables in the model
#' are known.The transformations performed on the explanatory variables to build
#' the model must be specified.
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
#' @param data Data frame of all the explanatory variables (EVs) included in the
#'   model, with column names matching EV names. See \code{\link{readData}}.
#' @param transformation Full pathway of the 'transformations.Rdata' file
#'   containing the transformations used to build the model. This file is saved
#'   as a result of the \code{\link{deriveVars}} function. Equivalently, the
#'   second item in the list returned by \code{\link{deriveVars}} can be used
#'   directly.
#' @param model Full pathway of the '.lambdas' file of the model in question.
#'   This file is saved as a result of \code{\link{selectEV}}.
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
#'   \item A data frame showing the range of \code{data} compared to the
#'   training data, on a 0-1 scale.}
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
#' @examples
#' \dontrun{
#' modeloutput <- projectModel(newdat,
#'    transformation = "D:/path/to/modeling/directory/deriveVars/transformations.Rdata",
#'    model = "D:/path/to/modeling/directory/selectEV/round/model/1.lambdas")
#' }
#'
#' proj <- projectModel(toydata_sp1po, toydata_dvs$transformations,
#'    system.file("extdata/sommerfeltia", "1.lambdas", package = "MIAmaxent"))
#' proj
#'
#' \dontrun{
#' # From vignette:
#' grasslandPrediction <- projectModel(grasslandPO,
#'    transformation = grasslandDVs[[2]],
#'    model = system.file("extdata", "1.lambdas", package = "MIAmaxent"))
#' head(grasslandPrediction$output)
#' grasslandPrediction$ranges
#'
#' # From vignette:
#' library(raster)
#' contfiles <- list.files(system.file("extdata", "EV_continuous", package = "MIAmaxent"),
#'    full.names = TRUE)
#' catfiles <- list.files(system.file("extdata", "EV_categorical", package = "MIAmaxent"),
#'    full.names = TRUE)
#' stack <- raster::stack(c(contfiles, catfiles))
#' stackpts <- rasterToPoints(stack)
#' spatialPrediction <- projectModel(stackpts,
#'    transformation = grasslandDVs[[2]],
#'    model = system.file("extdata", "1.lambdas", package = "MIAmaxent"))
#' Predictionraster <- raster(stack, layer=0)
#' Predictionraster <- rasterize(spatialPrediction$output[, c("x", "y")], Predictionraster,
#'    field = spatialPrediction$output$PRO)
#' plot(Predictionraster, colNA="black")
#' }
#'
#' @export


projectModel <- function(data, transformation, model, clamping = FALSE,
                         rescale = FALSE, raw = FALSE) {

  lambdas <- utils::read.csv(model, header = FALSE)
  dvnames <- as.character(lambdas[1:(nrow(lambdas)-4), 1])
  dvnamesni <- dvnames[grep(":", dvnames, invert = TRUE)]
  dvnamesi <- dvnames[grep(":", dvnames)]

  .check.dvs.in.data(dvnamesni, data)
  evnames <- unique(sub("_.*", "", dvnamesni))

  alltransf <- .load.transf(transformation)
  .check.dvs.in.transf(dvnamesni, alltransf)

  Ranges <- lapply(evnames, function(x) {
    evdata <- data[, x]
    xnull <- environment(alltransf[x == sub("_.*", "", names(alltransf))][[1]])$xnull
    if (class(xnull) %in% c("numeric", "integer")) {
      L <- (evdata - range(xnull)[1]) / diff(range(xnull))
      return(range(L))
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
  names(dvdatani) <- dvnamesni

  prodlist <- strsplit(dvnamesi, ":")
  dvdatai <- lapply(prodlist, function(x) {
    dvdatani[[x[1]]] * dvdatani[[x[2]]]
  })
  names(dvdatai) <- dvnamesi

  dvdf <- data.frame(c(dvdatani, dvdatai), check.names = FALSE)
  modelfunction <- modelFromLambdas(model, raw)
  modeloutput <- modelfunction(dvdf)[, 1]
  if (rescale == TRUE) {
    if (raw == TRUE) {
      modeloutput <- modeloutput/sum(modeloutput)
    } else {
      modeloutput <- (modeloutput/sum(modeloutput))*length(modeloutput)
    }
  }

  if (raw == TRUE) {
    Output <- as.data.frame(cbind(raw = modeloutput, data))
  } else {
    Output <- as.data.frame(cbind(PRO = modeloutput, data))
  }

  return(list(output = Output, ranges = Ranges))
}
