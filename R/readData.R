#' Read in data from files.
#'
#' \code{readData} reads in occurence data in CSV format and environmental data
#' in ASCII raster format and produces a data object which can be used as the
#' starting point for many analyses in this package. This function is intended
#' to make reading in data easy for users familiar with the maxent.jar program.
#' It emphasized that important considerations for data preparation (e.g.
#' cleaning, sampling bias removal, etc.) are not treated in this package and
#' must be dealt with separately!
#'
#' Background points are randomly selected from the full extent of the grid
#' files, and may include presence locations. Only cells which contain data for
#' all variables are selected.
#'
#' @param presence Full pathway of the '.csv' file of occurrence data. The first
#'   column codes presence of the model target, while the second and third
#'   columns should contain X and Y coordinates cooresponding to the ASCII
#'   raster coordinate system. The first row is read as a header column.
#' @param contEV Pathway to a directory containing continuous environmental
#'   variables in '.asc' file format.
#' @param catEV Pathway to a directory containing categorical environmental
#'   variables in '.asc' file format.
#' @param maxbkg Integer. Maximum number of grid cells selected as unknown
#'   background points for the response varaible. Default is 10,000.
#'
#' @return Data frame with the Response Variable (RV) in the first column, and
#'   Explanatory Variables (EVs) in subsequent columns.
#'
#' @export


readData <- function(presence, contEV, catEV, maxbkg = 10000) {
  pres <- utils::read.csv(presence, header = TRUE)
  contfiles <- list.files(contEV, pattern = "\\.asc$", full.names = TRUE)
  catfiles <- list.files(catEV, pattern = "\\.asc$", full.names = TRUE)
  stack <- raster::stack(c(contfiles, catfiles))
  presdf <- data.frame(RV = 1, raster::extract(stack, pres[, 2:3]))
  bkg <- raster::as.matrix(stack)
  bkg <- bkg[!apply(bkg, 1, function(x) any(is.na(x))), ]
  if (nrow(bkg) > maxbkg) {
    bkgdf <- data.frame(RV = NA, bkg[sample(nrow(bkg), maxbkg), ])
  } else {
    bkgdf <- data.frame(RV = NA, bkg)
  }
 alldf <- rbind(presdf, bkgdf)
 catindex <- seq(ncol(alldf) - length(catfiles) + 1, ncol(alldf))
 alldf[catindex] <- lapply(alldf[catindex], function(x) as.factor(x))
 return(alldf)
}