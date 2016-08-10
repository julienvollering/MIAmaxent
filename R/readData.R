#' Read in data object from files.
#'
#' \code{readData} reads in occurrence data in CSV format and environmental data
#' in ASCII raster format and produces a data object which can be used as the
#' starting point for the functions in this package. This function is intended
#' to make reading in data easy for users familiar with the maxent.jar program.
#' It is emphasized that important considerations for data preparation (e.g.
#' cleaning, sampling bias removal, etc.) are not treated in this package and
#' must be dealt with separately!
#'
#' When \code{occurrence} represents presence-only data (\code{PA = FALSE}), all
#' rows are treated as presences, and the values in column 1 of the CSV are
#' irrelevant. Conversely, when \code{occurrence} represents presence/absence
#' data (\code{PA = TRUE}), rows with value 0 in column 1 of the CSV are treated
#' as absences, and all other rows are treated as presences.
#'
#' For presence-only occurrence data, background points are randomly selected
#' from the full extent of the grid files, and may include presence locations.
#' Only cells which contain data for all variables are selected.
#'
#' The names of the ASCII raster files are used as the names of the explanatory
#' variables, so these files should be uniquely named, and the names must not
#' contain spaces or colons.
#'
#' @param occurrence Full pathway of the '.csv' file of occurrence data. The
#'   first column of the CSV should code occurrence (see Details), while the
#'   second and third columns should contain X and Y coordinates cooresponding
#'   to the ASCII raster coordinate system. The first row is read as a header
#'   row.
#' @param contEV Pathway to a directory containing continuous environmental
#'   variables in '.asc' file format.
#' @param catEV Pathway to a directory containing categorical environmental
#'   variables in '.asc' file format.
#' @param maxbkg Integer. Maximum number of grid cells selected as unknown
#'   background points for the response variable. Default is 10,000.
#' @param PA Logical. Does \code{occurrence} represent presence/absence data?
#'   When \code{TRUE}, rows in \code{occurrence} with the value 0 in the first
#'   column are treated as absences, and all others are treated as presences.
#'
#' @return Data frame with the Response Variable (RV) in the first column, and
#'   Explanatory Variables (EVs) in subsequent columns. When \code{PA = FALSE},
#'   RV values are 1/NA, and when \code{PA = TRUE}, RV values are 1/0.
#'
#'   With presence-only occurrence data, the returned output can be used as the
#'   \code{data} argument for \code{\link{plotFOP}}, \code{\link{deriveVars}},
#'   \code{\link{selectDVforEV}}, \code{\link{selectEV}}, and
#'   \code{\link{plotResp}}. With presence/absence occurrence data, the returned
#'   output can be used as the \code{data} argument for \code{\link{testAUC}}.
#'   Output from \code{readData} can also be used as the \code{data} argument
#'   for \code{\link{plotResp2}}, and \code{\link{projectModel}}, but for these
#'   functions the values of RV are irrelevant.
#'
#' @export


readData <- function(occurrence, contEV, catEV, maxbkg = 10000, PA = FALSE) {
  occ <- utils::read.csv(occurrence, header = TRUE)
  contfiles <- list.files(contEV, pattern = "\\.asc$", full.names = TRUE)
  catfiles <- list.files(catEV, pattern = "\\.asc$", full.names = TRUE)
  stack <- raster::stack(c(contfiles, catfiles))

  if (PA == FALSE) {
    presdf <- data.frame(RV = 1, raster::extract(stack, occ[, 2:3]))
    bkg <- raster::as.matrix(stack)
    bkg <- bkg[!apply(bkg, 1, function(x) any(is.na(x))), ]
    if (nrow(bkg) > maxbkg) {
      bkgdf <- data.frame(RV = NA, bkg[sample(nrow(bkg), maxbkg), ])
    } else {
      bkgdf <- data.frame(RV = NA, bkg)
    }
    data <- rbind(presdf, bkgdf)
    catindex <- seq(ncol(data) - length(catfiles) + 1, ncol(data))
    data[catindex] <- lapply(data[catindex], function(x) as.factor(x))
  }

  if (PA == TRUE) {
    data <- data.frame(RV = NA, raster::extract(stack, occ[, 2:3]))
    absindex <- occ[, 1] == 0
    data$RV <- as.numeric(!absindex)
  }

 return(data)
}