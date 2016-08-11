#' Read in data object from files.
#'
#' \code{readData} reads in occurrence data in CSV file format and environmental
#' data in ASCII raster file format and produces a data object which can be used
#' as the starting point for the functions in this package. This function is
#' intended to make reading in data easy for users familiar with the maxent.jar
#' program. It is emphasized that important considerations for data preparation
#' (e.g. cleaning, sampling bias removal, etc.) are not treated in this package
#' and must be dealt with separately!
#'
#' When \code{occurrence} represents presence-only data (\code{PA = FALSE}), all
#' rows with values other than 'NA' in column 1 of the CSV file are treated as
#' presence locations. If column 1 contains any values of 'NA', these rows are
#' treated as the unknown background locations for the response variable. Thus,
#' 'NA' can be used to specify a specific set of background locations if
#' desired. Otherwise background points are randomly selected from the full
#' extent of the raster cells which are not already included as presence
#' locations. Only cells which contain data for all environmental variables are
#' selected as background locations.
#'
#' When \code{occurrence} represents presence/absence data (\code{PA = TRUE}),
#' rows with value '0' in column 1 of the CSV are treated as absence locations,
#' rows with value 'NA' are excluded, and all other rows are treated as
#' presences.
#'
#' The names of the ASCII raster files are used as the names of the explanatory
#' variables, so these files should be uniquely named, and the names must not
#' contain spaces or colons.
#'
#' @param occurrence Full pathway of the '.csv' file of occurrence data. The
#'   first column of the CSV should code occurrence (see Details), while the
#'   second and third columns should contain X and Y coordinates corresponding
#'   to the ASCII raster coordinate system. The first row is read as a header
#'   row.
#' @param contEV Pathway to a directory containing continuous environmental
#'   variables in '.asc' file format.
#' @param catEV Pathway to a directory containing categorical environmental
#'   variables in '.asc' file format.
#' @param maxbkg Integer. Maximum number of grid cells randomly selected as
#'   unknown background points for the response variable. Default is 10,000.
#'   Irrelevant for presence/absence data (\code{PA = TRUE}) and ignored for
#'   presence-only data (\code{PA = FALSE}) if \code{occurrence} contains 'NA'
#'   values.
#' @param PA Logical. Does \code{occurrence} represent presence/absence data?
#'   This argument affects how the values in \code{occurrence} are interpreted,
#'   and controls what type of data object is produced. See details.
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
  occ <- utils::read.csv(occurrence, header = TRUE, na.strings = "NA")
  contfiles <- list.files(contEV, pattern = "\\.asc$", full.names = TRUE)
  catfiles <- list.files(catEV, pattern = "\\.asc$", full.names = TRUE)
  stack <- raster::stack(c(contfiles, catfiles))

  if (PA == FALSE) {
    if (any(is.na(occ[, 1]))) {
      presev <- raster::extract(stack, occ[!is.na(occ[, 1]), 2:3])
      presdf <- data.frame(RV = rep(1, nrow(presev)), presev)
      bkgdf <- data.frame(RV = NA, raster::extract(stack,
        occ[is.na(occ[, 1]), 2:3]))
    } else {
      pres <- raster::extract(stack, occ[, 2:3], cellnumbers = TRUE)
      presdf <- data.frame(RV = 1, pres[,-1])
      bkg <- raster::extract(stack, raster::extent(stack), cellnumbers = TRUE)
      bkg <- bkg[!apply(bkg, 1, function(x) any(is.na(x))), ]
      bkg <- bkg[!(bkg[, 1] %in% pres[, 1]), ]
      if (nrow(bkg) > maxbkg) {
        bkgdf <- data.frame(RV = NA, bkg[sample(nrow(bkg), maxbkg), -1])
      } else {
        bkgdf <- data.frame(RV = NA, bkg[, -1])
      }
    }
    data <- rbind(presdf, bkgdf)
  }

  if (PA == TRUE) {
    presabs <- occ[!is.na(occ[, 1]), ]
    absindex <- presabs[, 1] == 0
    presev <- raster::extract(stack, presabs[!absindex, 2:3])
    presdf <- data.frame(RV = rep(1, nrow(presev)), presev)
    absev <- raster::extract(stack, presabs[absindex, 2:3])
    absdf <- data.frame(RV = rep(0, nrow(absev)), absev)
    data <- rbind(presdf, absdf)
  }

  catindex <- seq(ncol(data) - length(catfiles) + 1, ncol(data))
  data[catindex] <- lapply(data[catindex], function(x) as.factor(x))
  return(data)
}