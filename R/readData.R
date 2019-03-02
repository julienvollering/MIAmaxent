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
#' treated as the uninformed background locations. Thus, 'NA' can be used to
#' specify a specific set of uninformed background locations if desired.
#' Otherwise uninformed background locations are randomly selected from the full
#' extent of the raster cells which are not already included as presence
#' locations. Only cells which contain data for all environmental variables are
#' retained as presence locations or selected as uninformed background
#' locations.
#'
#' When \code{occurrence} represents presence/absence data (\code{PA = TRUE}),
#' rows with value '0' in column 1 of the CSV are treated as absence locations,
#' rows with value 'NA' are excluded, and all other rows are treated as
#' presences.
#'
#' The names of the ASCII raster files are used as the names of the explanatory
#' variables, so these files should be uniquely named. \code{readData} replaces
#' underscores '_', spaces ' ' and other special characters not allowed in names
#' with periods '.'. In MIAmaxent, underscores and colons are reserved to denote
#' derived variables and interaction terms, respectively.
#'
#' @param occurrence Full pathway of the '.csv' file of occurrence data. The
#'   first column of the CSV should code occurrence (see Details), while the
#'   second and third columns should contain X and Y coordinates corresponding
#'   to the ASCII raster coordinate system. The first row of the csv is read as
#'   a header row.
#' @param contEV Pathway to a directory containing continuous environmental
#'   variables in '.asc' file format.
#' @param catEV Pathway to a directory containing categorical environmental
#'   variables in '.asc' file format.
#' @param maxbkg Integer. Maximum number of grid cells randomly selected as
#'   uninformed background locations for the response variable. Default is
#'   10,000. Irrelevant for presence/absence data (\code{PA = TRUE}) and ignored
#'   for presence-only data (\code{PA = FALSE}) if \code{occurrence} contains
#'   'NA' values. See Details.
#' @param PA Logical. Does \code{occurrence} represent presence/absence data?
#'   This argument affects how the values in \code{occurrence} are interpreted,
#'   and controls what type of data object is produced. See Details.
#' @param XY Logical. Include XY coordinates in the output. May be useful for
#'   spatial plotting. Note that coordinates included in the training data used
#'   to build the model will be treated as explanatory variables.
#'
#' @return Data frame with the Response Variable (RV) in the first column, and
#'   Explanatory Variables (EVs) in subsequent columns. When \code{PA = FALSE},
#'   RV values are 1/NA, and when \code{PA = TRUE}, RV values are 1/0.
#'
#' @examples
#' toydata_sp1po <- readData(system.file("extdata/sommerfeltia", "Sp1.csv", package = "MIAmaxent"),
#'    contEV = system.file("extdata/sommerfeltia", "EV_continuous", package = "MIAmaxent"))
#' toydata_sp1po
#'
#' \dontrun{
#' # From vignette:
#' grasslandPO <- readData(
#'  occurrence=system.file("extdata", "occurrence_PO.csv", package="MIAmaxent"),
#'   contEV=system.file("extdata", "EV_continuous", package="MIAmaxent"),
#'   catEV=system.file("extdata", "EV_categorical", package="MIAmaxent"),
#'   maxbkg=20000)
#' str(grasslandPO)
#'
#' # From vignette:
#' grasslandPA <- readData(
#'   occurrence = system.file("extdata", "occurrence_PA.csv", package="MIAmaxent"),
#'   contEV = system.file("extdata", "EV_continuous", package="MIAmaxent"),
#'   catEV = system.file("extdata", "EV_categorical", package="MIAmaxent"),
#'   PA = TRUE, XY = TRUE)
#' head(grasslandPA)
#' tail(grasslandPA)
#' }
#'
#' @export
#'


readData <- function(occurrence, contEV = NULL, catEV = NULL, maxbkg = 10000,
                     PA = FALSE, XY = FALSE) {

  if (is.null(contEV) && is.null(catEV)) {
    stop("Please specify at least 1 explanatory variable directory.
       (contEV, catEV, or both)", call. = FALSE)
  }

  occ <- utils::read.csv(occurrence, header = TRUE, na.strings = "NA")
  if (is.null(contEV)) {contfiles <- NULL} else {
    contfiles <- list.files(contEV, pattern = "\\.asc$", full.names = TRUE)
  }
  if (is.null(catEV)) {catfiles <- NULL} else {
    catfiles <- list.files(catEV, pattern = "\\.asc$", full.names = TRUE)
  }
  stack <- raster::stack(c(contfiles, catfiles))
  names(stack) <- gsub(".asc", "", basename(c(contfiles, catfiles)))

  if (PA == FALSE) {
    if (any(is.na(occ[, 1]))) {
      pres <- raster::extract(stack, occ[!is.na(occ[, 1]), 2:3], cellnumbers = TRUE)
      pres <- pres[!duplicated(pres[, 1]), ]
      bkg <- raster::extract(stack, occ[is.na(occ[, 1]), 2:3], cellnumbers =TRUE)
      bkg <- bkg[!duplicated(bkg[, 1]), ]
      cells <- c(pres[, 1], bkg[, 1])
      presbkg <- rbind(pres[, -1], bkg[, -1])
    } else {
      pres <- raster::extract(stack, occ[, 2:3], cellnumbers = TRUE)
      pres <- pres[!duplicated(pres[, 1]), ]
      bkg <- raster::as.data.frame(stack, na.rm = TRUE)
      bkg <- bkg[!(rownames(bkg) %in% pres[,1]),]
      if (nrow(bkg) > maxbkg) {bkg <- bkg[sample(nrow(bkg), maxbkg), ]}
      cells <- c(pres[,1], as.numeric(rownames(bkg)))
      presbkg <- rbind(pres[, -1], bkg)
    }
    xy <- raster::xyFromCell(stack, cells)
    data <- data.frame("RV"=c(rep(1, nrow(pres)), rep(NA, nrow(bkg))), xy,
                       presbkg)
  }

  if (PA == TRUE) {
    presabs <- occ[!is.na(occ[, 1]), ]
    absindex <- presabs[, 1] == 0
    pres <- raster::extract(stack, presabs[!absindex, 2:3], cellnumbers = TRUE)
    pres <- pres[!duplicated(pres[, 1]), ]
    abs <- raster::extract(stack, presabs[absindex, 2:3], cellnumbers = TRUE)
    abs <- abs[!duplicated(abs[, 1]), ]
    xy <- raster::xyFromCell(stack, c(pres[, 1], abs[, 1]))
    presabs <- rbind(pres[, -1], abs[, -1])
    data <- data.frame("RV"=c(rep(1, nrow(pres)), rep(0, nrow(abs))), xy,
                       presabs)
  }

  if (XY == FALSE) {
    data <- data[, -c(2:3)]
  }
  if (!is.null(catEV)) {
    catindex <- seq(ncol(data) - length(catfiles) + 1, ncol(data))
    data[catindex] <- lapply(data[catindex], function(x) as.factor(x))
  }
  colnames(data) <- make.names(colnames(data), allow_ = FALSE)

  data <- data[apply(data[,2:ncol(data)], 1, function(x) !any(is.na(x))), ]
  rownames(data) <- seq(length=nrow(data))

  return(data)
}