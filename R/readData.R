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
#' DESCRIPTION: Import maptools (readAsciiGrid), and sp (spatial classes)
#'
#' @param presence Full pathway of the '.csv' file of occurrence data. The first
#'   column codes presence of the model target, while the second and third
#'   columns should contain X and Y coordinates cooresponding to the ASCII
#'   raster coordinate system. The first row is read as a header column.
#' @param contEV Pathway to a directory containing continuous environmental
#'   variables is '.asc' file format.
#' @param catEV Pathway to a directory containing categorical environmental
#'   variables is '.asc' file format.
#'
#' @return Data frame with the Response Variable (RV) in the first column, and
#'   Explanatory Variables (EVs) in subsequent columns.
#'
#' @export


readData <- function(presence, contEV, catEV) {
  presence <- read.csv(presence, header = TRUE)
}