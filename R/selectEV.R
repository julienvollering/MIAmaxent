#' Select parsimonious set of explanatory variables.
#'
#' \code{selectEV} selects the parsimonious set of explanatory variables (EVs)
#' which best explains variation in a given response variable (RV). Each EV can
#' be represented by 1 or more derived variables (see \code{\link{deriveVars}}).
#' The function uses a process of forward selection based on comparison of
#' nested models by the F-test. An EV is selected for inclusion when, during
#' nested model comparison, it accounts for a significant amount of remaining
#' variation, under the alpha value specified by the user.
#'
#' The F-statistic that \code{selectEV} uses for nested model comparison is
#' calculated using equation 59 in Halvorsen (2013). See Halvorsen et al. (2015)
#' for a more detailed explanation of the forward selection procedure.
#'
#' When \code{interaction = TRUE}, the forward selection procedure selects a
#' parsimonious group of individual EVs first, and then tests interactions
#' between EVs included in the model afterwards. Therefore, interactions are
#' only explored between terms which are individually explain a significant
#' amount of variation. When \code{interaction = FALSE}, interactions are not
#' considered.
#'
#' If \code{trainmax} reduces the number of uninformed background points in the
#' training data, a new \code{data} object is returned as part of the function
#' output. This \code{data} object shows which of the uninformed background
#' points were randomly selected, and should be used together with the selected
#' EVs in \code{\link{plotResp}} if plotting single-effect model response.
#'
#' Variables should be uniquely named, and the names must not contain spaces or
#' colons. Colons are reserved as the designator for interaction terms.
#'
#' @param data Data frame containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent presence/background data, coded as: 1/NA. See
#'   \code{\link{readData}}.
#' @param dvdata List of data frames, with each data frame containing a
#'   parsimonious group of derived variables for a given explanatory variable
#'   (e.g. the first item in the list returned by \code{\link{selectDVforEV}}).
#' @param alpha Alpha-level used in F-test comparison of models. Default is
#'   0.01.
#' @param interaction Logical. Allows interaction terms between pairs of EVs.
#'   Default is \code{TRUE}.
#' @param dir Directory to which files will be written during subset selection
#'   of explanatory variables. Defaults to the working directory.
#' @param trainmax Integer. Maximum number of uninformed background points to be
#'   used to train the models. May be used to reduce computation time for data
#'   sets with very large numbers of points. Default is no maximum. See Details
#'   for more information.
#'
#' @return List of 2 (3): \enumerate{ \item A list of data frames, with one data
#'   frame for each \emph{selected} EV. This item is recommended as input for
#'   \code{dvdata} in \code{\link{plotResp}}. \item A data frame showing the
#'   trail of forward selection of individual EVs (and interaction terms if
#'   necessary). \item (If \code{trainmax} reduces the number of uninformed
#'   background points) a new \code{data} object. See details.}
#'
#' @references Halvorsen, R. (2013). A strict maximum likelihood explanation of
#'   MaxEnt, and some implications for distribution modelling. Sommerfeltia, 36,
#'   1-132.
#' @references Halvorsen, R., Mazzoni, S., Bryn, A., & Bakkestuen, V. (2015).
#'   Opportunities for improved distribution modelling practice via a strict
#'   maximum likelihood interpretation of MaxEnt. Ecography, 38(2), 172-183.
#'
#' @export


selectEV <- function(data, dvdata, alpha = 0.01, interaction = TRUE, dir = NULL,
                     trainmax = NULL) {

  rv <- data[, 1]
  .binaryrvcheck(rv)

  if (is.null(dir)) { dir <- getwd()}

  fdir <- file.path(dir, "selectEV")
  if (file.exists(fdir)) {
    stop("The specified dir already contains a selection of EVs.
Please specify a different dir. \n ", call. = FALSE)
  } else {
    dir.create(fdir)
  }

  returndata <- FALSE
  if (!is.null(trainmax) && trainmax < sum(is.na(rv))) {
    nub <- min(sum(is.na(rv)), trainmax)
    trainindex <- c(which(!is.na(rv)), sample(which(is.na(rv)), nub))
    data <- data[trainindex, ]
    dvdata <- lapply(dvdata, function(x) {x[trainindex, , drop = FALSE]})
    rv <- rv[trainindex]
    returndata <- TRUE
  }

  message(paste0("Forward selection of ", length(dvdata), " EVs"))

  result <- .parsevs(rv, dvdata, alpha, interaction, fdir)
  utils::write.csv(result[[2]], file = file.path(fdir, "evselection.csv"),
    row.names = FALSE)

  if (returndata == TRUE) {
    Result <- list(selectedEV = result[[1]], selection = result[[2]], data = data)
  } else {
    Result <- list(selectedEV = result[[1]], selection = result[[2]])
  }
  selectedEV <- result[[1]]
  save(selectedEV, file = file.path(fdir, "selectedEV.Rdata"))


  return(Result)
}
