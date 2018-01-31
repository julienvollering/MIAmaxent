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
#' Explanatory variables should be uniquely named, and the names must not
#' contain spaces, underscores, or colons. Underscores and colons are reserved
#' to denote derived variables and interaction terms repectively.
#'
#' @param data Data frame containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent presence/background data, coded as: 1/NA. See
#'   \code{\link{readData}}.
#' @param dvdata List of data frames, with each data frame containing
#'   \emph{selected} derived variables for a given explanatory variable (e.g.
#'   the first item in the list returned by \code{\link{selectDVforEV}}).
#' @param alpha Alpha-level used in F-test comparison of models. Default is
#'   0.01.
#' @param interaction Logical. Allows interaction terms between pairs of EVs.
#'   Default is \code{FALSE}.
#' @param dir Directory to which files will be written during subset selection
#'   of explanatory variables. Defaults to the working directory.
#' @param trainmax Integer. Maximum number of uninformed background points to be
#'   used to train the models. May be used to reduce computation time for data
#'   sets with very large numbers of points. Default is no maximum. See Details
#'   for more information.
#' @param formula A model formula (in the form y ~ x + ...) specifying a
#'   starting point for forward model selection. The independent terms in the
#'   formula will be included in the model regardless of explanatory power, and
#'   must be represented in \code{dvdata}, while the remaining explanatory
#'   variables in \code{dvdata} are candidates for selection. The first column
#'   in \code{data} is still taken as the response variable, regardless of
#'   \code{formula}. Default is \code{NULL}, meaning that forward selection
#'   starts with zero selected variables.
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
#' @examples
#' \dontrun{
#' selectedevs <- selectEV(dat, selectedderiveddat, alpha = 0.0001,
#'    dir = "D:/path/to/modeling/directory", interaction = TRUE)
#'
#' # From vignette:
#' grasslandEVselect <- selectEV(grasslandPO, grasslandDVselect[[1]], alpha = 0.001,
#'    interaction = TRUE)
#' summary(grasslandDVselect[[1]])
#' length(grasslandDVselect[[1]])
#' summary(grasslandEVselect[[1]])
#' length(grasslandEVselect[[1]])
#' plot(grasslandEVselect$selection$round, grasslandEVselect$selection$addedFVA)
#' }
#'
#' @export


selectEV <- function(data, dvdata, alpha = 0.01, interaction = FALSE, dir = NULL,
                     trainmax = NULL, formula = NULL) {

  rv <- data[, 1]
  .binaryrvcheck(rv)

  if (is.null(dir)) { dir <- getwd()}

  fdir <- file.path(dir, "selectEV")
  if (file.exists(fdir)) {
    yn <- readline("The specified dir already contains a selectEV result. Overwrite this result? (y/n) ")
    if (yn == "y") {
      unlink(fdir, recursive = TRUE)
      Sys.sleep(1)
    }
    if (yn != "y") { return(message("Overwrite declined")) }
  }
  dir.create(fdir, recursive = TRUE)

  returndata <- FALSE
  if (!is.null(trainmax) && trainmax < sum(is.na(rv))) {
    nub <- min(sum(is.na(rv)), trainmax)
    trainindex <- c(which(!is.na(rv)), sample(which(is.na(rv)), nub))
    data <- data[trainindex, ]
    dvdata <- lapply(dvdata, function(x) {x[trainindex, , drop = FALSE]})
    rv <- rv[trainindex]
    returndata <- TRUE
  }

  if (!is.null(formula)) {
    .formulacheck(formula, dvdata)
  }

  if (!is.null(formula) && length(labels(stats::terms(formula))) != 0) {
    nterms <- length(labels(stats::terms(formula)))
    message(paste0("Forward selection of ", length(dvdata) - nterms, " EVs"))
  } else {
    message(paste0("Forward selection of ", length(dvdata), " EVs"))
  }

  result <- .parsevs(rv, dvdata, alpha, interaction, fdir, formula)
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
