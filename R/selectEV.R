#' Select parsimonious sets of explanatory variables.
#'
#' \code{selectEV} selects the parsimonious set of explanatory variables (EVs)
#' which best explains variation in a given response variable (RV). Each EV can
#' be represented by 1 or more derived variables (see \code{\link{deriveVars}}).
#' The function uses a process of forward selection based on comparison of
#' nested models by the F-test, where the F-statistic is calculated using
#' equation 59 in Halvorsen (2013). See Halvorsen et al. (2015) for an
#' explanation of the forward selection procedure.
#'
#' When \code{interaction = TRUE}, the forward selection procedure selects a
#' parsimonious group of individual EVs first, and then tests interactions
#' between EVs included in the model afterwards. Therefore, interactions are
#' only explored between terms which are individually explain a significant
#' amount of variation. When \code{interaction = FALSE}, interactions are not
#' considered.
#'
#' Variables should be uniquely named, and the names should not contain spaces
#' or colons.
#'
#' @param rv Response variable vector. The RV should represent
#'   presence/background data, coded as: 1/NA.
#' @param ev Named list of data frames, with each data frame containing 1 or
#'   more DVs for a given EV. E.g. output [[1]] of \code{selectDVforEV}.
#' @param alpha Alpha-level used in F-test comparison of models. Default is
#'   0.01.
#' @param interaction Logical. Allows interaction terms between pairs of EVs.
#'   Default is \code{TRUE}.
#' @param writedir The directory to which Maxent files will be written during
#'   subset selection of DVs. Defaults to the working directory.
#' @param jarpath The pathway to the maxent.jar executable jar file. If
#'   unspecified, the function looks for the file in the writedir.
#' @param trainmax Integer. Maximum number of uninformed background points to be
#'   used to train the models. May be used to reduce computation time for data
#'   sets with very large numbers of points. Default is no maximum.
#'
#' @return List of length two. The first item is a list of data frames, with one
#'   data frame for each \emph{selected} EV. The second item is a data frame
#'   showing the trail of forward selection of individual EVs (and interaction
#'   terms if necessary).
#'
#' @references Halvorsen, R. (2013). A strict maximum likelihood explanation of
#'   MaxEnt, and some implications for distribution modelling. Sommerfeltia, 36,
#'   1-132.
#' @references Halvorsen, R., Mazzoni, S., Bryn, A., & Bakkestuen, V. (2015).
#'   Opportunities for improved distribution modelling practice via a strict
#'   maximum likelihood interpretation of MaxEnt. Ecography, 38(2), 172-183.
#'
#' @export


selectEV <- function(rv, ev, alpha = 0.01, interaction = TRUE, writedir = NULL,
                     jarpath = NULL, trainmax = NULL) {

  .binaryrvcheck(rv)

  if (is.null(writedir)) {
    writedir <- getwd()
  }

  if (is.null(jarpath)) {
    jarpath <- paste(writedir, "\\maxent.jar", sep="")
  }

  if (file.exists(jarpath) == F) {
    stop("maxent.jar file must be present in writedir, or its pathway must be
specified by the jarpath argument. \n ", call. = FALSE)
  }

  dir <- paste(writedir, "\\selectEV", sep="")
  if (file.exists(dir)) {
    stop("The specified writedir already contains a selection of EVs.
Please specify a different writedir. \n ", call. = FALSE)
  } else {
    dir.create(dir)
  }

  if (!is.null(trainmax)) {
    ub <- min(sum(is.na(rv)), trainmax)
    trainindex <- c(which(!is.na(rv)), sample(which(is.na(rv)), ub))
    rv <- rv[trainindex]
    ev <- lapply(ev, function(x) {x[trainindex, , drop = FALSE]})
  }

  message(paste0("Forward selection of ", length(ev), " EVs"))

  result <- altrMaxent:::.parsevs(rv, ev, alpha, interaction, dir, jarpath)
  write.csv(result[[2]], file = paste(dir, "evselection.csv", sep="\\"),
    row.names = FALSE)

  Result <- list(selectedEV = result[[1]], selection = result[[2]])

  return(Result)
}
