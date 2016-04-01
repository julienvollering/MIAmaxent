#' Select parsimonious sets of derived variables.
#'
#' For each explanatory variable, \code{selectDV} selects the parsimonious set
#' of derived variables which best explains variation in a given response
#' variable. The function uses a process of forward selection based on
#' comparison of nested models by the F-test.
#'
#' If the derived variables were created using \code{\link{deriveVars}}, the
#' same response variable should be used in \code{selectDV}, as the deviation
#' and spline transformations produced by \code{deriveVars} are RV-specific.
#'
#' @param rv Response variable vector. The RV should represent
#'   presence/background data, coded as: 1/NA.
#' @param dv List of data frames, with each data frame containing DVs for a
#'   given EV. E.g. output of \code{deriveVars}.
#' @param alpha Alpha-level used in F-test comparison of models.
#' @param writedir The directory to which Maxent files will be written during
#'   subset selection of DVs. Defaults to the working directory.
#' @param jarpath The pathway to the maxent.jar executable jar file. If
#'   unspecified, the function looks for the file in the writedir.
#'
#' @return List of data frames, with each data frame containing \emph{selected}
#'   DVs for a given EV. Dimensions: [[EV]] [N, <=DV].
#'
#' @references Halvorsen, R., Mazzoni, S., Bryn, A., & Bakkestuen, V. (2015).
#'   Opportunities for improved distribution modelling practice via a strict
#'   maximum likelihood interpretation of MaxEnt. Ecography, 38(2), 172-183.
#'
#' @export


selectDV <- function(rv, dv, alpha = 0.01, writedir = NULL, jarpath = NULL) {

  altrMaxent:::.binaryrvcheck(rv)

  if (is.null(writedir)) {
    writedir <- getwd()
  }

  if (is.null(jarpath)) {
    jarpath <- paste(writedir, "\\maxent.jar", sep="")
  }

  if (file.exists(jarpath) == F) {
    stop("maxent.jar file must be present in writedir, or its pathway must be
specified by the jarpath parameter. \n ")
  }

  dir <- paste(writedir, "\\selectDV", sep="")
  if (file.exists(dir)) {
    stop("The specified writedir already contains a selection of DVs.
Please specify a different writedir. \n ")
  } else {
    dir.create(dir)
  }

  EVDV <- list()
  trail <- list()
  for (i in 1:length(dv)) {
    evname <- names(dv)[i]
    evdir <- paste(dir, "\\", evname, sep="")
    dir.create(evdir)
    df <- dv[[i]]
    result <- altrMaxent:::.parsdvs(rv, df, alpha, evdir, jarpath)
    EVDV[[i]] <- result[[1]]
    trail[[i]] <- result[[2]]
  }
  names(EVDV) <- names(dv)
  names(trail) <- names(dv)

  Result <- list(selectedDV = EVDV, selection = trail)

  return(Result)
}
