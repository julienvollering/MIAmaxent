#' Select parsimonious sets of derived variables.
#'
#' For each explanatory variable, \code{selectDV} selects the set of derived
#' variables which best explains variation in a given response variable.
#'
#' If the derived variables were created using \code{\link{deriveVars}}, the
#' same response variable should be used in \code{selectDV}, as the deviation
#' and spline transformations produced by \code{deriveVars} are RV-specific.
#'
#' @param dv List of data frames, with each data frame containing DVs for a
#'   given EV.
#' @param rv Response variable vector. The RV should represent
#'   presence/background data, coded as: 1/NA.
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


selectDV <- function(dv, rv, writedir = NULL, jarpath = NULL) {

  if (any(c("HF", "HR", "T") %in% transformtype) && allsplines == F) {
    altrMaxent:::.binaryrvcheck(data[,1])

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

    dir <- paste(writedir, "\\deriveVars", sep="")
    if (file.exists(dir)) {
      stop("The specified writedir already contains a selection of spline DVs.
Please specify a different writedir. \n ")
    } else {
      dir.create(dir)
    }
  }

  EVDV <- list()
  for (i in 2:ncol(data)) {
    df <- data[,c(1,i)]
    EVDV[[i-1]] <- altrMaxent:::.dvfromev(df, transformtype, allsplines,
      dir, jarpath)
    names(EVDV)[i-1] <- colnames(data)[i]
  }

  return(EVDV)
}
