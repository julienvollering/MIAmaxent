#' Select parsimonious sets of derived variables.
#'
#' For each explanatory variable (EV), \code{selectDV} selects the parsimonious
#' set of derived variables (DV) which best explains variation in a given
#' response variable. The function uses a process of forward selection based on
#' comparison of nested models by the F-test, where the F-statistic is
#' calculated using equation 58 in Halvorsen (2013). See Halvorsen et al. (2015)
#' for an explanation of the forward selection procedure.
#'
#' If the derived variables were created using \code{\link{deriveVars}}, the
#' same response variable should be used in \code{selectDV}, as the deviation
#' and spline transformations produced by \code{deriveVars} are RV-specific.
#'
#' DVs must be uniquely named, and the names must not contain spaces.
#'
#' @param rv Response variable vector. The RV should represent
#'   presence/background data, coded as: 1/NA.
#' @param dv List of data frames, with each data frame containing DVs for a
#'   given EV. E.g. output of \code{deriveVars}.
#' @param alpha Alpha-level used in F-test comparison of models. Default is
#'   0.01.
#' @param writedir The directory to which Maxent files will be written during
#'   subset selection of DVs. Defaults to the working directory.
#' @param jarpath The pathway to the maxent.jar executable jar file. If
#'   unspecified, the function looks for the file in the writedir.
#'
#' @return List of length two. The first item is a list of data frames, with
#'   each data frame containing \emph{selected} DVs for a given EV. The second
#'   item is also a list of data frames, where each data frame shows the trail
#'   of forward selection for a given EV.
#'
#' @references Halvorsen, R. (2013). A strict maximum likelihood explanation of
#'   MaxEnt, and some implications for distribution modelling. Sommerfeltia, 36,
#'   1-132.
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

  message(paste0("Forward selection of DVs for ", length(dv), " EVs"))
  pb <- txtProgressBar(min = 0, max = length(dv), style = 3)

  for (i in 1:length(dv)) {
    evname <- names(dv)[i]
    evdir <- paste(dir, "\\", evname, sep="")
    dir.create(evdir)
    df <- dv[[i]]
    result <- altrMaxent:::.parsdvs(rv, df, alpha, evdir, jarpath)
    write.csv(result[[2]], file = paste(evdir, "dvselection.csv", sep="\\"),
      row.names = FALSE)
    EVDV[[i]] <- result[[1]]
    trail[[i]] <- result[[2]]
    setTxtProgressBar(pb, i)
  }
  names(EVDV) <- names(dv)
  names(trail) <- names(dv)

  Result <- list(selectedDV = EVDV, selection = trail)

  return(Result)
}
