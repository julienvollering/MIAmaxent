#' Select parsimonious sets of derived variables.
#'
#' For each explanatory variable (EV), \code{selectDVforEV} selects the
#' parsimonious set of derived variables (DV) which best explains variation in a
#' given response variable. The function uses a process of forward selection
#' based on comparison of nested models by the F-test, where the F-statistic is
#' calculated using equation 59 in Halvorsen (2013). See Halvorsen et al. (2015)
#' for an explanation of the forward selection procedure.
#'
#' If the derived variables were created using \code{\link{deriveVars}}, the
#' same response variable should be used in \code{selectDVforEV}, as the
#' deviation and spline transformations produced by \code{deriveVars} are
#' RV-specific.
#'
#' Variables should be uniquely named, and the names should not contain spaces
#' or colons.
#'
#' @param data Data frame containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent presence/background data, coded as: 1/NA.
#' @param dvdata List of data frames, with each data frame containing derived
#'   variables for a given explanatory variable (e.g. the first item in the list
#'   returned by \code{\link{deriveVars}}: \code{deriveVars(...)[[1]]})
#' @param alpha Alpha-level used in F-test comparison of models. Default is
#'   0.01.
#' @param dir Directory to which files will be written during subset selection
#'   of derived variables. Defaults to the working directory.
#' @param jar Pathway to the 'maxent.jar' executable jar file. If unspecified,
#'   the function looks for the file in \code{dir}.
#' @param trainmax Integer. Maximum number of uninformed background points to be
#'   used to train the models. May be used to reduce computation time for data
#'   sets with very large numbers of points. Default is no maximum.
#'
#' @return List of 2: \enumerate{ \item A list of data frames, with each data
#'   frame containing \emph{selected} DVs for a given EV. \item A list of data
#'   frames, where each data frame shows the trail of forward selection for a
#'   given EV.}
#'
#' @references Halvorsen, R. (2013). A strict maximum likelihood explanation of
#'   MaxEnt, and some implications for distribution modelling. Sommerfeltia, 36,
#'   1-132.
#' @references Halvorsen, R., Mazzoni, S., Bryn, A., & Bakkestuen, V. (2015).
#'   Opportunities for improved distribution modelling practice via a strict
#'   maximum likelihood interpretation of MaxEnt. Ecography, 38(2), 172-183.
#'
#' @export


selectDVforEV <- function(data, dvdata, alpha = 0.01, dir = NULL, jar = NULL,
                          trainmax = NULL) {

  rv <- data[, 1]
  .binaryrvcheck(rv)

  if (is.null(dir)) { dir <- getwd()}
  jar <- .jar.check(dir, jar)

  fdir <- file.path(dir, "selectDVforEV")
  if (file.exists(fdir)) {
    stop("The specified dir already contains a selection of DVs.
Please specify a different dir. \n ")
  } else {
    dir.create(fdir)
  }

  if (!is.null(trainmax)) {
    ub <- min(sum(is.na(rv)), trainmax)
    trainindex <- c(which(!is.na(rv)), sample(which(is.na(rv)), ub))
    rv <- rv[trainindex]
    dvdata <- lapply(dvdata, function(x) {x[trainindex, , drop = FALSE]})
  }

  EVDV <- list()
  trail <- list()

  message(paste0("Forward selection of DVs for ", length(dvdata), " EVs"))
  pb <- txtProgressBar(min = 0, max = length(dvdata), style = 3)

  for (i in 1:length(dvdata)) {
    evname <- names(dvdata)[i]
    evdir <- .dirpath.create(fdir, evname)
    df <- dvdata[[i]]
    result <- .parsdvs(rv, df, alpha, evdir, jar)
    write.csv(result[[2]], file = file.path(evdir, "dvselection.csv"),
      row.names = FALSE)
    EVDV[[i]] <- result[[1]]
    trail[[i]] <- result[[2]]
    setTxtProgressBar(pb, i)
  }
  names(EVDV) <- names(dvdata)
  names(trail) <- names(dvdata)
  EVDV <- EVDV[sapply(EVDV, function(x) {dim(x)[2] != 0})]

  Result <- list(selectedDV = EVDV, selection = trail)

  return(Result)
}
