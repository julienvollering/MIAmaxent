#' Select parsimonious sets of derived variables.
#'
#' For each explanatory variable (EV), \code{selectDVforEV} selects the
#' parsimonious set of derived variables (DV) which best explains variation in a
#' given response variable. The function uses a process of forward selection
#' based on comparison of nested models using inference tests. A DV is selected
#' for inclusion when, during nested model comparison, it accounts for a
#' significant amount of remaining variation, under the alpha value specified by
#' the user.
#'
#' The F-statistic that \code{selectDVforEV} uses for nested model comparison is
#' calculated using equation 59 in Halvorsen (2013). See Halvorsen et al. (2015)
#' for a more detailed explanation of the forward selection procedure.
#'
#' If the derived variables were created using \code{\link{deriveVars}}, the
#' same response variable should be used in \code{selectDVforEV}, because the
#' deviation and spline transformations produced by \code{deriveVars} are
#' RV-specific.
#'
#' If using binary-type derived variables from \code{\link{deriveVars}}, be
#' aware that a model including all of these DVs will be considered equal to the
#' the closest nested model, due to perfect multicollinearity (i.e. the dummy
#' variable trap).
#'
#' Explanatory variables should be uniquely named, and the names must not
#' contain spaces, underscores, or colons. Underscores and colons are reserved
#' to denote derived variables and interaction terms repectively.
#'
#' @param data Data frame containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent presence/background data, coded as: 1/NA. See
#'   \code{\link{readData}}.
#' @param dvdata List of data frames, with each data frame containing derived
#'   variables for a given explanatory variable (e.g. the first item in the list
#'   returned by \code{\link{deriveVars}}).
#' @param alpha Alpha-level used for inference testing in nested model
#'   comparison. Default is 0.01.
#' @param test Character string matching either "Chisq" or "F" to determine
#'   which inference test is used in nested model comparison. The Chi-squared
#'   test is implemented as in stats::anova, while the F-test is implemented as
#'   described in Halvorsen (2013, 2015). Default is "Chisq".
#' @param dir Directory to which files will be written during subset selection
#'   of derived variables. Defaults to the working directory.
#'
#' @return List of 2: \enumerate{ \item A list of data frames, with each data
#'   frame containing \emph{selected} DVs for a given EV. This item is
#'   recommended as input for \code{dvdata} in \code{\link{selectEV}}. \item A
#'   list of data frames, where each data frame shows the trail of forward
#'   selection of DVs for a given EV. }
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
#' selecteddvs <- selectDVforEV(dat, deriveddat, alpha = 0.0001,
#'    dir = "D:/path/to/modeling/directory")
#'
#' # From vignette:
#' grasslandDVselect <- selectDVforEV(grasslandPO, grasslandDVs[[1]], alpha = 0.001)
#' summary(grasslandDVs$EVDV)
#' sum(sapply(grasslandDVs$EVDV, length))
#' summary(grasslandDVselect$selectedDV)
#' sum(sapply(grasslandDVselect$selectedDV, length))
#' }
#'
#' @export


selectDVforEV <- function(data, dvdata, alpha = 0.01, test="Chisq",
                          dir = NULL) {

  rv <- data[, 1]
  .binaryrvcheck(rv)

  if (is.null(dir)) { dir <- getwd()}

  fdir <- file.path(dir, "selectDVforEV")
  if (file.exists(fdir)) {
    yn <- readline("The specified dir already contains a selectDVforEV result. Overwrite this result? (y/n) ")
    if (yn == "y") {
      unlink(fdir, recursive = TRUE)
      Sys.sleep(1)
    }
    if (yn != "y") { return(message("Overwrite declined")) }
  }
  dir.create(fdir, recursive = TRUE)

  EVDV <- list()
  trail <- list()

  message(paste0("Forward selection of DVs for ", length(dvdata), " EVs"))
  pb <- utils::txtProgressBar(min = 0, max = length(dvdata), style = 3)

  for (i in 1:length(dvdata)) {
    evname <- names(dvdata)[i]
    df <- data.frame("RV"=rv, dvdata[[i]])
    result <- .parsdvs(df, alpha, test=test)
    utils::write.csv(result[[2]],
                     file=file.path(fdir, paste0(evname, "_dvselection.csv")),
                     row.names = FALSE)
    EVDV[[i]] <- result[[1]]
    trail[[i]] <- result[[2]]
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  names(EVDV) <- names(dvdata)
  names(trail) <- names(dvdata)
  EVDV <- EVDV[sapply(EVDV, function(x) {dim(x)[2] != 0})]

  Result <- list(selectedDV = EVDV, selection = trail)

  return(Result)
}
