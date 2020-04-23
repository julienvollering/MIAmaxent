#' Select parsimonious sets of derived variables.
#'
#' For each explanatory variable (EV), \code{selectDVforEV} selects the
#' parsimonious set of derived variables (DV) which best explains variation in a
#' given response variable. The function uses a process of forward selection
#' based on comparison of nested models using inference tests. A DV is selected
#' for inclusion when, during nested model comparison, it accounts for a
#' significant amount of remaining variation, under the alpha value specified by
#' the user. See Halvorsen et al. (2015) for a more detailed explanation of the
#' forward selection procedure.
#'
#' The F-test available in \code{selectDVforEV} is calculated using equation 59
#' in Halvorsen (2013).
#'
#' If using binary-type derived variables from \code{\link{deriveVars}}, be
#' aware that a model including all of these DVs will be considered equal to the
#' the closest nested model, due to perfect multicollinearity (i.e. the dummy
#' variable trap).
#'
#' The maximum entropy algorithm ("maxent") --- which is implemented in
#' MIAmaxent as an infinitely-weighted logistic regression with presences added
#' to the background --- is conventionally used with presence-only occurrence
#' data. In contrast, standard logistic regression (algorithm = "LR"), is
#' conventionally used with presence-absence occurrence data.
#'
#' Explanatory variables should be uniquely named. Underscores ('_') and colons
#' (':') are reserved to denote derived variables and interaction terms
#' respectively, and \code{selectDVforEV} will replace these --- along with
#' other special characters --- with periods ('.').
#'
#' @param dvdata List containing first the response variable, followed by data
#'   frames of derived variables produced for each explanatory variable (e.g.
#'   the first item in the list returned by \code{\link{deriveVars}}).
#' @param alpha Alpha-level used for inference testing in nested model
#'   comparison. Default is 0.01.
#' @param retest Logical. Test variables that do not meet the alpha criterion
#'   in a given round in subsequent rounds? Default is \code{FALSE}.
#' @param test Character string matching either "Chisq" or "F" to determine
#'   which inference test is used in nested model comparison. The Chi-squared
#'   test is implemented by stats::anova, while the F-test is implemented as
#'   described in Halvorsen (2013, 2015). Default is "Chisq".
#' @param algorithm Character string matching either "maxent" or "LR", which
#'   determines the type of model used during forward selection. Default is
#'   "maxent".
#' @param write Logical. Write the trail of forward selection for each EV to
#'   .csv file? Default is \code{FALSE}.
#' @param dir Directory for file writing if \code{write = TRUE}. Defaults to the
#'   working directory.
#' @param quiet Suppress progress bar?
#'
#' @return List of 2: \enumerate{ \item dvdata: A list containing first the
#'   response variable, followed by data frames of \emph{selected} DVs for each
#'   EV. EVs with zero selected DVs are dropped. This item is recommended as
#'   input for \code{dvdata} in \code{\link{selectEV}}. \item selection: A list
#'   of data frames, where each data frame shows the trail of forward selection
#'   of DVs for a given EV. }
#'
#' @references Halvorsen, R. (2013). A strict maximum likelihood explanation of
#'   MaxEnt, and some implications for distribution modelling. Sommerfeltia, 36,
#'   1-132.
#' @references Halvorsen, R., Mazzoni, S., Bryn, A., & Bakkestuen, V. (2015).
#'   Opportunities for improved distribution modelling practice via a strict
#'   maximum likelihood interpretation of MaxEnt. Ecography, 38(2), 172-183.
#'
#' @examples
#' toydata_seldvs <- selectDVforEV(toydata_dvs$dvdata, alpha = 0.4)
#'
#' \dontrun{
#' # From vignette:
#' grasslandDVselect <- selectDVforEV(grasslandDVs$dvdata, alpha = 0.001)
#' summary(grasslandDVs$dvdata)
#' sum(sapply(grasslandDVs$dvdata[-1], length))
#' summary(grasslandDVselect$dvdata)
#' sum(sapply(grasslandDVselect$dvdata[-1], length))
#' grasslandDVselect$selection$terdem
#' }
#'
#' @export


selectDVforEV <- function(dvdata, alpha = 0.01, retest = FALSE, test = "Chisq",
                          algorithm = "maxent", write = FALSE, dir = NULL,
                          quiet = FALSE) {

  names(dvdata) <- make.names(names(dvdata), allow_ = FALSE)
  stopifnot(class(dvdata)=="list",
            length(dvdata)>1,
            all(lapply(dvdata[-1], class)=="data.frame"))
  rv <- dvdata[[1]]
  .binaryrvcheck(rv)
  evdv <- dvdata[-1]


  if (write == TRUE) {
    if (is.null(dir)) { dir <- getwd() }
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
  }

  sdvdata <- list()
  trail <- list()

  if (quiet == FALSE) {
    message(paste0("Forward selection of DVs for ", length(evdv), " EVs"))
    pb <- utils::txtProgressBar(min = 0, max = length(evdv), style = 3)
  }

  for (i in 1:length(evdv)) {
    evname <- names(evdv)[i]
    df <- data.frame("RV"=rv, evdv[[i]])
    result <- .parsdvs(df, alpha, retest, test, algorithm)
    if (write == TRUE) {
      utils::write.csv(result[[2]],
                       file=file.path(fdir, paste0(evname, "_dvselection.csv")),
                       row.names = FALSE)
    }
    sdvdata[[i]] <- result[[1]]
    trail[[i]] <- result[[2]]
    if (quiet == FALSE) {utils::setTxtProgressBar(pb, i)}
  }
  if (quiet == FALSE) {close(pb)}

  names(sdvdata) <- names(evdv)
  names(trail) <- names(evdv)
  sdvdata <- sdvdata[sapply(sdvdata, function(x) {dim(x)[2] != 0})]

  Result <- list(dvdata = c(list("RV"=rv), sdvdata), selection = trail)

  return(Result)
}
