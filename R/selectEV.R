#' Select parsimonious set of explanatory variables.
#'
#' \code{selectEV} selects the parsimonious set of explanatory variables (EVs)
#' which best explains variation in a given response variable (RV). Each EV can
#' be represented by 1 or more derived variables (see \code{\link{deriveVars}}
#' and \code{\link{selectDVforEV}}). The function uses a process of forward
#' selection based on comparison of nested models using inference tests. An EV
#' is selected for inclusion when, during nested model comparison, it accounts
#' for a significant amount of remaining variation, under the alpha value
#' specified by the user. See Halvorsen et al. (2015) for a more detailed
#' explanation of the forward selection procedure.
#'
#' The F-test available in \code{selectEV} is calculated using equation 59 in
#' Halvorsen (2013).
#'
#' When \code{interaction = TRUE}, the forward selection procedure selects a
#' parsimonious group of individual EVs first, and then tests interactions
#' between EVs included in the model afterwards. Therefore, interactions are
#' only explored between terms which are individually explain a significant
#' amount of variation. When \code{interaction = FALSE}, interactions are not
#' considered. Practically, interactions between EVs are represented by the
#' products of all combinations of their component DVs (Halvorsen, 2013).
#'
#' The maximum entropy algorithm ("maxent") --- which is implemented in
#' MIAmaxent as an infinitely-weighted logistic regression with presences added
#' to the background --- is conventionally used with presence-only occurrence
#' data. In contrast, standard logistic regression (algorithm = "LR"), is
#' conventionally used with presence-absence occurrence data.
#'
#' Explanatory variables should be uniquely named. Underscores ('_') and colons
#' (':') are reserved to denote derived variables and interaction terms
#' respectively, and \code{selectEV} will replace these --- along with other
#' special characters --- with periods ('.').
#'
#' @param dvdata List containing first the response variable, followed by data
#'   frames of \emph{selected} derived variables for a given explanatory
#'   variable (e.g. the first item in the list returned by
#'   \code{\link{selectDVforEV}}).
#' @param alpha Alpha-level used in F-test comparison of models. Default is
#'   0.01.
#' @param retest Logical. Test variables (or interaction terms) that do not meet
#'   the alpha criterion in a given round in subsequent rounds? Default is
#'   \code{FALSE}.
#' @param interaction Logical. Allow interaction terms between pairs of EVs?
#'   Default is \code{FALSE}.
#' @param formula A model formula (in the form y ~ x + ...) specifying a
#'   starting point for forward model selection. The independent terms in the
#'   formula will be included in the model regardless of explanatory power, and
#'   must be represented in \code{dvdata}, while the remaining explanatory
#'   variables in \code{dvdata} are candidates for selection. The first list
#'   item in \code{dvdata} is still taken as the response variable, regardless
#'   of \code{formula}. Default is \code{NULL}, meaning that forward selection
#'   starts with zero selected variables.
#' @param test Character string matching either "Chisq" or "F" to determine
#'   which inference test is used in nested model comparison. The Chi-squared
#'   test is implemented by stats::anova, while the F-test is implemented as
#'   described in Halvorsen (2013, 2015). Default is "Chisq".
#' @param algorithm Character string matching either "maxent" or "LR", which
#'   determines the type of model used during forward selection. Default is
#'   "maxent".
#' @param write Logical. Write the trail of forward selection to .csv file?
#'   Default is \code{FALSE}.
#' @param dir Directory for file writing if \code{write = TRUE}. Defaults to the
#'   working directory.
#' @param quiet Logical. Suppress progress messages from EV-selection?
#'
#' @return List of 3: \enumerate{ \item dvdata: A list containing first the
#'   response variable, followed by data frames of DVs for each \emph{selected}
#'   EV. \item selection: A data frame showing the trail of forward selection of
#'   individual EVs (and interaction terms if necessary). \item selectedmodel:
#'   the selected model under the given alpha value.}
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
#' # From vignette:
#' grasslandEVselect <- selectEV(grasslandDVselect$dvdata, alpha = 0.001,
#'                               interaction = TRUE)
#' summary(grasslandDVselect$dvdata)
#' length(grasslandDVselect$dvdata[-1])
#' summary(grasslandEVselect$dvdata)
#' length(grasslandEVselect$dvdata[-1])
#' grasslandEVselect$selectedmodel$formula
#' }
#'
#' @export


selectEV <- function(dvdata, alpha = 0.01, retest = FALSE, interaction = FALSE,
                     formula = NULL, test="Chisq", algorithm = "maxent",
                     write = FALSE, dir = NULL, quiet = FALSE) {

  names(dvdata) <- make.names(names(dvdata), allow_ = FALSE)
  stopifnot(class(dvdata)=="list",
            length(dvdata)>1,
            all(lapply(dvdata[-1], class)=="data.frame"))
  .binaryrvcheck(dvdata[[1]])

  if (write == TRUE) {
    if (is.null(dir)) { dir <- getwd() }
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
  }

  if (!is.null(formula)) {
    formula <- stats::as.formula(formula)
    .formulacheck(formula, dvdata)
  }

  if (quiet == F) {
    if (!is.null(formula) && length(labels(stats::terms(formula))) != 0) {
      nterms <- length(labels(stats::terms(formula)))
      message(paste0("Forward selection of ", length(dvdata[-1]) - nterms, " EVs"))
      } else {
      message(paste0("Forward selection of ", length(dvdata[-1]), " EVs"))
      }
  }

  result <- .parsevs(dvdata, alpha, retest, interaction, formula, test,
                     algorithm, quiet)
  if (write == TRUE) {
    utils::write.csv(result[[2]], file = file.path(fdir, "evselection.csv"),
                     row.names = FALSE)
  }

  Result <- list("dvdata"=c(dvdata[1], result[[1]]), "selection"=result[[2]],
                 "selectedmodel"=result[[3]])
  return(Result)
}
