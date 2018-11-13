#' Trains a model containing the explanatory variables specified.
#'
#' \code{chooseModel} trains a model based on the formula provided. The formula
#' specifies which explanatory variables (EVs) --- and potentially first-order
#' interactions between these --- should be included in the model. Each EV can
#' be represented by 1 or more derived variables (see \code{\link{deriveVars}}).
#' The function may be employed to choose a model from the selection pathway of
#' \code{\link{selectEV}} other than the model selected under the provided alpha
#' value.
#'
#' Explanatory variables should be uniquely named. Underscores ('_') and colons
#' (':') are reserved to denote derived variables and interaction terms
#' respectively, and \code{chooseModel} will replace these --- along with other
#' special characters --- with periods ('.').
#'
#' @param dvdata A list containing first the response variable, followed by data
#'   frames of \emph{selected} derived variables for a given explanatory
#'   variable (e.g. the first item in the list returned by
#'   \code{\link{selectDVforEV}}).
#' @param formula A model formula (in the form y ~ x + ...) specifying the
#'   independent terms (EVs) to be included in the model. The item in
#'   \code{dvdata} is still taken as the response variable, regardless of
#'   \code{formula}.
#' @param algorithm Character string matching either "maxent" or "LR", which
#'   determines the type of model built. Default is "maxent".
#'
#' @examples
#' \dontrun{
#' # From vignette:
#' grasslandmodel <- chooseModel(grasslandDVselect$dvdata,
#'                               formula("~ pr.bygall + geoberg + lcucor1 +
#'                                       tertpi09 + geolmja1"))
#' }
#'
#' @export


chooseModel <- function(dvdata, formula, algorithm = "maxent") {

  names(dvdata) <- make.names(names(dvdata), allow_ = FALSE)
  .binaryrvcheck(dvdata[[1]])
  algorithm <- match.arg(algorithm, choices = c("maxent", "LR"))

  formula <- stats::as.formula(formula)
  trms <- labels(stats::terms(formula))
  firstordertrms <- trms[attr(stats::terms(formula), "order")==1]
  secondordertrms <- trms[attr(stats::terms(formula), "order")==2]
  for (i in firstordertrms) {
    if (sum(names(dvdata) == i) != 1) {
      stop(paste(i, "must be represented in 'dvdata' (exactly once)"),
           call. = FALSE) }
  }

  dvdata[[1]] <- data.frame("RV"=dvdata[[1]])
  dvdata <- c(dvdata[1], dvdata[-1][names(dvdata[-1]) %in% firstordertrms])
  dfnames <-  unlist(lapply(dvdata, names))
  df <- data.frame(do.call(cbind, dvdata))
  names(df) <- dfnames

  maineffectdvs <- lapply(dvdata[-1], names)
  interactiondvs <- lapply(strsplit(secondordertrms, ":", fixed=TRUE),
                           function(x, data) {
                             dvs1 <- names(data[[x[1]]])
                             dvs2 <- names(data[[x[2]]])
                             apply(expand.grid(dvs1, dvs2), 1, function(y) {
                               paste(y, collapse=":")})
                           }, data=dvdata)
  names(interactiondvs) <- secondordertrms
  dvs <- c(maineffectdvs, interactiondvs)

  selectedsetdvs <- unlist(dvs[trms])
  formula <- stats::formula(paste(names(df)[1], "~",
                                  paste(selectedsetdvs, collapse=" + ")))
  if (algorithm == "maxent") {
    model <- .runIWLR(formula, df)
  } else {
    model <- .runLR(formula, df)
  }
  return(model)
}
