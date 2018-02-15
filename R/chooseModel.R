#' Trains a model containing the explanatory variables specified.
#'
#' \code{chooseModel} trains a model based on the formula provided. The formula
#' specifies which explanatory variables (EVs) -- and potentially first-order
#' interactions between these -- should be included in the model. Each EV can be
#' represented by 1 or more derived variables (see \code{\link{deriveVars}}).
#' The function may be employed to choose a model from the selection pathway of
#' \code{\link{selectEV}} other than the model selected under the provided alpha
#' value.
#'
#' Explanatory variables should be uniquely named, and the names must not
#' contain spaces, underscores, or colons. Underscores and colons are reserved
#' to denote derived variables and interaction terms repectively.
#'
#' @param data Data frame containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent presence/background data, coded as: 1/NA. See
#'   \code{\link{readData}}.
#' @param dvdata List of data frames, with each data frame containing selected
#'   derived variables for a given explanatory variable (e.g. the first item in
#'   the list returned by \code{\link{selectDVforEV}}).
#' @param formula A model formula (in the form y ~ x + ...) specifying the
#'   independent terms to be included in the model. The first column in
#'   \code{data} is still taken as the response variable, regardless of
#'   \code{formula}.
#'
#' @examples
#'
#' @export


chooseModel <- function(data, dvdata, formula) {

  .binaryrvcheck(data[, 1])

  formula <- stats::as.formula(formula)
  terms <- labels(stats::terms(formula))
  firstorderterms <- terms[attr(stats::terms(formula), "order")==1]
  secondorderterms <- terms[attr(stats::terms(formula), "order")==2]
  for (i in firstorderterms) {
    if (sum(names(dvdata) == i) != 1) {
      stop(paste(i, "must be represented in 'dvdata' (exactly once)"),
           call. = FALSE) }
  }

  names(data)[1] <- make.names(names(data)[1], allow_ = FALSE)
  datalist <- c(list("RV"=data[, 1, drop=FALSE]),
                dvdata[names(dvdata) %in% firstorderterms])
  dfnames <-  unlist(lapply(datalist, names))
  df <- data.frame(do.call(cbind, datalist))
  names(df) <- dfnames

  maineffectdvs <- lapply(datalist, names)
  interactiondvs <- lapply(strsplit(secondorderterms, ":", fixed=TRUE),
                           function(x, data) {
                             dvs1 <- names(data[[x[1]]])
                             dvs2 <- names(data[[x[2]]])
                             apply(expand.grid(dvs1, dvs2), 1, function(y) {
                               paste(y, collapse=":")})
                           }, data=datalist)
  names(interactiondvs) <- secondorderterms
  dvs <- c(maineffectdvs, interactiondvs)

  selectedsetdvs <- unlist(dvs[terms])
  formula <- stats::formula(paste(names(df)[1], "~",
                                  paste(selectedsetdvs, collapse=" + ")))
  model <- .runIWLR(formula, df)
  return(model)
}
