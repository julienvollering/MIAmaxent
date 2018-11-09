#' Forward selection of DVs
#'
#' @param df Data frame with response variable in first column and DVs to be
#'   selected from in subsequent columns.
#' @param alpha Alpha level for inference test.
#' @param test Character string matching either "Chisq" or "F".
#' @param algorithm Character string matching either "maxent" or "LR".
#'
#' @keywords internal
#' @noRd

.parsdvs <- function(df, alpha, test, algorithm) {

  selectedset <- character(length=0)
  remainingset <- names(df)[-1]
  evtable <- data.frame()
  roundnumber <- 0
  refformula <- stats::formula(paste(names(df)[1], "~ 1"))

  iterationexit <- FALSE
  while (iterationexit == FALSE) {

    roundnumber <- roundnumber + 1
    formulas <- lapply(remainingset, function(x) {
      stats::update.formula(refformula,  paste("~ . +", x))})
    ctable <- .compare(formulas, refformula, df, test, algorithm)
    ctable <- ctable[order(ctable$P,
                           -ctable[, match(test, names(ctable))]), ]
    evtable <- rbind(evtable, data.frame("round"=roundnumber, ctable),
                     make.row.names = FALSE)

    if (!is.na(ctable$P[1]) && ctable$P[1] < alpha) {
      selectedset <- unlist(strsplit(ctable$variables[1], split=" + ",
                                     fixed=TRUE))
      refformula <- stats::formula(paste(names(df)[1], "~",
                                         paste(selectedset, collapse=" + ")))
      testedDVs <- sapply(strsplit(ctable$variables[seq(nrow(ctable))[-1]],
                                   split=" + ", fixed=TRUE),
                          function(x) {x[length(selectedset)]})
      remainingset <- testedDVs[ctable$P[seq(nrow(ctable))[-1]] < alpha]
      remainingset <- remainingset[!is.na(remainingset)]
    }

    if (nrow(ctable) == 1 ||
        ctable$P[1] > alpha ||
        (all(ctable$P[seq(nrow(ctable))[-1]] >= alpha |
             is.na(ctable$P[seq(nrow(ctable))[-1]])))) {
      iterationexit <- TRUE
    }
  }

  return(list(df[,selectedset, drop = FALSE], evtable))
}