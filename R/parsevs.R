#' Forward selection of EVs
#'
#' @param data List of data frames with response variable (first element) and EVs to
#'   be selected from (remaining elements).
#' @param alpha Alpha level for F-test.
#' @param test Character string matching either "Chisq" or "F".
#' @param interaction Logical. Allows interaction terms.
#' @param formula Model formula specifying a starting point for model selection.
#'

.parsevs <- function(data, alpha, test, interaction, formula) {

  test <- match.arg(test, choices = c("Chisq", "F"))
  dfnames <-  unlist(lapply(data, names))
  df <- data.frame(do.call(cbind, data))
  names(df) <- dfnames
  roundnumber <- 0
  refformula <- stats::formula(paste(names(df)[1], "~ 1"))

  if (!is.null(formula) && length(labels(stats::terms(formula))) != 0) {
    selectedset <- labels(stats::terms(formula))
    remainingset <- names(data)[-1][!(names(data)[-1] %in% selectedset)]
    selectedsetdvs <- unlist(lapply(data[selectedset], names))
    formula <- stats::formula(paste(names(df)[1], "~",
                                    paste(selectedsetdvs, collapse=" + ")))
    ctable <- .compare(list(formula), refformula, df, test=test)
    ctable$variables <- paste(selectedset, collapse=" + ")
    modeltable <- data.frame("round"=roundnumber, ctable)
    refformula <- formula
    message("Round 0 complete.")

    if (length(remainingset) < 1) {
      selectedmodel <- .runIWLR(formula, df)
      return(list(data[selectedset], modeltable, selectedmodel))
    }

  } else {
    selectedset <- character(length=0)
    remainingset <- names(data)[-1]
    modeltable <- data.frame()
  }

  iterationexit <- FALSE
  while (iterationexit == FALSE) {

    roundnumber <- roundnumber + 1
    remainingsetdvs <- lapply(data[remainingset], names)
    formulas <- lapply(remainingsetdvs, function(x) {
      stats::update.formula(refformula,
                            paste("~ . +", paste(x, collapse = " + ")))})
    ctable <- .compare(formulas, refformula, df, test=test)
    if (length(selectedset)==0) {
      variables <- remainingset
    } else { variables <- paste(paste(selectedset, collapse=" + "),
                                remainingset, sep=" + ")}
    ctable$variables <- variables
    ctable <- ctable[order(ctable$P,
                           -ctable[, match(test, names(ctable))]), ]
    modeltable <- rbind(modeltable, data.frame("round"=roundnumber, ctable),
                     make.row.names = FALSE)

    if (!is.na(ctable$P[1]) && ctable$P[1] < alpha) {
      selectedset <- unlist(strsplit(ctable$variables[1], split=" + ",
                                     fixed=TRUE))
      selectedsetdvs <- unlist(lapply(data[selectedset], names))
      refformula <- stats::formula(paste(names(df)[1], "~",
                                         paste(selectedsetdvs, collapse=" + ")))
      testedEVs <- sapply(strsplit(ctable$variables[seq(nrow(ctable))[-1]],
                                   split=" + ", fixed=TRUE),
                          function(x) {x[length(selectedset)]})
      remainingset <- testedEVs[ctable$P[seq(nrow(ctable))[-1]] < alpha]
      remainingset <- remainingset[!is.na(remainingset)]
    }

    message(paste("Round", roundnumber, "complete."))

    if (nrow(ctable) == 1 ||
        ctable$P[1] > alpha ||
        (all(ctable$P[seq(nrow(ctable))[-1]] >= alpha |
             is.na(ctable$P[seq(nrow(ctable))[-1]])))) {
      iterationexit <- TRUE
    }
  }

  if (interaction == FALSE || length(selectedset) < 2) {
    selectedmodel <- .runIWLR(refformula, df)
    return(list(data[selectedset], modeltable, selectedmodel))
  }

  message(paste("Forward selection of interaction terms between",
    length(selectedset), "EVs"))

  remainingset <- apply(utils::combn(selectedset, 2), 2,
                        function(x) {paste(x, collapse=":")})
  maineffectdvs <- lapply(data, names)
  interactiondvs <- lapply(strsplit(remainingset, ":", fixed=TRUE),
                            function(x, data) {
                              dvs1 <- names(data[[x[1]]])
                              dvs2 <- names(data[[x[2]]])
                              apply(expand.grid(dvs1, dvs2), 1, function(y) {
                                paste(y, collapse=":")})
                            }, data=data)
  names(interactiondvs) <- remainingset
  dvs <- c(maineffectdvs, interactiondvs)

  iterationexit <- FALSE
  while (iterationexit == FALSE) {

    roundnumber <- roundnumber + 1
    remainingsetdvs <- dvs[remainingset]
    formulas <- lapply(remainingsetdvs, function(x) {
      stats::update.formula(refformula,
                            paste("~ . +", paste(x, collapse = " + ")))})
    ctable <- .compare(formulas, refformula, df, test=test)
    variables <- paste(paste(selectedset, collapse=" + "),
                                remainingset, sep=" + ")
    ctable$variables <- variables
    ctable <- ctable[order(ctable$P,
                           -ctable[, match(test, names(ctable))]), ]
    modeltable <- rbind(modeltable, data.frame("round"=roundnumber, ctable),
                        make.row.names = FALSE)

    if (!is.na(ctable$P[1]) && ctable$P[1] < alpha) {
      selectedset <- unlist(strsplit(ctable$variables[1], split=" + ",
                                     fixed=TRUE))
      selectedsetdvs <- unlist(dvs[selectedset])
      refformula <- stats::formula(paste(names(df)[1], "~",
                                         paste(selectedsetdvs, collapse=" + ")))
      testedEVs <- sapply(strsplit(ctable$variables[seq(nrow(ctable))[-1]],
                                   split=" + ", fixed=TRUE),
                          function(x) {x[length(selectedset)]})
      remainingset <- testedEVs[ctable$P[seq(nrow(ctable))[-1]] < alpha]
      remainingset <- remainingset[!is.na(remainingset)]
    }

    message(paste("Round", roundnumber, "complete."))

    if (nrow(ctable) == 1 ||
        ctable$P[1] > alpha ||
        (all(ctable$P[seq(nrow(ctable))[-1]] >= alpha |
             is.na(ctable$P[seq(nrow(ctable))[-1]])))) {
      iterationexit <- TRUE
    }
  }

  selectedmodel <- .runIWLR(refformula, df)
  return(list(data[selectedset], modeltable, selectedmodel))

}