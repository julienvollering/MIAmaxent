#' Forward selection of EVs
#'
#' @param dvdata List with response variable (vector or single-column data
#'   frame) followed by EVs to be selected from (data frames).
#' @param alpha Alpha level for F-test.
#' @param retest Test rejected variables in subsequent rounds?
#' @param test Character string matching either "Chisq" or "F".
#' @param algorithm Character string matching either "maxent" or "LR".
#' @param interaction Logical. Allows interaction terms.
#' @param formula Model formula specifying a starting point for model selection.
#' @param quiet Logical. Suppress progress messages?
#'
#' @keywords internal
#' @noRd

.parsevs <- function(dvdata, alpha, retest, interaction, formula, test,
                     algorithm, quiet) {

  dvdata[[1]] <- data.frame("RV"=dvdata[[1]])
  test <- match.arg(test, choices = c("Chisq", "F"))
  algorithm <- match.arg(algorithm, choices = c("maxent", "LR"))
  dfnames <-  unlist(lapply(dvdata, names))
  df <- data.frame(do.call(cbind, dvdata))
  names(df) <- dfnames
  roundnumber <- 0
  refformula <- stats::formula(paste(names(df)[1], "~ 1"))

  if (!is.null(formula) && length(labels(stats::terms(formula))) != 0) {
    selectedset <- labels(stats::terms(formula))
    remainingset <- names(dvdata)[-1][!(names(dvdata)[-1] %in% selectedset)]
    selectedsetdvs <- unlist(lapply(dvdata[selectedset], names))
    formula <- stats::formula(paste(names(df)[1], "~",
                                    paste(selectedsetdvs, collapse=" + ")))
    ctable <- .compare(list(formula), refformula, df, test, algorithm)
    ctable$variables <- paste(selectedset, collapse=" + ")
    modeltable <- data.frame("round"=roundnumber, ctable)
    refformula <- formula
    if (quiet == FALSE) {message("Round 0 complete.")}

    if (length(remainingset) < 1) {
      if (algorithm == "maxent") {
        selectedmodel <- .runIWLR(formula, df)
      } else {
        selectedmodel <- .runLR(formula, df)
      }
      return(list(dvdata[selectedset], modeltable, selectedmodel))
    }

  } else {
    selectedset <- character(length=0)
    remainingset <- names(dvdata)[-1]
    modeltable <- data.frame()
  }

  iterationexit <- FALSE
  while (iterationexit == FALSE) {

    roundnumber <- roundnumber + 1
    remainingsetdvs <- lapply(dvdata[remainingset], names)
    formulas <- lapply(remainingsetdvs, function(x) {
      stats::update.formula(refformula,
                            paste("~ . +", paste(x, collapse = " + ")))})
    ctable <- .compare(formulas, refformula, df, test, algorithm)
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
      selectedsetdvs <- unlist(lapply(dvdata[selectedset], names))
      refformula <- stats::formula(paste(names(df)[1], "~",
                                         paste(selectedsetdvs, collapse=" + ")))
      testedEVs <- sapply(strsplit(ctable$variables[seq(nrow(ctable))[-1]],
                                   split=" + ", fixed=TRUE),
                          function(x) {x[length(selectedset)]})
      if (retest) {
        remainingset <- testedEVs
      } else {
        remainingset <- testedEVs[ctable$P[seq(nrow(ctable))[-1]] < alpha]
      }
      remainingset <- remainingset[!is.na(remainingset)]
    }

    if (quiet == FALSE) {message(paste("Round", roundnumber, "complete."))}

    if (nrow(ctable) == 1 || ctable$P[1] > alpha) { iterationexit <- TRUE }
    if (!iterationexit && !retest) {
      if (all(ctable$P[seq(nrow(ctable))[-1]] >= alpha |
              is.na(ctable$P[seq(nrow(ctable))[-1]]))) { iterationexit <- TRUE }
    }

  }

  if (interaction == FALSE || length(selectedset) < 2) {
    if (algorithm == "maxent") {
      selectedmodel <- .runIWLR(refformula, df)
    } else {
      selectedmodel <- .runLR(refformula, df)
    }
    return(list(dvdata[selectedset], modeltable, selectedmodel))
  }

  if (quiet == FALSE) {
    message(paste("Forward selection of interaction terms between",
                  length(selectedset), "EVs"))
  }

  remainingset <- apply(utils::combn(selectedset, 2), 2,
                        function(x) {paste(x, collapse=":")})
  maineffectdvs <- lapply(dvdata[-1], names)
  interactiondvs <- lapply(strsplit(remainingset, ":", fixed=TRUE),
                            function(x, dvdata) {
                              dvs1 <- names(dvdata[[x[1]]])
                              dvs2 <- names(dvdata[[x[2]]])
                              apply(expand.grid(dvs1, dvs2), 1, function(y) {
                                paste(y, collapse=":")})
                            }, dvdata=dvdata)
  names(interactiondvs) <- remainingset
  dvs <- c(maineffectdvs, interactiondvs)

  iterationexit <- FALSE
  while (iterationexit == FALSE) {

    roundnumber <- roundnumber + 1
    remainingsetdvs <- dvs[remainingset]
    formulas <- lapply(remainingsetdvs, function(x) {
      stats::update.formula(refformula,
                            paste("~ . +", paste(x, collapse = " + ")))})
    ctable <- .compare(formulas, refformula, df, test, algorithm)
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
      if (retest) {
        remainingset <- testedEVs
      } else {
        remainingset <- testedEVs[ctable$P[seq(nrow(ctable))[-1]] < alpha]
      }
      remainingset <- remainingset[!is.na(remainingset)]
    }

    if (quiet == FALSE) {message(paste("Round", roundnumber, "complete."))}

    if (nrow(ctable) == 1 || ctable$P[1] > alpha) { iterationexit <- TRUE }
    if (!iterationexit && !retest) {
      if (all(ctable$P[seq(nrow(ctable))[-1]] >= alpha |
              is.na(ctable$P[seq(nrow(ctable))[-1]]))) { iterationexit <- TRUE }
    }
  }

  if (algorithm == "maxent") {
    selectedmodel <- .runIWLR(refformula, df)
  } else {
    selectedmodel <- .runLR(refformula, df)
  }
  selectedsetni <- selectedset[selectedset %in% names(dvdata)]
  return(list(dvdata[selectedsetni], modeltable, selectedmodel))

}