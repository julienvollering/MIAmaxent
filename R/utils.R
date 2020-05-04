if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("EV", "int", "PRO", "RV", "n"))
}



#' checks the validity of RV values
#'
#' Presence-only data should be coded as: 1/NA (preferred) or 1/0 (danger of
#' misinterpretation as presence/absence data)
#'
#' @param rv Vector of response variable values
#' @keywords internal
#' @noRd

.binaryrvcheck <- function(rv) {
  if (class(rv) != "numeric" && class(rv) != "integer") {
    stop("The response variable must be numeric or integer class: presence (1)
         and either background or absence (NA/0)", call. = FALSE)
  }
  if (anyNA(rv) && !all(levels(as.factor(rv)) %in% "1")) {
    stop("The response variable must contain exactly 2 levels: presence (1)
       and either background or absence (NA/0)", call. = FALSE)
  }
  if (!anyNA(rv) && !all(levels(as.factor(rv)) %in% c("1", "0"))) {
    stop("The response variable must contain exactly 2 levels: presence (1)
       and either background or absence (NA/0)", call. = FALSE)
  }
}



#' checks representation of dvs in tranformations
#'
#' @param dvnamesni Names of DVs in model (no interaction terms)
#' @param alltransf List of transformation functions
#' @keywords internal
#' @noRd

.check.dvs.in.transf <- function(dvnamesni, alltransf) {
  for (i in dvnamesni) {
    a <- paste0(i, "_transf")
    if (sum(names(alltransf) == a) != 1) {
      stop(paste(i, "must be represented in 'transformations' (exactly once)"),
        call. = FALSE)
    }
  }
}



#' calculates optimum EV value based on FOP
#'
#' The optimum that is returned is based on the loess-smoothed data (for
#' continuous variables).
#'
#' @param data Dataframe containing the response variable in the first column
#'   and explanatory variables in the second column. The response variable
#'   should represent presence or background, coded as: 1/NA.
#' @param span The proportion of FOP points included in the local regression
#'   neighborhood. Should be between 0 and 1. Irrelevant for categorical EVs.
#' @param intervals Number of intervals into which the continuous EV is divided.
#'   Defaults to the minimum of N/10 and 100. Irrelevant for categorical EVs.
#'
#' @return the EV value at which FOP is highest (\code{EVoptimum})
#' @keywords internal
#' @noRd

.fopoptimum <- function(df, span = 0.5, intervals = NULL) {

  df <- data.frame(RV = df[, 1], EV = df[, 2])
  .binaryrvcheck(df[, 1])
  df[, 1][is.na(df[, 1])] <- 0

  if (class(df[, 2]) %in% c("numeric", "integer")) {
    if (is.null(intervals)) {intervals <- min(c(ceiling(nrow(df)/10), 100))}
    df$int <- cut(df[, 2], breaks=max(2, intervals))
    grouped <- dplyr::group_by(df, int)
    FOPdf <- as.data.frame(dplyr::summarise(grouped, n = dplyr::n(),
                                            intEV = mean(EV),
                                            intRV = mean(RV, na.rm=FALSE)))
    FOPdf$loess <- stats::predict(stats::loess(intRV~intEV, FOPdf,
                                               weights=FOPdf$n, span=span))
    if (any(is.na(FOPdf$loess))) {
      EVoptimum <- FOPdf$intEV[which.max(FOPdf$intRV)]
    } else { EVoptimum <- FOPdf$intEV[which.max(FOPdf$loess)]  }

  }

  if (class(df[, 2]) %in% c("factor", "character")) {
    grouped <- dplyr::group_by(df, EV)
    FOPdf <- as.data.frame(dplyr::summarise(grouped, n = dplyr::n(),
                                            lvlRV = mean(RV, na.rm=FALSE)))
    EVoptimum <- FOPdf$EV[which.max(FOPdf$lvlRV)]
  }

  return(EVoptimum)
}



#' checks the validity of formulas
#'
#' @param formula Formula entered as selection start point
#' @param dvdata List of data frames containing EVs
#' @keywords internal
#' @noRd

.formulacheck <- function(formula, dvdata) {
  if (any(attr(stats::terms(formula), "order") != 1)) {
    stop("The provided formula may contain first-order explanatory variables
      only (no interactions)", call. = FALSE)
  }
  trms <- labels(stats::terms(formula))
  for (i in trms) {
    if (sum(names(dvdata) == i) != 1) {
      stop(paste(i, "must be represented in 'dvdata' (exactly once)"),
        call. = FALSE)
    }
  }
}


#' Loads a transformations object
#'
#' From .Rdata file or from existing object
#'
#' @param transformations transformations object produced by deriveVars
#' @keywords internal
#' @noRd

.load.transf <- function(transformations) {
  if (class(transformations) == "character") {
    alltransf <- get(load(transformations))
  } else {
    alltransf <- transformations
  }
  if (!all(sapply(alltransf[-1], class) == "function")) {
    stop("'transformations' should be a transformations object returned by 'deriveVars'",
         call. = FALSE)
  }
  return(alltransf)
}



#' calculates skewness of a vector
#'
#' Also calculates the constant 'c' needed for zero-skewness transformation in
#' \code{scalex}
#'
#' @param x Vector of data. Must have scale [0,1]!
#' @keywords internal
#' @noRd

.minskew <- function(x) {
  cmin <- min(x)-10*(max(x)-min(x))
  cmax <- max(x)+10*(max(x)-min(x))
  if(e1071::skewness(x, na.rm=TRUE, type=2) >= 0 && cmin < -min(x)) {
    cmin <- -min(x)
  }
  cmid <- (cmin + cmax) / 2
  skew <- e1071::skewness(.scalex(x, x, cmid), na.rm=TRUE)
  while (abs(skew) > 1 * 10^-05 && min(abs(c(cmax, cmin)-cmid)) > 10^-10) {
    sleft <- e1071::skewness(.scalex(x, x, (cmid + cmin) / 2), na.rm = TRUE,
      type = 2)
    sright <- e1071::skewness(.scalex(x, x, (cmid + cmax) / 2), na.rm = TRUE,
      type = 2)
    if (abs(sleft) < abs(skew) && abs(sleft) < abs(sright)) {
      cmax <- cmid
      skew <- sleft
    }
    else if (abs(sright) < abs(skew)) {
      cmin <- cmid
      skew <- sright
    }
    else {
      cmin <- (cmin + cmid) / 2
      cmax <- (cmax + cmid) / 2
    }
    cmid <- (cmin + cmax) / 2
  }
  return(list(c = cmid, skew = skew))
}



#' Plotting helper function for testAUC
#'
#' @param fpr false positive rate vector
#' @param tpr true positive rate vector
#' @param AUC AUC value
#' @param x PRO = 1 x-coordinate
#' @param y PRO = 1 y-coordinate
#' @keywords internal
#' @noRd

.plotROC <- function(fpr, tpr, AUC, PROpt, x, y, ...) {
  args1 <- list(xlab="1 - specificity (false positive rate)",
                ylab="Sensitivity (true positive rate)", col="red",
                main=paste("AUC = ", signif(AUC, 3)))
  inargs <- list(...)
  args1[names(inargs)] <- inargs
  do.call(graphics::plot, c(list(x=fpr, y=tpr, xlim=c(0,1), ylim=c(0,1),
                                 type="l"), args1))

  graphics::abline(0, 1, lty=3)

  if (PROpt == TRUE) {
    args2 <- list(cex=0.8, col="#999999", pch=19)
    inargs <- list(...)
    args2[names(inargs)] <- inargs
    do.call(graphics::points, c(list(x=x, y=y), args2))

    args3 <- list(cex=0.8, col="#999999")
    inargs <- list(...)
    args3[names(inargs)] <- inargs
    do.call(graphics::text, c(list(x=x, y=y, labels="PRO = 1", pos=4), args3))
  }
}



#' Reminders when using devtools::release
#'
#' @keywords internal

release_questions <- function() {
  c(
    "Have you reknitted the static vignette and copied the html file into /vignettes?"
  )
}



#' skewness transformation using constant c
#'
#' @param x Vector of data.
#' @param c Constant
#' @keywords internal
#' @noRd

.scalex <- function(xnull, x, c) {
  if(e1071::skewness(xnull, na.rm = TRUE, type = 2) < 0) {
    return(exp(c * x))
  } else {
    return(log(x + c))
  }
}


#' Return output from \code{\link{selectDVforEV}} as if it had been produced
#' under a stricter (lower) alpha. Results will match \code{selectDVforEV(...,
#' alpha = stricter, retest = TRUE)} if \code{list} was also produced with
#' retest = TRUE.
#'
#' @param list Output list from selectDVforEV()
#' @param alpha Stricter alpha than used to produce \code{list}
#' @keywords internal
#' @noRd
#' @importFrom rlang .data

.stricterselectDVforEV <- function(list, alpha) {
  dvdata <- list()
  selection <- list()

  for (i in seq_along(list$selection)) {
    evname <- names(list$selection)[i]
    drop <- FALSE
    ctable <- list$selection[[i]]
    bests <- ctable[!duplicated(ctable$round),]
    if (any(bests$P < alpha)) {
      selectedmod <- utils::tail(dplyr::filter(bests, .data$P < alpha), 1)
      lastround <- min(selectedmod$round + 1, max(bests$round))
      } else {
        lastround <- 1
        drop <- TRUE
      }
    selection[[i]] <- dplyr::filter(ctable, round <= lastround)
    names(selection)[i] <- evname
    if (!drop) {
      selectedset <- unlist(strsplit(selectedmod$variables, split=" + ",
                                     fixed=TRUE))
      dvdata[[i]] <- list$dvdata[[evname]][, selectedset, drop = FALSE]
      names(dvdata)[i] <- evname
    }
  }
  RV <- list(RV = list$dvdata$RV)
  dvdata <- c(RV, Filter(Negate(is.null), dvdata))
  return(list(dvdata = dvdata, selection = selection))
}
