if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("EV", "int", "PRO", "RV", "n"))
}



#' checks the validity of RV values
#'
#' Presence-only data should be coded as: 1/NA (preferred) or 1/0 (danger of
#' misinterpretation as presence/absence data)
#'
#' @param rv Vector of response variable values

.binaryrvcheck <- function(rv) {
  if (length(levels(as.factor(rv))) > 2) {
    stop("The response variable must contain 2 levels only: presence (1)
       and background (NA/0)", call. = FALSE)
  }
  if (anyNA(rv) && length(levels(as.factor(rv))) > 1) {
    stop("The response variable must contain 2 levels only: presence (1)
       and background (NA/0)", call. = FALSE)
  }
  if (class(rv) != "numeric" && class(rv) != "integer") {
    stop("The response variable must be numeric or integer class: presence (1)
       and background (NA/0)", call. = FALSE)
  }
}



#' checks representation of dvs in tranformations
#'
#' @param dvnamesni Names of DVs in model (no interaction terms)
#' @param alltransf List of transformation functions

.check.dvs.in.transf <- function(dvnamesni, alltransf) {
  for (i in dvnamesni) {
    a <- paste0(i, "_transf")
    if (sum(names(alltransf) == a) != 1) {
      stop(paste(i, "must be represented in 'transformations' (exactly once)"),
        call. = FALSE)
    }
  }
}



#' Name and create directory
#'
#' Simultaneuosly pastes arguments into pathway and creates the directory
#'
#' @param ... Arguments to be pasted together into directory pathway

.dirpath.create <- function(...) {
  path <- file.path(...)
  dir.create(path)
  return(path)
}



#' calculates optimum EV value based on FOP
#'
#' The optimum that is returned is based on the smoothed data, unless a maximum
#' exists at the extremes of the EV (outside the 5-interval smoothing window).
#'
#' @param data Dataframe containing the response variable in the first column and
#'   explanatory variables in the second column. The response variable should
#'   represent presence or background, coded as: 1/NA.
#' @param smoothwindow Width of the smoothing window (in an exponentially
#'   weighted moving average). Irrelevant for categorical EVs.
#' @param intervals Number of intervals into which the continuous EV is divided.
#'   Defaults to the minimum of N/50 and 100. Irrelevant for categorical EVs.
#'
#' @return the EV value at which FOP is highest (\code{EVoptimum})

.fopoptimum <- function(data, smoothwindow = 5, intervals = NULL) {

  df <- data.frame(RV = data[, 1], EV = data[, 2])
  .binaryrvcheck(df[, 1])
  df[, 1][is.na(df[, 1])] <- 0

  if (!class(df[, 2]) %in% c("numeric", "integer")) {
    stop("EVoptimum is calculated for numeric or integer class EVs only",
      call. = F)
  }

  if (is.null(intervals)) {intervals <- min(c(ceiling(nrow(df) / 50), 100))}
  df$int <- cut(df[, 2], max(2, intervals))

  grouped <- dplyr::group_by(df, int)
  FOPdf <- dplyr::summarise(grouped, intEV = mean(EV), intRV = mean(RV, na.rm=F))

  if (length(FOPdf$intRV) > smoothwindow) {
    FOPdf$smoothRV <- .ewma(FOPdf$intRV, smoothwindow)
  } else { FOPdf$smoothRV <- NA }

  maxRV <- FOPdf$smoothRV
  maxRV[is.na(maxRV)] <- FOPdf$intRV[is.na(maxRV)]
  EVoptimum <- FOPdf$intEV[which(maxRV == max(maxRV))]

  while (length(EVoptimum) > 1) {
    intervals <- intervals - 1
    df$int <- cut(df[, 2], max(2, intervals))
    grouped <- dplyr::group_by(df, int)
    FOPdf <- as.data.frame(dplyr::summarise(grouped, n = n(),
      intEV = mean(EV), intRV = mean(RV, na.rm=F)))
    if (length(FOPdf$intRV) > smoothwindow) {
      FOPdf$smoothRV <- .ewma(FOPdf$intRV, smoothwindow)
    } else { FOPdf$smoothRV <- NA }
    maxRV <- FOPdf$smoothRV
    maxRV[is.na(maxRV)] <- FOPdf$intRV[is.na(maxRV)]
    EVoptimum <- FOPdf$intEV[which(maxRV == max(maxRV))]
  }

  return(EVoptimum)
}

#' checks the validity of formulas
#'
#' @param formula Formula entered as selection start point
#' @param dvdata List of data frames containing EVs

.formulacheck <- function(formula, dvdata) {
  if (any(attr(stats::terms(formula), "order") != 1)) {
    stop("The provided formula may contain first-order explanatory variables
      only (no interactions)", call. = FALSE)
  }
  terms <- labels(stats::terms(formula))
  for (i in terms) {
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

.plotROC <- function(fpr, tpr, AUC, x, y, ...) {
  args1 <- list(xlab="1 - specificity (false positive rate)",
                ylab="Sensitivity (true positive rate)", col="red",
                main=paste("AUC = ", signif(AUC, 3)))
  inargs <- list(...)
  args1[names(inargs)] <- inargs
  do.call(graphics::plot, c(list(x=fpr, y=tpr, xlim=c(0,1), ylim=c(0,1),
                                 type="l"), args1))

  graphics::abline(0, 1, lty=3)

  args2 <- list(cex=0.8, col="#999999", pch=19)
  inargs <- list(...)
  args2[names(inargs)] <- inargs
  do.call(graphics::points, c(list(x=x, y=y), args2))

  args3 <- list(cex=0.8, col="#999999")
  inargs <- list(...)
  args3[names(inargs)] <- inargs
  do.call(graphics::text, c(list(x=x, y=y, labels="PRO = 1", pos=4), args3))
}

#' Reminders when using devtools::release
#'
#' @keywords internal

release_questions <- function() {
  c(
    "Have you reknitted the static vignette and copied the html file into /vignettes?",
    "Have you removed the vignitte-produced directories?",
    "Have you removed 'maxent.jar' from inst/java?"
  )
}



#' skewness transformation using constant c
#'
#' @param x Vector of data.
#' @param c Constant

.scalex <- function(xnull, x, c) {
  if(e1071::skewness(xnull, na.rm = TRUE, type = 2) < 0) {
    return(exp(c * x))
  } else {
    return(log(x + c))
  }
}
