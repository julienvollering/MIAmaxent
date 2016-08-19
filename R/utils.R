if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("EV", "int", "PRO", "RV"))
}



#' Best prefix match.
#'
#' Return only the string in 'a' which best matches any of the strings in 'b'
#'
#' @param a Vector of character strings to be matched
#' @param b vector of character strings to match
#'
#' @return single character string

.best.match <- function(a, b) {
  score <- sapply(a, function(a) {
    reg <- regexpr(a, b)
    max(attr(reg, "match.length"))
  })
  a[which.max(score)]
}



#' Index of best prefix match.
#'
#' Return the index of the string in 'a' which best matches any of the strings
#' in 'b'
#'
#' @param a Vector of character strings to be matched
#' @param b vector of character strings to match
#'
#' @return single character string

.best.match.ind <- function(a, b) {
  score <- sapply(a, function(a) {
    reg <- regexpr(b, a)
    max(attr(reg, "match.length"))
  })
  maxs <- score[score == max(score)]
  best <- names(maxs[nchar(names(maxs)) == min(nchar(names(maxs)))])
  best.ind <- which(a == best)
  return(best.ind)
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



#' calculates exponentially weighted moving average
#'
#' @param x numeric. Vector across which the moving average is to be applied.
#' @param n integer. Width of the moving average window. Should be odd,
#'   otherwise the window will be uncentered.
#'
#' @return vector of moving average values

.ewma <- function(x, n) {
  if (missing(n)) {
    stop("Specify the width of the moving average window (n)", call. = FALSE)
  }
  if (n < 3) {
    stop("Width of window should be at least 3", call. = FALSE)
  }

  if (n %% 2 != 0) {
    expwindow <- stats::dexp(c(((n-1)/2):0,1:((n-1)/2)))
  } else {
    expwindow <- stats::dexp(c((n/2):0,1:((n-2)/2)))
  }
  weights <- expwindow/sum(expwindow)
  as.numeric(stats::filter(x, weights, sides=2))
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
  df$int <- .reg.interval(df[, 2], intervals)

  grouped <- dplyr::group_by(df, int)
  FOPdf <- dplyr::summarise(grouped, intEV = mean(EV), intRV = mean(RV, na.rm=F))

  FOPdf$smoothRV <- .ewma(FOPdf$intRV, smoothwindow)
  maxRV <- FOPdf$smoothRV
  maxRV[is.na(maxRV)] <- FOPdf$intRV[is.na(maxRV)]
  EVoptimum <- FOPdf$intEV[which(maxRV == max(maxRV))]

  while (length(EVoptimum) > 1) {
    intervals <- intervals - 1
    df$int <- .reg.interval(df[, 2], intervals)
    grouped <- dplyr::group_by(df, int)
    FOPdf <- dplyr::summarise(grouped, intEV = mean(EV), intRV = mean(RV, na.rm=F))
    FOPdf$smoothRV <- .ewma(FOPdf$intRV, smoothwindow)
    maxRV <- FOPdf$smoothRV
    maxRV[is.na(maxRV)] <- FOPdf$intRV[is.na(maxRV)]
    EVoptimum <- FOPdf$intEV[which(maxRV == max(maxRV))]
  }
  return(EVoptimum)
}



#' Loads a transformation object
#'
#' From .Rdata file or from existing object
#'
#' @param transformation Vector of data. Must have scale [0,1]!

.load.transf <- function(transformation) {
  if (class(transformation) == "character") {
    alltransf <- get(load(transformation))
  } else {
    alltransf <- transformation
  }
  if (!all(sapply(alltransf, class) == "function")) {
    stop("transformation argument should contain functions only", call. = FALSE)
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



#' Make regular intervals
#'
#' @param a Numeric vector
#' @param b number of intervals
#'
#' @return factor variable with 1 level for each interval

.reg.interval <- function(a, b) {
  intwidth <- (max(a) - min(a)) / b
  cutpts <- seq(min(a), max(a), by = intwidth)
  Hmisc::cut2(a, cuts = cutpts, oneval = FALSE)
}



#' Executes basic maxent.jar run from R
#'
#' @param rv Vector of response variable values
#' @param ev Data frame of explanatory variables
#' @param maxbkg Maximum number of uninformed background points to use for
#'   training
#' @param dir Directory to which Maxent files will be written

.runjar <- function(rv, ev, maxbkg = 10000, dir) {
  jarpath <- system.file("java/maxent.jar", package = "maxentmodelselectr")
  if (file.exists(jarpath) == FALSE) {
    stop("Missing 'maxent.jar' file. Place this file in the java folder of the
       package (see System Requirements in package description).", call. = FALSE)}
  df <- data.frame("RV" = rv, "X" = -9999, "Y" = -9999, ev, check.names = FALSE)
  samplesdf <- stats::na.omit(df)
  environlayersdf <- df
  csvfiles <- file.path(dir, c("samples.csv", "environlayers.csv"))
  utils::write.csv(samplesdf, csvfiles[1], row.names = F)
  utils::write.csv(environlayersdf, csvfiles[2], row.names = F)

  jarflags1 <- " removeduplicates=FALSE addsamplestobackground=FALSE"
  jarflags2 <- " autofeature=FALSE betamultiplier=0"
  jarflags3 <- " quadratic=FALSE product=FALSE hinge=FALSE threshold=FALSE"
  jarflags4 <- " outputformat=raw writebackgroundpredictions=TRUE"
  jarflags5 <- " outputgrids=FALSE pictures=FALSE"
  jarflags6 <- " extrapolate=FALSE writemess=FALSE plots=FALSE"
  jarflags7 <- " doclamp=FALSE writeclampgrid=FALSE"
  jarflags8 <- " autorun=TRUE threads=2 visible=FALSE warnings=FALSE"
  jarflags <- paste0(jarflags1, jarflags2, jarflags3, jarflags4, jarflags5,
    jarflags6, jarflags7, jarflags8)

  command <- paste0("java -mx512m -jar ",
    "\"", jarpath, "\"",
    jarflags, " maximumbackground=", maxbkg,
    " samplesfile=","\"", csvfiles[1], "\"",
    " environmentallayers=", "\"", csvfiles[2], "\"",
    " outputdirectory=", "\"", dir, "\\", "\"")
  javacommand <- gsub("\\\\","/", command)
  system(paste(javacommand), wait = TRUE)
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
