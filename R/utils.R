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



#' checks representation of dvs in data
#'
#' @param dvnamesni Names of DVs in model (no interaction terms)
#' @param data Data frame with EV column names

.check.dvs.in.data <- function(dvnamesni, data) {
  for (i in dvnamesni) {
    a <- sub("_.*", "", i)
    if (sum(colnames(data) == a) != 1) {
      stop(paste(a, "must be represented in 'data' (exactly once)"),
        call. = FALSE)
    }
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
      stop(paste(i, "must be represented in 'transformation' (exactly once)"),
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



#' Executes basic maxent.jar run from R
#'
#' @param rv Vector of response variable values
#' @param ev Data frame of explanatory variables
#' @param maxbkg Maximum number of uninformed background points to use for
#'   training
#' @param dir Directory to which Maxent files will be written

.runjar <- function(rv, ev, maxbkg = 10000, dir) {
  jarpath <- system.file("java/maxent.jar", package = "MIAmaxent")
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



#' Creates infinitely weighted logistic regression model (equivalent to maxent)
#' without regularization
#'
#' @param formula Object of class "formula": a symbolic description of the model
#'   to be fitted. Do not use '.' term, as weights are concatenated to the data
#'   object.
#' @param data Data frame or list containing the variables in the model.
#'   Response variable as (1/NA).

.runIWLR <- function(formula, data) {
  RV <- all.vars(formula)[1]
  data[,RV][is.na(data[,RV])] <- 0
  padd <- data[data[,RV]==1, ]
  padd[, RV] <- 0
  padddata <- rbind(data, padd)
  # Code below this line was modified from the MIT-licensed 'maxnet' library
  wgts <- padddata[,1]+(1-padddata[,1])*100
  glmdata <- cbind(padddata, wgts)
  model <- stats::glm(formula=formula, family=binomial, data=glmdata,
                      weights=wgts)
  if (any(is.na(model$coefficients))) {
    nacoef <- names(model$coefficients)[is.na(model$coefficients)]
    model$betas <- model$coefficients[-1][!is.na(model$coefficients[-1])]
    formula <- stats::update(formula, paste("~ . -", nacoef))
  } else {
    model$betas <- model$coefficients[-1]
    }
  bkg <- model.matrix(stats::update(formula, ~. -1),
                      padddata[padddata[, RV]==0, ])
  model$alpha <- 0
  link <- (bkg %*% model$betas) + model$alpha
  rr <- exp(link)
  raw <- rr / sum(rr)
  model$entropy <- -sum(raw * log(raw), na.rm = TRUE)
  model$alpha <- -log(sum(rr))
  return(model)
  # Code above this line was modified from the MIT-licensed 'maxnet' library
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
