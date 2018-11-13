#' Derive variables by transformation.
#'
#' \code{deriveVars} produces derived variables from explanatory variables by
#' transformation, and returns a list of dataframes. The available
#' transformation types are as follows, described in Halvorsen et al. (2015): L,
#' M, D, HF, HR, T (for continuous EVs), and B (for categorical EVs). For spline
#' transformation types (HF, HR, T),  a subset of possible DVs is pre-selected
#' by the criteria described under Details.
#'
#' The linear transformation "L" is a simple rescaling to the range [0, 1].
#'
#' The monotonous transformation "M" performed is a zero-skew transformation
#' (Oekland et al. 2001).
#'
#' The deviation transformation "D" is performed around an optimum EV value that
#' is found by looking at frequency of presence (see \code{\link{plotFOP}}).
#' Three deviation transformations are created with different steepness and
#' curvature around the optimum.
#'
#' For spline transformations ("HF", "HR", and "T"), DVs are created around 20
#' different break points (knots) which span the range of the EV. Only DVs which
#' satisfy all of the following criteria are retained: \enumerate{ \item 3 <=
#' knot <= 18 (DVs with knots at the extremes of the EV are never retained).
#' \item Chi-square test of the single-variable model from the given DV compared
#' to the null model gives a p-value < 0.05. \item The single-variable model
#' from the given DV shows a local maximum in fraction of variation explained
#' (D^2, sensu Guisan & Zimmerman, 2000) compared to DVs from the neighboring 4
#' knots.} The models used in this pre-selection procedure may be maxent models
#' (algorithm="maxent") or standard logistic regression models (algorithm="LR").
#'
#' For categorical variables, 1 binary derived variable (type "B") is created
#' for each category.
#'
#' The maximum entropy algorithm ("maxent") --- which is implemented in
#' MIAmaxent as an infinitely-weighted logistic regression with presences added
#' to the background --- is conventionally used with presence-only occurrence
#' data. In contrast, standard logistic regression (algorithm = "LR"), is
#' conventionally used with presence-absence occurrence data.
#'
#' Explanatory variables should be uniquely named. Underscores ('_') and colons
#' (':') are reserved to denote derived variables and interaction terms
#' respectively, and \code{deriveVars} will replace these --- along with other
#' special characters --- with periods ('.').
#'
#' @param data Data frame containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent either presence and background (coded as 1/NA) or presence
#'   and absence (coded as 1/0). The explanatory variable data should be
#'   complete (no NAs). See \code{\link{readData}}.
#' @param transformtype Specifies the types of transformations types to be
#'   performed. Default is the full set of the following transformation types: L
#'   (linear), M (monotonous), D (deviation), HF (forward hinge), HR (reverse
#'   hinge), T (threshold), and B (binary).
#' @param allsplines Logical. Keep all spline transformations created, rather
#'   than pre-selecting particular splines based on fraction of total variation
#'   explained.
#' @param algorithm Character string matching either "maxent" or "LR", which
#'   determines the type of model used for spline pre-selection. See Details.
#' @param write Logical. Write the transformation functions to .Rdata file?
#'   Default is \code{FALSE}.
#' @param dir Directory for file writing if \code{write = TRUE}. Defaults to the
#'   working directory.
#' @param quiet Logical. Suppress progress messages from spline pre-selection?
#'
#' @return List of 2: \enumerate{ \item dvdata: List containing first the
#'   response variable, followed data frames of derived variables produced for
#'   each explanatory variable. This item is recommended as input for
#'   \code{dvdata} in \code{\link{selectDVforEV}}. \item transformations: List
#'   containing first the response variable, followed by all the transformation
#'   functions used to produce the derived variables. }
#'
#' @references Guisan, A., & Zimmermann, N. E. (2000). Predictive habitat
#'   distribution models in ecology. Ecological modelling, 135(2-3), 147-186.
#' @references Halvorsen, R., Mazzoni, S., Bryn, A., & Bakkestuen, V. (2015).
#'   Opportunities for improved distribution modelling practice via a strict
#'   maximum likelihood interpretation of MaxEnt. Ecography, 38(2), 172-183.
#' @references Oekland, R.H., Oekland, T. & Rydgren, K. (2001).
#'   Vegetation-environment relationships of boreal spruce swamp forests in
#'   Oestmarka Nature Reserve, SE Norway. Sommerfeltia, 29, 1-190.
#'
#' @examples
#' toydata_dvs <- deriveVars(toydata_sp1po, c("L", "M", "D", "HF", "HR", "T", "B"))
#' str(toydata_dvs$dvdata)
#' summary(toydata_dvs$transformations)
#'
#' \dontrun{
#' # From vignette:
#' grasslandDVs <- deriveVars(grasslandPO,
#'                            transformtype = c("L","M","D","HF","HR","T","B"))
#' summary(grasslandDVs$dvdata)
#' head(summary(grasslandDVs$transformations))
#' length(grasslandDVs$transformations)
#' plot(grasslandPO$terslpdg, grasslandDVs$dvdata$terslpdg$terslpdg_D2, pch=20,
#'      ylab="terslpdg_D2")
#' plot(grasslandPO$terslpdg, grasslandDVs$dvdata$terslpdg$terslpdg_M, pch=20,
#'      ylab="terslpdg_M")
#' }
#'
#' @export


deriveVars <- function(data,
                       transformtype = c("L", "M", "D", "HF", "HR", "T", "B"),
                       allsplines = FALSE, algorithm = "maxent", write = FALSE,
                       dir = NULL, quiet = FALSE) {

  colnames(data) <- make.names(colnames(data), allow_ = FALSE)

    if (any(c("HF", "HR", "T") %in% transformtype) && allsplines == F) {
    .binaryrvcheck(data[, 1])
  }

  if (any(!stats::complete.cases(data[,-1]))) {
    warning(paste(sum(!stats::complete.cases(data[,-1])),
                  "rows in 'data' were dropped due to missing EV values."), call. = FALSE)
    data <- data[stats::complete.cases(data[,-1]), ]
  }

  if (write == TRUE) {
    if (is.null(dir)) { dir <- getwd() }
    fdir <- file.path(dir, "deriveVars")
    if (file.exists(fdir)) {
      yn <- readline("The specified dir already contains a deriveVars result. Overwrite this result? (y/n) ")
      if (yn == "y") {
        unlink(fdir, recursive = TRUE)
        Sys.sleep(1)
      }
      if (yn != "y") { return(message("Overwrite declined")) }
    }
    dir.create(fdir, recursive = TRUE)
  }

  transformations <- list("RV"=data[, 1])
  dvdata <- list("RV"=data[, 1])

  for (i in 2:ncol(data)) {
    df <- data[, c(1,i)]
    result <- .dvsfromev(df, transformtype, allsplines, algorithm, quiet)
    transformations <- c(transformations, result$storage)
    dvdata[[names(data)[i]]] <- result$evdv
  }

  dvdata <- c(dvdata[1], dvdata[-1][sapply(dvdata[-1], function(x) {dim(x)[2] != 0})])

  if (write == TRUE) {
    save(transformations, file = file.path(fdir, "transformations.Rdata"))
  }
  return(list("dvdata" = dvdata, "transformations" = transformations))
}
