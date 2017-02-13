#' Derive variables by transformation.
#'
#' \code{deriveVars} produces derived variables from explanatory variables by
#' transformation, and returns a list of dataframes. The available
#' transformation types are as follows, described in Halvorsen et al. (2015): L,
#' M, D, HF, HR, T (for continuous EVs), and B (for categorical EVs). For spline
#' transformation types (HF, HR, T),  a subset of possible DVs is selected by
#' the criteria described under Details.
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
#' \item F-test of the single-variable Maxent model from the given DV gives a
#' p-value < 0.05. \item The single-variable Maxent model from the given DV
#' shows a local maximum in fraction of variation explained (FVA) compared to
#' DVs from the neighboring 4 knots.}
#'
#' For categorical variables, 1 binary derived variable (type "B") is created
#' for each category.
#'
#' Explanatory variables should be uniquely named, and the names must not
#' contain spaces, underscores, or colons. Underscores and colons are reserved
#' to denote derived variables and interaction terms repectively.
#'
#' @param data Data frame containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent presence/background data, coded as: 1/NA. See
#'   \code{\link{readData}}.
#' @param transformtype Specifies the types of transformations types to be
#'   performed. Default is the full set of the following transfomation types: L
#'   (linear), M (monotonous), D (deviation), HF (forward hinge), HR (reverse
#'   hinge), T (threshold), and B (binary).
#' @param allsplines Logical. Keep all spline transformations created, rather
#'   than selecting particular splines based on fraction of total variation
#'   explained.
#' @param dir Directory to which files will be written during selection of
#'   spline-type derived variables. Defaults to the working directory.
#'
#' @return List of 2: \enumerate{ \item A list of data frames, with each
#'   containing the derived variables produced for a given explanatory variable.
#'   This item is recommended as input for \code{dvdata} in
#'   \code{\link{selectDVforEV}}. \item A list of all the transformation
#'   functions used to produce the derived variables. This item is recommended
#'   as input for \code{transformation} in \code{\link{plotResp2}},
#'   \code{\link{testAUC}}, and \code{\link{projectModel}}. }
#'
#' @references Halvorsen, R., Mazzoni, S., Bryn, A., & Bakkestuen, V. (2015).
#'   Opportunities for improved distribution modelling practice via a strict
#'   maximum likelihood interpretation of MaxEnt. Ecography, 38(2), 172-183.
#' @references Oekland, R.H., Oekland, T. & Rydgren, K. (2001).
#'   Vegetation-environment relationships of boreal spruce swamp forests in
#'   Oestmarka Nature Reserve, SE Norway. Sommerfeltia, 29, 1-190.
#'
#' @examples
#' \dontrun{
#' deriveddat <- deriveVars(dat, transformtype = c("HF", "HR", "T"), allsplines = TRUE,
#'    dir = "D:/path/to/modeling/directory")
#' }
#'
#' toydata_dvs <- deriveVars(toydata_sp1po, transformtype = c("L", "M", "D", "B"))
#' str(toydata_dvs$EVDV)
#' summary(toydata_dvs$transformations)
#'
#' \dontrun{
#' # From vignette:
#' grasslandDVs <- deriveVars(grasslandPO, transformtype = c("L", "M", "D", "HF", "HR", "T", "B"))
#' summary(grasslandDVs$EVDV) # alternatively: summary(grasslandDVs[[1]])
#' head(summary(grasslandDVs$transformations)) # alternatively: head(summary(grasslandDVs[[2]]))
#' length(grasslandDVs$transformations)
#' plot(grasslandPO$terslpdg, grasslandDVs$EVDV$terslpdg$terslpdg_D2, pch = 20)
#' plot(grasslandPO$terslpdg, grasslandDVs$EVDV$terslpdg$terslpdg_HR4, pch = 20)
#' }
#'
#' @export


deriveVars <- function(data,
                       transformtype = c("L", "M", "D", "HF", "HR", "T", "B"),
                       allsplines = FALSE, dir = NULL) {

  if (any(c("HF", "HR", "T") %in% transformtype) && allsplines == F) {
    .binaryrvcheck(data[, 1])
  }

  if (is.null(dir)) { dir <- getwd()}

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

  transformations <- list()
  EVDV <- list()
  for (i in 2:ncol(data)) {
    df <- data[, c(1,i)]
    result <- .dvsfromev(df, transformtype, allsplines, fdir)
    transformations <- c(transformations, result$storage)
    EVDV[[colnames(df)[2]]] <- result$evdv
  }

  save(EVDV, file = file.path(fdir, "EVDV.Rdata"))
  save(transformations, file = file.path(fdir, "transformations.Rdata"))
  return(list("EVDV" = EVDV, "transformations" = transformations))
}
