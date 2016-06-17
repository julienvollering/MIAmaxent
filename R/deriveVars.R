#' Derive variables by transformation.
#'
#' \code{deriveVars} produces derived variables from explanatory variables by
#' transformation, and returns a list of dataframes. The available
#' transformation types are as follows, described in Halvorsen et al. (2015): L,
#' M, D, HF, HR, T (for continuous EVs), and B (for categorical EVs). For spline
#' transformation types (HF, HR, T),  a subset of possible DVs is selected by
#' the criteria described under Details.
#'
#' The monotonous transformation "M" performed in this function is a zero-skew
#' transformation (Oekland et al. 2001).
#'
#' For spline transformations, DVs are created around 20 different break points
#' (knots). Only DVs which satisfy all of the following criteria are retained:
#' \enumerate{ \item 3 <= knot <= 18 (DVs with knots at the extremes of the EV
#' are never retained). \item F-test of the single-variable Maxent model from
#' the given DV gives a p-value < 0.05. \item The single-variable Maxent model
#' from the given DV shows a local maximum in fraction of variation explained
#' (FVA) compared to DVs from the neighboring 4 knots.}
#'
#' @param data Dataframe containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent presence/background data, coded as: 1/NA.
#' @param transformtype Specifies the types of transformations types to be
#'   performed. Default is the full set of the following transfomation types: L
#'   (linear), M (monotonous), D (deviation), HF (forward hinge), HR (reverse
#'   hinge), T (threshold), and B (binary).
#' @param allsplines Logical. Keep all spline transformations created, rather
#'   than selecting particular splines based on fraction of total variation
#'   explained.
#' @param writedir The directory to which Maxent files will be written during
#'   subset selection of spline variables. Defaults to the working directory.
#' @param jarpath The pathway to the maxent.jar executable jar file. If
#'   unspecified, the function looks for the file in the writedir.
#'
#' @return List of length 2. \enumerate{ \item A list of data frames, with each
#'   containing the derived variables produced for a given explanatory variable.
#'   This item may be used as input for \code{\link{selectDV}}. \item A list of
#'   all the transformation functions used to produce the derived variables.
#'   This item may be used as input for \code{\link{plotResp}},
#'   \code{\link{testAUC}}, and \code{\link{projectModel}} }
#'
#' @references Halvorsen, R., Mazzoni, S., Bryn, A., & Bakkestuen, V. (2015).
#'   Opportunities for improved distribution modelling practice via a strict
#'   maximum likelihood interpretation of MaxEnt. Ecography, 38(2), 172-183.
#' @references Oekland, R.H., Oekland, T. & Rydgren, K. (2001)
#'   Vegetation-environment relationships of boreal spruce swamp forests in
#'   Oestmarka Nature Reserve, SE Norway. Natural History Museums and Botanical
#'   Garden, University of Oslo, Oslo.

#'
#' @export


deriveVars <- function(data,
                       transformtype = c("L", "M", "D", "HF", "HR", "T", "B"),
                       allsplines = FALSE, writedir = NULL, jarpath = NULL) {

  if (any(c("HF", "HR", "T") %in% transformtype) && allsplines == F) {
    altrMaxent:::.binaryrvcheck(data[,1])
  }

  if (is.null(writedir)) {
    writedir <- getwd()
  }

  if (is.null(jarpath)) {
    jarpath <- paste(writedir, "\\maxent.jar", sep="")
  }

  if (file.exists(jarpath) == F) {
    stop("maxent.jar file must be present in writedir, or its pathway must be
specified by the jarpath parameter. \n ")
  }

  dir <- paste(writedir, "\\deriveVars", sep="")
  if (file.exists(dir)) {
    stop("The specified writedir already contains a selection of spline DVs.
Please specify a different writedir. \n ", call. = FALSE)
  } else {
    dir.create(dir)
  }


  Storage <- list()
  EVDV <- list()
  for (i in 2:ncol(data)) {
    df <- data[,c(1,i)]
    result <- .dvfromevstored(df, transformtype, allsplines, dir, jarpath)
    Storage <- c(Storage, result$storage)
    EVDV[[colnames(df)[2]]] <- result$evdv
  }

  save(Storage, file = paste0(dir, "\\transf_storage.Rdata"))
  return(list("EVDV" = EVDV, "transformations" = Storage))
}
