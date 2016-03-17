#' Derive variables by transformation.
#'
#' \code{deriveVars} produces derived variables from explanatory variables by
#' the following transformation types as described in Halvorsen et al. (2015):
#' L, M, D, HF, HR, T (for continuous EVs), and B (for categorical EVs). For
#' spline transformation types (HF, HR, T),  a subset of possible DVs is
#' selected by the criteria described in Halvorsen et al. (2015) p. 179.
#'
#' (Details)
#'
#' @param data Dataframe containing the response variable in the first column
#'   and explanatory variables in subsequent columns. The response variable
#'   should represent presence/background data, coded as: 1/NA.
#' @param writedir The directory to which Maxent files will be written during
#'   subset selection of spline variables. Defaults to the working directory.
#' @param transformtype Specifies the types of transformations types to be
#'   performed. Default is the full set of the following transfomation types: L
#'   (linear), M (monotonous), D (deviation), HF (forward hinge), HR (reverse
#'   hinge), Th (threshold), and B (binary).
#' @param allsplines Logical. Keep all spline transformations created, rather
#'   than selecting particular splines based on faction of total variation
#'   explained.
#'
#' @return List of dataframes, with each data frame containing the derived
#'   variables produced for a given explanatory variable by transformation.
#'   Dimensions: [[EV]] [N, DV].
#'
#' @references Halvorsen, R., Mazzoni, S., Bryn, A., & Bakkestuen, V. (2015).
#'   Opportunities for improved distribution modelling practice via a strict
#'   maximum likelihood interpretation of MaxEnt. Ecography, 38(2), 172-183.
#'
#' @export


deriveVars <- function(data, writedir = NULL,
  transformtype = c("L", "M", "D", "HF", "HR", "T", "B"), allsplines = FALSE) {

  if (any(c("HF", "HR", "T") %in% transformtype) && allsplines == F) {
    altrMaxent:::.binaryrvcheck(data[,1])
  }
  if (is.null(writedir)) {
    writedir <- getwd()
  }

  EVDV <- list()
  for (i in 2:ncol(data)) {
    df <- data[,c(1,i)]
    EVDV[[i]] <- altrMaxent:::.dvfromev(df, writedir, transformtype, allsplines)
    names(EVDV[[i]]) <- colnames(data)[i]
  }

  return(EVDV)
}
