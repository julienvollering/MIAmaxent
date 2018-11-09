#' selects a subset of spline dvs based on Maxent FTVE
#'
#' @param rv Vector of response variable values
#' @param dvs List of spline dvs to be selected from (HF, HR, or Th)
#' @param algorithm Character string matching either "maxent" or "LR".
#' @keywords internal
#' @noRd

.splselect <- function(rv, dvs, algorithm) {

  data <- data.frame("RV"=rv, do.call(cbind, dvs))
  names(data)[-1] <- names(dvs)
  formulas <- lapply(names(dvs), function(x) {
    stats::formula(paste("RV ~", x))})
  ctable <- .compare(formulas, stats::formula("RV ~ 1"), data, "Chisq", algorithm)

  selected <- character()
  for (i in 3:(nrow(ctable)-2)) {
    if (ctable$Chisq[i] >= ctable$Chisq[i-2] &&
        ctable$Chisq[i] >= ctable$Chisq[i-1] &&
        ctable$Chisq[i] >= ctable$Chisq[i+1] &&
        ctable$Chisq[i] >= ctable$Chisq[i+2] &&
        ctable$P[i] < 0.05) {
      selected <- append(selected, ctable$variables[i])
    }
  }

  return(selected)
}