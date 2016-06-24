#' Calculate model AUC with test data
#'
#' For a given model, \code{testAUC} calculates the Area Under the Curve (AUC)
#' of the Receiver Operating Characteristic (ROC) as a threshold-independent
#' measure of binary classification performance. This function is intended to be
#' used with occurence data that is independent from the data used to train the
#' model, to obtain an unbiased measure of model performance.
#'
#' \code{DESCRIPTION Imports}: PresenceAbsence
#'
#' @param data Data frame containing test data in the first column and
#'   corresponding explanatory variables in subsequent columns. The test data
#'   should be coded as: 1/0/NA, representing presence, absence, and unknown.
#' @param transformation Pathway to the 'transformation.Rdata' file containing
#'   the transformations used to build the model. This file is saved as a result
#'   of the \code{\link{deriveVars}} function. Equivalently, the second item in
#'   the list returned by \code{\link{deriveVars}} can be used directly.
#' @param model Pathway to the '.lambdas' file of the model in question. This
#'   file is saved as a result of \code{\link{selectEV}}.
#'
#' @return In addition to returning the testAUC value, graphical output showing
#'   the corresponding ROC plot is produced.
#'
#' @export


testAUC <- function(data, transformation, model) {

  data <- na.omit(data)
  test <- data[, 1]
  PRO <- projectModel(data, transformation, model)[[1]][, 1]
  rangedoutput <- (PRO - min(PRO)) / diff(range(PRO))
  PRO1 <- (1 - min(PRO)) / diff(range(PRO))
  df <- data.frame(ID = 1:nrow(data), test = test, output = rangedoutput)

  AUC <- PresenceAbsence::auc(df, st.dev = FALSE)

  orderedpreds <- df$output[order(df$output)]
  breaks <- orderedpreds[-length(orderedpreds)] + diff(orderedpreds)/2
  PresenceAbsence::auc.roc.plot(df, threshold = breaks, add.legend = FALSE,
    line.type = FALSE, color = "red")

  PRO1thresh <- PresenceAbsence::roc.plot.calculate(df, PRO1)
  x <- 1 - PRO1thresh$specificity
  y <- PRO1thresh$sensitivity
  points(x, y, pch = 19, col = "#999999")
  text(x, y, "PRO = 1", pos = 3, col = "#999999", cex = 0.9)

  return(AUC)
}
