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
#' @param predictions Pathway to the "backgroundPredictions.csv" file of the
#'   model in question.
#' @param test Vector of test data (1/0/NA), with identical length and order as
#'   training data.
#'
#' @export


testAUC <- function(predictions = NULL, test) {

  logistic <- read.csv(predictions)[,5]
  df <- data.frame(1:length(logistic), test, logistic)
  colnames(df) <- c("ID", "test", "logistic")
  testdf <- na.omit(df)
  orderedpreds <- df$logistic[order(df$logistic)]
  breaks <- orderedpreds[-length(orderedpreds)] + diff(orderedpreds)/2

  AUC <- PresenceAbsence::auc(testdf, st.dev = FALSE, na.rm=T)
  suppressWarnings(PresenceAbsence::auc.roc.plot(df, threshold = breaks, na.rm = TRUE,
    add.legend = FALSE, line.type = FALSE, color = "red"))

  return(AUC)
}
