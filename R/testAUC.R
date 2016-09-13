#' Calculate model AUC with test data.
#'
#' For a given model, \code{testAUC} calculates the Area Under the Curve (AUC)
#' of the Receiver Operating Characteristic (ROC) as a threshold-independent
#' measure of binary classification performance. This function is intended to be
#' used with occurence data that is independent from the data used to train the
#' model, to obtain an unbiased measure of model performance.
#'
#' @param data Data frame containing test occurrence data in the first column
#'   and corresponding explanatory variables in subsequent columns. The test
#'   data should be coded as: 1/0/NA, representing presence, absence, and
#'   unknown. See \code{\link{readData}}.
#' @param transformation Full pathway of the 'transformation.Rdata' file
#'   containing the transformations used to build the model. This file is saved
#'   as a result of the \code{\link{deriveVars}} function. Equivalently, the
#'   second item in the list returned by \code{\link{deriveVars}} can be used
#'   directly.
#' @param model Full pathway of the '.lambdas' file of the model in question.
#'   This file is saved as a result of \code{\link{selectEV}}.
#'
#' @return In addition to returning the testAUC value, graphical output showing
#'   the corresponding ROC plot is produced. The point along the ROC curve where
#'   the discrimination threshold is PRO = 1 is shown for reference.
#'
#' @examples
#' \dontrun{
#' AUC <- testAUC(testdat,
#'    transformation = "D:/path/to/modeling/directory/deriveVars/transformations.Rdata",
#'    model = "D:/path/to/modeling/directory/selectEV/round/model/1.lambdas")
#' }
#'
#' sp1pa <- toydata_sp1po
#' sp1pa$RV[is.na(sp1pa$RV)] <- 0
#' sp1pa[, 1] <- sample(sp1pa$RV, 40)
#' auc <- testAUC(sp1pa, toydata_dvs$transformations,
#'    system.file("extdata/sommerfeltia", "1.lambdas", package = "MIAmaxent"))
#' auc
#'
#' \dontrun{
#' From vignette:
#' grasslandAUC <- testAUC(grasslandPA, transformation = grasslandDVs[[2]],
#'    model = system.file("extdata", "1.lambdas", package = "MIAmaxent"))
#' grasslandAUC
#' }
#'
#' @export


testAUC <- function(data, transformation, model) {

  data <- stats::na.omit(data)
  test <- data[, 1]
  PRO <- projectModel(data, transformation, model)[[1]][, 1]

  cont <- as.matrix(table(PRO, test))
  cont <- cont[order(as.numeric(rownames(cont)), decreasing = T), ]
  falspos <- c(0, unname(cumsum(cont[, "0"])))
  truepos <- c(0, unname(cumsum(cont[, "1"])))
  fpr <- falspos/sum(cont[, "0"])
  tpr <- truepos/sum(cont[, "1"])

  graphics::plot(fpr, tpr, type="l", col="red", cex = 0.5, xlim=c(0,1), ylim=c(0,1),
    xlab="1 - specificity (false positive rate)",
    ylab="Sensitivity (true positive rate)")
  graphics::abline(0,1, lty=3)

  PRO1fp <- sum(cont[as.numeric(rownames(cont)) > 1, "0"])
  PRO1tp <- sum(cont[as.numeric(rownames(cont)) > 1, "1"])
  x <- PRO1fp/sum(cont[, "0"])
  y <- PRO1tp/sum(cont[, "1"])

  graphics::points(x, y, pch = 19, col = "#999999")
  graphics::text(x, y, "PRO = 1", pos = 3, col = "#999999", cex = 0.9)

  hgtl <- tpr[-length(tpr)]
  hgtr <- tpr[-1]
  wdth <- diff(fpr)
  AUC <- sum(((hgtl + hgtr)/2) * wdth)

  return(AUC)
}
