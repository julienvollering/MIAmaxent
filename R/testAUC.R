#' Calculate model AUC with test data.
#'
#' For a given model, \code{testAUC} calculates the Area Under the Curve (AUC)
#' of the Receiver Operating Characteristic (ROC) as a threshold-independent
#' measure of binary classification performance. This function is intended to be
#' used with occurrence data that is independent from the data used to train the
#' model, to obtain an unbiased measure of model performance.
#'
#' If plotted, the point along the ROC curve where the discrimination threshold
#' is PRO = 1, is shown for reference.
#'
#' @param model The model to be projected, represented by an object of class
#'   'glm'. This may be the object returned by \code{\link{chooseModel}}, or the
#'   'selectedmodel' returned by \code{\link{selectEV}}.
#' @param transformations Transformation functions used to create the derived
#'   variables in the model. I.e. the 'transformations' returned by
#'   \code{\link{deriveVars}}. Equivalently, the full file pathway of the
#'   'transformations.Rdata' file saved as a result of \code{\link{deriveVars}}.
#' @param data Data frame containing test occurrence data in the first column
#'   and corresponding explanatory variables in the model in subsequent columns.
#'   The test data should be coded as: 1/0/NA, representing presence, absence,
#'   and unknown. See \code{\link{readData}}.
#' @param plot Logical. Plot the ROC curve?
#' @param ... Arguments to be passed to \code{plot} to control the appearance of
#'   the ROC plot. For example: \itemize{ \item \code{lwd} for line width \item
#'   \code{main} for plot title \item \code{cex} for plot text and symbol size }
#'   Note that some graphical parameters may return errors or warnings if they
#'   cannot be changed or correspond to multiple elements in the plot.
#'
#' @examples
#' \dontrun{
#' # From vignette:
#' grasslandPA <- readData(
#'   occurrence = system.file("extdata", "occurrence_PA.csv", package="MIAmaxent"),
#'   contEV = system.file("extdata", "EV_continuous", package="MIAmaxent"),
#'   catEV = system.file("extdata", "EV_categorical", package="MIAmaxent"),
#'   PA = TRUE, XY = TRUE)
#' head(grasslandPA)
#' tail(grasslandPA)
#' testAUC(model = grasslandmodel, transformations = grasslandDVs$transformations,
#'         data = grasslandPA)
#' }
#'
#' @export


testAUC <- function(model, transformations, data, plot = TRUE, ...) {

  dvnamesni <- names(model$betas)[grep(":", names(model$betas), invert = TRUE)]
  evnames <- unique(sub("_.*", "", dvnamesni))
  for (i in evnames) {
    if (sum(colnames(data) == i) != 1) {
      stop(paste(i, "must be represented in 'data' (exactly once)"),
           call. = FALSE) }
  }

  if (!all(levels(as.factor(data[, 1])) %in% c("1", "0"))) {
    stop("The test data is not coded as (1/0/NA)", call. = FALSE)
  }
  if (all(levels(as.factor(data[, 1])) == "1") && anyNA(data[, 1])) {
    warning("The test data consist of 1/NA only, so NA is treated as absence.
Be aware of implications for the interpretation of the AUC value.", call. = FALSE)
    data[is.na(data[, 1]), 1] <- 0
  }

  data <- stats::na.omit(data)
  test <- data[, 1]
  PRO <- projectModel(model, transformations, data)[[1]][, 1]

  cont <- as.matrix(table(PRO, test))
  cont <- cont[order(as.numeric(rownames(cont)), decreasing = T), ]
  falspos <- c(0, unname(cumsum(cont[, "0"])))
  truepos <- c(0, unname(cumsum(cont[, "1"])))
  fpr <- falspos/sum(cont[, "0"])
  tpr <- truepos/sum(cont[, "1"])

  hgtl <- tpr[-length(tpr)]
  hgtr <- tpr[-1]
  wdth <- diff(fpr)
  AUC <- sum(((hgtl + hgtr)/2) * wdth)

  if (class(model)[1] == "iwlr") {
    PROpt <- TRUE
    PRO1fp <- sum(cont[as.numeric(rownames(cont)) > 1, "0"])
    PRO1tp <- sum(cont[as.numeric(rownames(cont)) > 1, "1"])
    x <- PRO1fp/sum(cont[, "0"])
    y <- PRO1tp/sum(cont[, "1"])
  } else {
    PROpt <- FALSE
    x <- NULL
    y <- NULL
  }

  if (plot == TRUE) {
    .plotROC(fpr, tpr, AUC, PROpt, x, y, ...)
  }

  return(AUC)
}
