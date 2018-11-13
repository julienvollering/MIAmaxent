#' MIAmaxent: A Modular, Integrated Approach to Maximum Entropy Distribution
#' Modeling
#'
#' Tools for training, selecting, and evaluating maximum entropy (and standard
#' logistic regression) distribution models. This package provides tools for
#' user-controlled transformation of explanatory variables,  selection of
#' variables by nested model comparison, and flexible model evaluation and
#' projection. It follows principles based on the maximum-likelihood
#' interpretation of maximum entropy modeling (Halvorsen et al., 2015), and uses
#' infinitely-weighted logistic regression for model fitting. (Fithian & Hastie,
#' 2013).
#'
#' MIAmaxent is intended primarily for maximum entropy distribution modeling
#' (Phillips et al., 2006; Phillips et al., 2017), but provides an alternative
#' to the standard methodology for training, selecting, and using models. The
#' major advantage in this alternative methodology is greater user control -- in
#' variable transformations, in variable selection, and in model output.
#' Comparisons also suggest that this methodology results in simpler models with
#' equally good predictive ability, and reduces the risk of overfitting
#' (Halvorsen et al., 2016).
#'
#' The predecessor to this package is the MIA Toolbox, which is described in
#' detail in Mazzoni et al. (2015).
#'
#' @section User Workflow: This diagram outlines a common workflow for users of
#'   this package. Functions are shown in red.
#'
#'   \if{html}{\figure{workflow-flowchart.png}{options: width="70\%"
#'   alt="Figure: workflow-flowchart.png"}}
#'   \if{latex}{\figure{workflow-flowchart.pdf}{options: width=12cm}}
#'
#' @references Fithian, W., & Hastie, T. (2013). Finite-sample equivalence in
#'   statistical models for presence-only data. The annals of applied
#'   statistics, 7(4), 1917.
#' @references Halvorsen, R., Mazzoni, S., Bryn, A. & Bakkestuen, V. (2015)
#'   Opportunities for improved distribution modelling practice via a strict
#'   maximum likelihood interpretation of MaxEnt. Ecography, 38, 172-183.
#' @references Halvorsen, R., Mazzoni, S., Dirksen, J.W., Naesset, E., Gobakken,
#'   T. & Ohlson, M. (2016) How important are choice of model selection method
#'   and spatial autocorrelation of presence data for distribution modelling by
#'   MaxEnt? Ecological Modelling, 328, 108-118.
#' @references Mazzoni, S., Halvorsen, R. & Bakkestuen, V. (2015) MIAT: Modular
#'   R-wrappers for flexible implementation of MaxEnt distribution modelling.
#'   Ecological Informatics, 30, 215-221.
#' @references Phillips, S.J., Anderson, R.P., Dudik, M., Schapire, R.E., &
#'   Blair, M.E. (2017). Opening the black box: an open-source release of
#'   Maxent. Ecography, 40(7), 887-893.
#' @references Phillips, S.J., Anderson, R.P. & Schapire, R.E. (2006) Maximum
#'   entropy modeling of species geographic distributions. Ecological Modelling,
#'   190, 231-259.
#'
#' @docType package
#' @name MIAmaxent
#' @keywords internal
#'
NULL