#' Calculates variable contributions (FTVA)
#'
#' Calculates the Fraction of Total Variation Accounted for (Halvorsen et al.
#' 2015), for the selected model or a chosen model from the results of
#' \code{\link{selectEV}}.
#'
#' @param selectedEV The list returned by \code{selectEV}.
#' @param formula  If null, FTVA is calculated for the selected model in
#'   \code{selectedEV}. Otherwise, a model formula (in the form ~ x + ...)
#'   specifying a model in the trail of forward selection
#'   (\code{selectEV$selection}) for which to calculate FTVA. Response variable
#'   is irrelevant.
#'
#' @references Halvorsen, R., Mazzoni, S., Bryn, A., & Bakkestuen, V. (2015).
#'   Opportunities for improved distribution modelling practice via a strict
#'   maximum likelihood interpretation of MaxEnt. Ecography, 38(2), 172-183.
#'
#' @examples
#' \dontrun{
#' # From vignette:
#' calculateFTVA(grasslandEVselect, formula("~ prbygall + geoberg + lcucor1 +
#' tertpi09 + geolmja1"))
#' }
#'
#' @export


calculateFTVA <- function(selectedEV, formula = NULL) {
  tab <- selectedEV$selection
  if (is.character(formula) && !grepl(pattern = '~', x = formula)) {
    formula <- paste("~", formula)
  }
  if (is.null(formula)) {
    formula <- selectedEV$selectedmodel$formula
  }
  formula <- stats::as.formula(formula)
  trms <- labels(stats::terms(formula))
  evterms <- unique(gsub("_.*", "", trms))
  trailterms <- sapply(tab$variables, FUN=strsplit, split = " + ", fixed=TRUE)
  equal <- sapply(trailterms, setequal, y = evterms)
  if (!any(equal)) {
    stop("The specified formula does not correspond to any model in the selection table",
         call. = FALSE)
  }
  modrow <- which(equal)
  modround <- tab[modrow, match("round", names(tab))]
  roundtab <- tab[c(match(seq_len(modround - 1), tab$round), modrow),]

  Vt <- roundtab[nrow(roundtab), match("Dsq", names(roundtab))]
  previous <- c(0, roundtab$Dsq[-nrow(roundtab)])
  df <- data.frame(variable = unlist(unname(trailterms[modrow])),
                   FTVA = round(((roundtab$Dsq - previous) / Vt), digits = 3))
  return(df)
}
