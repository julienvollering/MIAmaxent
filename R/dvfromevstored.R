#' produces dvs from a given ev by different transformations.
#'
#' In MIAT, "D" transformation is only performed if the optimum occurs in the
#' middle 80\% of the EV range. \code{dvfromev} does not currently specify this
#' condition.
#'
#' @param df Dataframe with 2 columns: response variable and explanatory
#'   variable (in that order). Column names are used as identifiers. The
#'   response variable represents presence or background, coded as: 1/NA. The
#'   explanatory variable may be continuous or categorical.
#' @param transformtype Set of transformation types to be used.
#' @param allsplines Logical. Keep all spline transformations.
#' @param dir Directory to which Maxent runs of spline transformations are
#'   written
#' @param jarpath Pathway to maxent.jar
#'
#' @return Dataframe with one column for each DV.


# IMPORTANT: This file is for testing purposes only (e.g. transformation
# storage)! Editing for inclusion in the package must proceed in the _dvfromev
# file as it contains the latest bug fixes.

.dvfromev <- function(df, transformtype, allsplines, dir, jarpath) {

  rv <- df[,1]
  ev <- df[,2]
  evname <- colnames(df)[2]
  storage <- list()

  if (class(ev) == "numeric" || class(ev) == "integer") {

    if (any(c("HF", "HR", "T") %in% transformtype) && allsplines == F) {
      evdir <- paste(dir, "\\", evname, sep="")
      dir.create(evdir)
    }

    if ("L" %in% transformtype) {
      tfunction <- .transfL(ev)
      storage[[paste0(evname, "_L_transf")]] <- tfunction
    }

    if ("M" %in% transformtype) {
      tfunction <- .transfM(ev)
      storage[[paste0(evname, "_M_transf")]] <- tfunction
    }

    if ("D" %in% transformtype) {
      for (i in c(0.5, 1, 2)) {
        tfunction <- .transfD(ev, rv, i)
        storage[[paste0(evname, "_D", i, "_transf")]] <- tfunction
      }
    }

    if ("HF" %in% transformtype) {
      splall <- list()
      knots <- 20
      for (i in 1:knots) {
        k <- (2 * i - 1) / (2 * knots)
        tfunction <- .transfSpline(ev, k, type = "HF")
        splall[[paste0(evname, "_HF", i, "_transf")]] <- tfunction
      }
      if (allsplines == T) {
        storage <- c(storage, splall)
      } else {
        hfdir <- paste(evdir, "\\HF", sep="")
        dir.create(hfdir)
        message(paste0("Selecting forward hinge transformations of ", evname))
        dvs <- lapply(splall, function(x) {x(ev)})
        selected <- altrMaxent:::.splselect(rv, dvs, hfdir, jarpath)
        storage <- c(storage, splall[selected])
      }
    }
  }

  if (class(ev) == "factor" || class(ev) == "character") {
    if ("B" %in% transformtype) {
      B <- stats::model.matrix( ~ ev - 1, data=df )
      colnames(B) <- paste(evname, "_B", levels(ev), sep="")
      evdv <- cbind(evdv, B)
    }
  }

  evdv <- data.frame(lapply(storage, function(x) {x(ev)}))
  colnames(evdv) <- gsub("_transf", "", names(storage))
  return(evdv)
}
