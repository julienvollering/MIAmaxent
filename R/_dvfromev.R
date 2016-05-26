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


.dvfromev <- function(df, transformtype, allsplines, dir, jarpath) {

  rv <- df[,1]
  ev <- df[,2]
  evname <- colnames(df)[2]
  evdv <- data.frame(df[,2, drop=F])

  if (class(ev) == "numeric" || class(ev) == "integer") {

    if (any(c("HF", "HR", "T") %in% transformtype) && allsplines == F) {
      evdir <- paste(dir, "\\", evname, sep="")
      dir.create(evdir)
    }

    if ("L" %in% transformtype) {
      L <- data.frame((ev - range(ev)[1])/diff(range(ev)))
      colnames(L) <- paste(evname, "_L", sep = "")
      evdv <- cbind(evdv, L)
    }

    if ("M" %in% transformtype) {
      L <- (ev - range(ev)[1])/diff(range(ev))
      ZSk <- altrMaxent:::.scalex(L, altrMaxent:::.minskew(L)$c)
      M <- data.frame((ZSk - range(ZSk)[1])/diff(range(ZSk)))
      colnames(M) <- paste(evname, "_M", sep = "")
      evdv <- cbind(evdv, M)
    }

    if ("D" %in% transformtype) {
      L <- (ev - range(ev)[1])/diff(range(ev))
      opt <- altrMaxent:::.fopoptimum(data.frame(rv, L))
      devexp = c(0.5, 1, 2)
      D <- altrMaxent:::.devtransf(L, opt, devexp)
      colnames(D) <- paste(evname, "_D", devexp, sep = "")
      evdv <- cbind(evdv, D)
    }

    if ("HF" %in% transformtype) {
      L <- (ev - range(ev)[1])/diff(range(ev))
      knots <- 20
      hf <- altrMaxent:::.spltransf(L, knots, type = "HF")
      colnames(hf) <- paste(evname, "_HF", 1:knots, sep = "")
      if (allsplines == T) {
        HF <- hf
      } else {
        hfdir <- paste(evdir, "\\HF", sep="")
        dir.create(hfdir)
        message(paste("Selecting forward hinge transformations of ", evname,
          sep = ""))
        HF <- altrMaxent:::.splselect(rv, hf, hfdir, jarpath)
      }
      evdv <- cbind(evdv, HF)
    }

    if ("HR" %in% transformtype) {
      L <- (ev - range(ev)[1])/diff(range(ev))
      knots <- 20
      hr <- altrMaxent:::.spltransf(L, knots, type="HR")
      colnames(hr) <- paste(evname, "_HR", 1:knots, sep = "")
      if (allsplines == T) {
        HR <- hr
      } else {
        hrdir <- paste(evdir, "\\HR", sep="")
        dir.create(hrdir)
        message(paste("Selecting reverse hinge transformations of ", evname,
          sep = ""))
        HR <- altrMaxent:::.splselect(rv, hr, hrdir, jarpath)
      }
      evdv <- cbind(evdv, HR)
    }

    if ("T" %in% transformtype) {
      L <- (ev - range(ev)[1])/diff(range(ev))
      knots <- 20
      th <- altrMaxent:::.spltransf(L, knots, type = "T")
      colnames(th) <- paste(evname, "_T", 1:knots, sep = "")
      if (allsplines == T) {
        Th <- th
      } else {
        thdir <- paste(evdir, "\\T", sep="")
        dir.create(thdir)
        message(paste("Selecting threshold transformations of ", evname,
          sep = ""))
        Th <- altrMaxent:::.splselect(rv, th, thdir, jarpath)
      }
      evdv <- cbind(evdv, Th)
    }

  }


  if (class(ev) == "factor" || class(ev) == "character") {
    if ("B" %in% transformtype) {
      B <- stats::model.matrix( ~ ev - 1, data=df )
      colnames(B) <- paste(evname, "_B", levels(as.factor(ev)), sep="")
      evdv <- cbind(evdv, B)
    }
  }

  evdv <- evdv[,-1]
  return(evdv)
}
