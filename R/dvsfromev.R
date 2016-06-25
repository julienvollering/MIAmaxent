#' Produces DVs from a given EV by different transformations.
#'
#' In MIAT, "D" transformation is only performed if the optimum occurs in the
#' middle 80\% of the EV range. \code{dvsfromev} does not currently specify this
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


.dvsfromev <- function(df, transformtype, allsplines, dir, jarpath) {

  rv <- df[, 1]
  ev <- df[, 2]
  evname <- colnames(df)[2]
  storage <- list()

  if (class(ev) %in% c("numeric", "integer")) {

    if (any(c("HF", "HR", "T") %in% transformtype) && allsplines == F) {
      evdir <- .dirpath.create(dir, evname)
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
        hfdir <- .dirpath.create(evdir, "HF")
        message(paste0("Selecting forward hinge transformations of ", evname))
        dvs <- lapply(splall, function(x) {x(ev)})
        names(dvs) <- gsub("_transf", "", names(splall))
        selected <- .splselect(rv, dvs, hfdir, jarpath)
        if (length(selected) > 0) {
          storage <- c(storage, splall[paste0(selected, "_transf")])
        }
      }
    }

    if ("HR" %in% transformtype) {
      splall <- list()
      knots <- 20
      for (i in 1:knots) {
        k <- (2 * i - 1) / (2 * knots)
        tfunction <- .transfSpline(ev, k, type = "HR")
        splall[[paste0(evname, "_HR", i, "_transf")]] <- tfunction
      }
      if (allsplines == T) {
        storage <- c(storage, splall)
      } else {
        hrdir <- .dirpath.create(evdir, "HR")
        message(paste0("Selecting reverse hinge transformations of ", evname))
        dvs <- lapply(splall, function(x) {x(ev)})
        names(dvs) <- gsub("_transf", "", names(splall))
        selected <- .splselect(rv, dvs, hrdir, jarpath)
        if (length(selected) > 0) {
          storage <- c(storage, splall[paste0(selected, "_transf")])
        }
      }
    }

    if ("T" %in% transformtype) {
      splall <- list()
      knots <- 20
      for (i in 1:knots) {
        k <- (2 * i - 1) / (2 * knots)
        tfunction <- .transfSpline(ev, k, type = "T")
        splall[[paste0(evname, "_T", i, "_transf")]] <- tfunction
      }
      if (allsplines == T) {
        storage <- c(storage, splall)
      } else {
        tdir <- .dirpath.create(evdir, "T")
        message(paste0("Selecting threshold transformations of ", evname))
        dvs <- lapply(splall, function(x) {x(ev)})
        names(dvs) <- gsub("_transf", "", names(splall))
        selected <- .splselect(rv, dvs, tdir, jarpath)
        if (length(selected) > 0) {
          storage <- c(storage, splall[paste0(selected, "_transf")])
        }
      }
    }
  }

  if (class(ev) %in% c("factor", "character")) {
    ev <- as.factor(ev)

    if ("B" %in% transformtype) {
      for (i in levels(ev)) {
        tfunction <- .transfB(ev, i)
        storage[[paste0(evname, "_B", i, "_transf")]] <- tfunction
      }
    }
  }

  evdv <- data.frame(lapply(storage, function(x) {x(ev)}))
  colnames(evdv) <- gsub("_transf", "", names(storage))
  return(list("storage" = storage, "evdv" = evdv))
}
