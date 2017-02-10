<!-- README.md is generated from README.Rmd. Please edit that file -->
MIAmaxent
=========

### Description

Training, selecting, and evaluating maximum entropy (Maxent) distribution models. This package provides tools for user-controlled transformation of explanatory variables, selection of variables by nested model comparison, and flexible model evaluation and projection. The methods implemented here are based on the strict maximum likelihood interpretation of maximum entropy modelling (Halvorsen, 2013, Halvorsen et al., 2015).

The predecessor to this package is the MIA Toolbox, which is described in detail in Mazzoni et al. (2015).

**MIAmaxent** is built around the highly popular MaxEnt distribution modeling program (Phillips et al., 2006), but provides an alternative methodology for training, selecting, and using models. The major advantage in this alternative methodology is greater user control -- in variable transformations, in variable selection, and in model output. Comparisons also suggest that this methodology results in simpler models with equally good predictive ability, and reduces the risk of overfitting (Halvorsen et al., 2016).

### Installation

Install the release version from CRAN:

``` r
install.packages("MIAmaxent")
```

Or the development version from github

``` r
# install.packages('devtools')
# install.packages('R.rsp')
devtools::install_github("julienvollering/MIAmaxent", build_vignettes=TRUE)
```

### System Requirements

The maximum entropy algorithm utilized in this package is provided by the MaxEnt Java program (Phillips et al., 2006). This software is freely available, but may not be distributed further. Therefore, you must download the MaxEnt program (v3.3.3k) from <https://www.cs.princeton.edu/~schapire/maxent/>, and place the 'maxent.jar' file in the 'java' folder of this package. This folder can be located by the following R command: system.file("java", package = "MIAmaxent").

You must have the Java Runtime Environment (JRE) installed on your computer for the MaxEnt program to function. You can check if you have Java installed, and download it if necessary, at <http://java.com/download>.

### User Workflow

This diagram outlines a common workflow for users of this package. Functions are shown in red.

![](https://raw.githubusercontent.com/julienvollering/MIAmaxent/master/man/figures/workflow-flowchart.png)

### References

Halvorsen, R. (2013) A strict maximum likelihood explanation of MaxEnt, and some implications for distribution modelling. Sommerfeltia, 36, 1-132.

Halvorsen, R., Mazzoni, S., Bryn, A. & Bakkestuen, V. (2015) Opportunities for improved distribution modelling practice via a strict maximum likelihood interpretation of MaxEnt. Ecography, 38, 172-183.

Halvorsen, R., Mazzoni, S., Dirksen, J.W., Næsset, E., Gobakken, T. & Ohlson, M. (2016) How important are choice of model selection method and spatial autocorrelation of presence data for distribution modelling by MaxEnt? Ecological Modelling, 328, 108-118.

Mazzoni, S., Halvorsen, R. & Bakkestuen, V. (2015) MIAT: Modular R-wrappers for flexible implementation of MaxEnt distribution modelling. Ecological Informatics, 30, 215-221.

Phillips, S.J., Anderson, R.P. & Schapire, R.E. (2006) Maximum entropy modeling of species geographic distributions. Ecological Modelling, 190, 231-259.

------------------------------------------------------------------------

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/MIAmaxent)](https://cran.r-project.org/package=MIAmaxent) [![Travis-CI Build Status](https://travis-ci.org/julienvollering/MIAmaxent.svg?branch=master)](https://travis-ci.org/julienvollering/MIAmaxent) [![CRAN download rate](http://cranlogs.r-pkg.org/badges/MIAmaxent)](http://cran.rstudio.com/web/packages/MIAmaxent/index.html)
