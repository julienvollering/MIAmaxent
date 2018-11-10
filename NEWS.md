# MIAmaxent 1.0.0

## Major changes
* Model fitting implemented as infinitely-weighted logistic regression, so that all computation can be done natively in R (maxent.jar no longer required).
* Implements choice of algorithm: "maxent" for maximum entropy or "LR" for standard logistic regression (binomial GLM).
* No files written to system unless write = TRUE
* Choice of Chi-squared or F-test in nested model comparison
* More consistency in arguments across top-level functions
* Selection trail tables simplified and clarified
## Minor changes
* increased flexibility in graphics arguments passed to plotting functions
* quiet option added to top-level functions performing selection
* readData() automatically removes duplicates when two or more presences/absences fall in the same cell
* readData() discards presence locations with missing EV data
* formula argument to selectEV() function, to specifiy starting point for selection
* plotFOP() smoother changed to loess from exponentially weighted moving average
* plotFOP() plots data density behind FOP values
* plotFOP() returns plot data invisibly
* plotResp() and plotResp2() take identical arguments, the first of which is a model object
* projectModel() takes data in data.frame or raster classes, and plots output spatially in the case of the latter
* trainmax argument removed from selectDVforEV() and selectEV()
* testAUC() plotting optional

# MIAmaxent 0.4.0

* Model ranking within selection rounds based on p-value and then F-statistic (tiebreaker), rather than simply F-statistic
* Directories specified by 'dir' argument are created if they do not already exist.
* Existing results in directories specified by 'dir' argument are overwritten, if desired.
* Fixed bug in selectEV that occurred when the last round of model selection before interaction terms did not result in a significant variable.
* Unnecessary dependency on Hmisc removed.

# MIAmaxent 0.3.7

* Removed version minimums for dependencies which are default packages, to allow r-oldrel binary.
* Changed names of toy data used in examples for better organization.

# MIAmaxent 0.3.6

* First public release.



