## Resubmission
This is a resubmission. In this version I have:
* Added examples that are not wrapped in \dontrun{}, and therefore are executed during checks. 3 functions rely on external java software that the user must install separately, so executed examples are not included for these.

## Test environments
* local Windows install, R 3.3.1
* ubuntu 12.04 (on travis-ci), R 3.3.1
* os x 10.9, (on travis-ci), R 3.3.1
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* 1 note under R devel: 
    "checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Julien Vollering <julien.vollering@hisf.no>'
    New submission
    ...
    Possibly mis-spelled words in DESCRIPTION:
    Maxent (3:8, 11:6)"
* "Maxent" is not a mis-spelling.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---
