## Patch submission
I am submitting this patch to correct an error in the CRAN check results for r-oldrel-windows. 

## Test environments
* local Windows install, R 3.3.1
* ubuntu 12.04 (on travis-ci), R 3.3.1
* os x 10.9, (on travis-ci), R 3.3.1
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* 1 note under R devel: 
    "checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Julien Vollering <julien.vollering@hisf.no>'
    Days since last update: 2
    ...
    Possibly mis-spelled words in DESCRIPTION:
    Maxent (3:8, 11:6)"
* "Maxent" is not a mis-spelling.

## Reverse dependencies

There are currently no reverse dependencies.

---
