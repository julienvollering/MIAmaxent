This patch is submitted to address an ERROR and a NOTE found during R CMD check on CRAN. The CRAN team requested: "Please correct before 2025-12-15 to safely retain your package on CRAN."

## R CMD check results

0 errors | 1 warning | 1 note

❯ checking PDF version of manual ... WARNING
  LaTeX errors when creating PDF version.
  This typically indicates Rd problems.
  LaTeX errors found:

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    'MIAmaxent-manual.tex'

I believe the warning and note to be false positives. They do not appear on the r-devel check, only my local install. Even there, the PDF manual compiles successfully (23 pages, 185KB) with no actual LaTeX errors. The log shows only cosmetic "Overfull \hbox" warnings for code lines exceeding margins. The "PDF version of manual without index" check passes without issues.

## R CMD check environments

* local windows 11, r 4.5.2
* winbuilder, r-devel

## revdepcheck results

There are currently no reverse dependencies.
