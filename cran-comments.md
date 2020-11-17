## This is the third resubmission

* Actions taken regarding feedback from Prof. Brian Ripley:
  * Removed 'ncdFlow', 'flowWorkspaceData', and 'knitr' from Suggests
  * Created a random data set 'randCyto' and all documentation
  * Updated examples and testthat to use 'randCyto' data
  * Updated vignette with clearer language

## Test environments
* local OS X install, R 3.6.3
* win-builder, (devel, oldrelease, release)
* Rhub
  * Fedora Linux, R-devel, clang, gfortran
  * Ubuntu Linux 16.04 LTS, R-release, GCC
  * Windows Server 2008 R2 SP1, R-devel, 32⁄64 bit
  * Debian Linux, R-devel, GCC

## R CMD check results
0 errors | 0 warnings | 0 notes

* Duration
  * local OS X install, R 3.6.3: 98 seconds
  * win-builder, devel: 337 seconds
  * win-builder, oldrelease: 302 seconds
  * win-builder, release: 345 seconds

## Submitted by Maintainer
