## This is the tenth resubmission

* Updates since previous submission:
  * Updated dependencies `spatstat.core` and `spatstat.linnet` packages based on feedback from the Spatstat Team (Adrian Baddeley and Ege Rubak). All random generators in `spatstat.core` were moved to a new package `spatstat.random`
    * `spatstat.geom`, `spatstat.core`, `spatstat.linnet`, and `spatstat (>=2.0-0)` are no longer Depends
    * `spatstat.geom` is now Imports
  * Fixed annotation typos in the vignette. Removed packages no longer used in the vignette 
  * `dplyr`, `ncdfFlow`, `flowWorkspaceData`, and `usethis` now Suggests (for generating random data set `randCyto`)

* Documentation for `pval_correct()` references doi <https://doi.org/10.2307/2283989> and <https://doi.org/10.1214/aop/1176996176> that throw NOTES in win-builder, Fedora Linux, and Ubuntu Linux but these are valid URLs
  
## Test environments
* local OS X install, R 4.1.2
* win-builder, (devel, release, oldrelease)
* Rhub
  * Fedora Linux, R-devel, clang, gfortran
  * Ubuntu Linux 20.04.1 LTS, R-release, GCC
  * Windows Server 2022, R-devel, 64 bit
  * Windows Server 2008 R2 SP1, R-release, 32⁄64 bit
  * Oracle Solaris 10, x86, 32 bit, R-release
  * macOS 10.13.6 High Sierra, R-release, CRAN's setup

## R CMD check results
0 errors | 0 warnings | 0 notes

## Submitted by Maintainer
