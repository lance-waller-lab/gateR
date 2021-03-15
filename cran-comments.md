## This is the eighth resubmission

* Updates since previous submission:
  * Updated `spatstat` package to new subsetted packages based on feedback from the Spatstat Team (Adrian Baddeley and Ege Rubak). Now `spatstat.geom`, `spatstat.core`, `spatstat.linnet`, and `spatstat (>= 2.0-0)` are in Depends
  * Fixed check for `vars` in `dat` within the `gating()` function
  
* Documentation for `pval_correct()` references a doi <https://doi.org/10.2307/2283989> that throws a NOTE in win-builder but no other environment
  
## Test environments
* local OS X install, R 4.0.4
* win-builder, (devel, release, oldrelease)
* Rhub
  * Fedora Linux, R-devel, clang, gfortran
  * Ubuntu Linux 16.04 LTS, R-release, GCC
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Oracle Solaris 10, x86, 32 bit, R-release

## R CMD check results
0 errors | 0 warnings | 0 notes

## Submitted by Maintainer
