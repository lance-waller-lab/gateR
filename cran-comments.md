## This is the seventh resubmission

* Updates since previous submission:
  * Updated 'spatstat' dependency to include 'spatstat.core' and 'spatstat.geom' because 'spatstat' was recently subsetted
  * Added `p_correct` arugment to `lrren()` and `perlrren()` which calls a new, internal function `pval_correct()` that calculates seven types of corrections for multiple testing
  * Import 'lifecycle' package to document deprecated arguments `doplot` and `verbose` in `lotrrs()`, `rrs()`, and `gating()` functions
  
* Documentation for `pval_correct()` references a doi <https://doi.org/10.2307/2283989> that throws a NOTE in win-builder but no other environment
  
## Test environments
* local OS X install, R 3.6.3
* win-builder, (devel, oldrelease, release)
* Rhub
  * Fedora Linux, R-devel, clang, gfortran
  * Ubuntu Linux 16.04 LTS, R-release, GCC
  * Windows Server 2008 R2 SP1, R-devel, 32⁄64 bit

## R CMD check results
0 errors | 0 warnings | 0 notes

## Submitted by Maintainer
