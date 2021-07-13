## This is the tenth resubmission

* Updates since previous submission:
  * New function `fcsprocessor()` to convert a collection of Flow Cytometry Standard (FCS) files to a data frame the `gateR` package can read
  * Added new contributor, Surabhi Nair, who created the `fcsprocessor()` function
  * Now `flowCore` and `tools` are in Depends
  * ...
  
* Documentation for `pval_correct()` references doi <https://doi.org/10.2307/2283989> and <https://doi.org/10.1214/aop/1176996176> that throw NOTES in win-builder, Fedora Linux, and Ubuntu Linux but these are valid URLs
  
## Test environments
* local OS X install, R 4.1.0
* win-builder, (devel, release, oldrelease)
* Rhub
  * Fedora Linux, R-devel, clang, gfortran
  * Ubuntu Linux 16.04 LTS, R-release, GCC
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Oracle Solaris 10, x86, 32 bit, R-release

## R CMD check results
0 errors | 0 warnings | 0 notes

## Submitted by Maintainer
