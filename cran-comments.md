## This is the fourth resubmission

* Actions taken since previous resubmission:
  * Added documentation to lotrrs(), rrs(), and gating() about the levels of condition(s)
  * Fixed bug in lotrrs() that was mislabeling numerator and denominator levels of second condition
  * Added parameters 'c1n' and 'c2n' in lotrrs(), rrs(), and gating() to specify the numerator level
  
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
