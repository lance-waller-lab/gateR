## This is the first resubmission

* Actions taken regarding feedback from CRAN teams' auto-check service
  * Reduced overall checktime by changing example in gating() by changing from \donttest{} to if (interactive()) {}
  * Reduced overall checktime by streamlining the number of "[it] works" tests in testthat for gating(), lotrrs(), and rrs()
  * Note about "Possibly mis-spelled words in DESCRIPTION" can be safely ignored

## Test environments
* local OS X install, R 3.6.3
* win-builder, (devel, oldrelease, release)
* Rhub
  * Fedora Linux, R-devel, clang, gfortran
  * Ubuntu Linux 16.04 LTS, R-release, GCC
  * Windows Server 2008 R2 SP1, R-devel, 32⁄64 bit
  * Debian Linux, R-devel, GCC
  
There was 1 NOTE:

* Possibly mis-spelled words in DESCRIPTION:
  * Bithell (38:63, 40:24)
  * Cytometry (3:57)
  * al (37:57)
  * et (37:54)
  * immunologically (32:19)

## R CMD check results
0 errors | 0 warnings | 0 notes

* Duration
  * local OS X install, R 3.6.3: 99 seconds
  * win-builder, devel: 584 seconds (under 10 minutes)
  * win-builder, oldrelease: 530 seconds
  * win-builder, release: 620 seconds

## Submitted by Maintainer
