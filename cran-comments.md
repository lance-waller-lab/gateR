## This is the eleventh resubmission

* First resubmission after CRAN archived `gateR` on 2022-08-08 because a dependency `SpatialPack` and its dependency `fastmatrix` were archived. Since the warning from CRAN 2022-07-20, the `gateR` maintainer worked with the maintainer of both packages. The dependencies have returned to CRAN and this maintainer is now resubmitting `gateR`.

* Actions taken regarding feedback from CRAN teams' auto-check service:
  * Replaced `if()` conditions comparing `class()` to string with `inherits()`

* Updates since previous submission:
  * Updated package URL and BugReports to renamed GitHub account "lance-waller-lab" (previously "Waller-SUSAN")
  * `tools` is no longer Imports
  * `utils` is now Suggests because "zzz.R" calls the `packageDescription()` function
  * `ncdfFlow`, `flowWorkspaceData` are no longer Suggests because "Package suggested but not available for checking" in the following CRAN environments:
    * r-devel-linux-x86_64-fedora-clang
    * r-devel-linux-x86_64-fedora-gcc
    * r-devel-windows-x86_64-new-TK
    * r-release-linux-x86_64
    * r-release-macos-x86_64
    * r-oldrel-macos-x86_64
  * Added CITATION file
  * Fixed typos in documentation throughout

* Documentation for DESCRIPTION references the following DOIs that throw a NOTE in win-builder, Fedora Linux, and Ubuntu Linux but are valid URLs:
  * <https://doi.org/10.1002/sim.4780090616>
  * <https://doi.org/10.1002/sim.4780101112>
  * <https://doi.org/10.1002/sim.7577> 

* Documentation for `pval_correct()` references the following DOIs that throw a NOTE in win-builder, Fedora Linux, and Ubuntu Linux but are valid URLs:
  * <https://doi.org/10.2307/2283989>
  * <https://doi.org/10.1214/aop/1176996176>
  * <https://doi.org/10.1038/jcbfm.1991.122>
  * <https://doi.org/10.1111/j.2517-6161.1995.tb02031.x>
  
## Test environments
* local OS X install, R 4.2.1
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
