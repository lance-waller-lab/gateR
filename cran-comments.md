## This is the sixteenth resubmission

* Updates since previous submission:
  * Renamed 'package.R' to 'gateR-package.R' after Roxygen (>=7.3.0) update
  * Fixed Rd `\link{}` targets missing package within 'gateR-package.Rd'
  * Updated vignette with updated package version and dependencies

## Notes
* The win-builder oldrelease throws a NOTE that "Author field differs from that derived from Authors@R". The behavior is OK because ORCID has different formatting but same information
  
## Test environments
* local Windows install, R 4.5.1
* win-builder, (devel, release, oldrelease)
* R-CMD-check on GitHub
  * macos-latest (release)
  * windows-latest (release)
  * ubuntu-latest (devel)
  * ubuntu-latest (release)
  * ubuntu-latest (oldrel-1)
* Rhub v2
  * macos-15 on GitHub, ASAN + UBSAN on macOS (`m1-san`)
  * macos-13 on GitHub(`macos`)
  * Fedora Linux 40 (Container Image) (`gcc-asan`)
  * Ubuntu 22.04.5 LTS (`ubuntu-clang`)
  * Ubuntu 22.04.5 LTS (`ubuntu-gcc12`)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Submitted by Maintainer
