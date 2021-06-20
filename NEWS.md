# gateR (development version)

# gateR v0.1.9
* Now `rlang` is in Depends. 
* Output for `gating()` now includes a diagnostic message saved as a character string for reference in a slot called `note`.
* A diagnostic message can be viewed running `rlang::last_error()$out$note` after the unsuccessful run of `gating()`.
* Removed redundant `@importFrom fields image.plot` in 'package.R'

# gateR v0.1.8
  * Updated `spatstat` package to new subsetted packages based on feedback from the Spatstat Team (Adrian Baddeley and Ege Rubak). Now `spatstat.geom`, `spatstat.core`, `spatstat.linnet`, and `spatstat (>= 2.0-0)` are in Depends
  * Fixed check for `vars` in `dat` within the `gating()` function

# gateR v0.1.7
  * Updated `spatstat` package to new subsetted packages based on feedback from the Spatstat Team (Adrian Baddeley and Ege Rubak). `spatstat.geom` package replaces `spatstat` package in Imports
  * Added additional multiple testing corrections, including False Discovery Rate, spatially dependent Sidak correction, independent Sidak correction, and two corrections based on Random Field Theory (Adler and Hasofer or Friston)
  * The latter two corrections required a new argument `bandw` to be added to `gating()`, `lotrrs()`, and `rrs()` functions to allow users to specify bandwidth for the kernel density estimation
  * Updated the calculation of the spatial correlogram in internal `pval_correct()` function from the `correlog()` function in the `pgrimess` package to the `modified.ttest()` function in the `SpatialPack` package
  * Imports `lifecycle` package to document deprecated arguments `doplot` and `verbose` in `lotrrs()`, `rrs()`, and `gating()` functions
  * In `gating()` function, creates a categorized `im` based on critical p-value, assigns that value to every point in a `ppp` object, and subsets points by category
  * Removed `maptools` and `sp` packages from Imports
  * Updated links in 'gateR-package.Rd' for package updates

# gateR v0.1.6
  * Updated URLs in 'gateR-package.Rd'
  * Updated year in DESCRIPTION

# gateR v0.1.5
  * Added arguments `save_gate`, `name_gate`, and `path_gate` in `lotrrs()`, `rrs()`, and `gating()` to save plots as PNG files as output
  * Renamed `doplot` argument as `plot_gate` for consistency with new plotting arguments
  * Added a stop (and return no results) if no significant clusters detected during first gate
  * Deprecate `verbose` argument in `gating()`, `rrs()`, and `lotrrs()`
  * Added `try()` error catches in `rrs()` and `lotrrs()` for `c1n` and `c2n` arguments
  * Changed the `right` argument in `cut()` in `pval_plot()` to "TRUE" (the default for cut)
  * Removed fullstop in error messages
  * Added a `make.names()` check for `vars` and `colnames(dat)` in `gating()`

# gateR v0.1.4
  * Added documentation to `lotrrs()`, `rrs()`, and `gating()` about the levels of condition(s)
  * Fixed bug in `lotrrs()` that was mislabeling numerator and denominator levels of second condition
  * Added parameters `c1n` and `c2n` in `lotrrs()`, `rrs()`, and `gating()` to specify the numerator level

# gateR v0.1.3
  * Removed `ncdFlow`, `flowWorkspaceData`, and `knitr` packages from Suggests
  * Created a random data set `randCyto` and all documentation
  * Updated examples and testthat to use `randCyto` data
  * Updated vignette with clearer language
