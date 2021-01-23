# gateR (development version)

# gateR v0.1.7
  * Updated {spatstat} dependency to new packages based on feedback from the Spatstat Team (Adrian Baddeley and Ege Rubak). {spatstat.geom} replaces {spatstat} in Imports
  * Added additional multiple testing corrections, including False Discovery Rate, spatially dependent Sidak correction, and independent Sidak correction
  * Updated the calculation of the spatial correlogram in internal pval_correct() function from the correlog() function in the {pgrimess} package to the modified.ttest() function in the {SpatialPack} package

# gateR v0.1.6
  * Updated URLs in gateR-package.Rd
  * Updated year in DESCRIPTION

# gateR v0.1.5
  * Added arguments 'save_gate', 'name_gate', and 'path_gate' in lotrrs(), rrs(), and gating() to save plots as PNG files as output
  * Renamed 'doplot' argument as 'plot_gate' for consistency with new plotting arguments
  * Added a stop (and return no results) if no significant clusters detected during first gate
  * Removed 'verbose' argument in gating(), rrs(), and lotrrs()
  * Added try() error catches in rrs() and lotrrs() for 'c1n' and 'c2n' arguments
  * Changed the 'right' argument in cut() in pval_plot() to "TRUE" (the default)
  * Removed fullstop in error messages
  * Added a make.names() check for vars and colnames(dat) in gating()

# gateR v0.1.4
  * Added documentation to lotrrs(), rrs(), and gating() about the levels of condition(s)
  * Fixed bug in lotrrs() that was mislabeling numerator and denominator levels of second condition
  * Added parameters 'c1n' and 'c2n' in lotrrs(), rrs(), and gating() to specify the numerator level

# gateR v0.1.3
  * Removed 'ncdFlow', 'flowWorkspaceData', and 'knitr' from Suggests
  * Created a random data set 'randCyto' and all documentation
  * Updated examples and testthat to use 'randCyto' data
  * Updated vignette with clearer language
