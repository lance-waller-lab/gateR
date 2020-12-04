# gateR (development version)

# gateR v0.1.5.9000
  * Added arguments 'save_gate', 'name_gate', and 'path_gate' in lotrrs(), rrs(), and gating() to save plots as PNG files as output
  * Renamed 'doplot' argument as 'plot_gate' for consistency with new plotting arguments

# gateR v0.1.4
  * Added documentation to lotrrs(), rrs(), and gating() about the levels of condition(s)
  * Fixed bug in lotrrs() that was mislabeling numerator and denominator levels of second condition
  * Added parameters 'c1n' and 'c2n' in lotrrs(), rrs(), and gating() to specify the numerator level

# gateR v0.1.3
  * Removed 'ncdFlow', 'flowWorkspaceData', and 'knitr' from Suggests
  * Created a random data set 'randCyto' and all documentation
  * Updated examples and testthat to use 'randCyto' data
  * Updated vignette with clearer language
