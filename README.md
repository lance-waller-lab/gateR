gateR: Gating strategy for mass cytometry using kernel density estimation <img src="man/figures/gateR.png" width="120" align="right" />
===================================================

<h2 id="overview">

Overview

</h2>

The `gateR` package is a suite of `R` functions to identify significant spatial clustering of mass cytometry data used in immunological investigations. For a two group comparison we detect clusters using the kernel-based spatial relative risk function that is estimated using the [sparr](https://CRAN.R-project.org/package=sparr) package. The tests are conducted in two-dimenstional space comprised of two flourescent markers. 

Examples for a single condition:

1. Disease case v. healthy control
2. Time 1 v. Time 0

For a two group comparison for two conditions we estimate two relative risk surfaces for one condition and then a ratio of the relative risks. For example:

1. Estimate a relative risk surface for
A. Time 1: Disease case v. healthy control
B. Time 0: Disease case v. healthy control
2. Estimate  relative risk surface for Time 1 v. Time 2

Within areas where the relative risk exceeds an asymptotic normal assumption, the `gateR` package has functionality to examine the features of these cells. Basic visualization is also supported. 

<h2 id="install">

Installation

</h2>

To install the release version from CRAN:

    install.packages("gateR")

To install the development version from GitHub:

    devtools::install_github("Waller-SUSAN/gateR")

<h2 id="available-functions">

Available functions

</h2>

<table>
<colgroup>
<col width="30%" />
<col width="70%" />
</colgroup>
<thead>
<tr class="header">
<th>Function</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<td><code>gating</code></td>
<td>Main function. Conduct a gating strategy for mass cytometry data.</td>
</tr>
<td><code>rrs</code></td>
<td>Called within `gating`, one condition comparison.</td>
</tr>
<td><code>lotrrs</code></td>
<td>Called within `gating`, two condition comparison. </td>
</tr>
</tbody>
<table>

## Usage
``` r
# ------------------ #
# Necessary packages #
# ------------------ #

library(flowWorkspaceData)
library(ncdfFlow)
library(stats)

# ---------------- #
# Data preparation #
# ---------------- #

# Use 'extdata' from the {flowWorkspaceData} package
flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
fcsFiles <- list.files(pattern = "CytoTrol", flowDataPath, full = TRUE)
ncfs  <- ncdfFlow::read.ncdfFlowSet(fcsFiles)
fr1 <- ncfs[[1]]
fr2 <- ncfs[[2]]

## Comparison of two samples (single condition) "g1"
## Two gates (four markers) "CD4", "CD38", "CD8", and "CD3"
## Log10 Transformation for all markers
## Remove cells with NA and Inf values

# First sample
obs_dat1 <- data.frame("id" = seq(1, nrow(fr1@exprs), 1),
                       "g1" = rep(1, nrow(fr1@exprs)),
                       "log10_CD4" = log(fr1@exprs[ , 5], 10),
                       "log10_CD38" = log(fr1@exprs[ , 6], 10),
                       "log10_CD8" = log(fr1@exprs[ , 7], 10),
                       "log10_CD3" = log(fr1@exprs[ , 8], 10))
# Second sample
obs_dat2 <- data.frame("id" = seq(1, nrow(fr2@exprs), 1),
                       "g1" = rep(2, nrow(fr2@exprs)),
                       "log10_CD4" = log(fr2@exprs[ , 5], 10),
                       "log10_CD38" = log(fr2@exprs[ , 6], 10),
                       "log10_CD8" = log(fr2@exprs[ , 7], 10),
                       "log10_CD3" = log(fr2@exprs[ , 8], 10))
# Full set
obs_dat <- rbind(obs_dat1, obs_dat2)
obs_dat <- obs_dat[complete.cases(obs_dat), ] # remove NAs
obs_dat <- obs_dat[is.finite(rowSums(obs_dat)),] # remove Infs
obs_dat$g1 <- as.factor(obs_dat$g1) # set "g1" as binary factor

# --------- #
# Run gateR #
# --------- #

# Single Condition
## A p-value uncorrected for multiple testing
test_gate <- gateR::gating(dat = obs_dat,
                           vars = c("log10_CD4", "log10_CD38", "log10_CD8", "log10_CD3"),
                           n_condition = 1,
                           doplot = TRUE,
                           p_cor = "none")

# -------------------- #
# Post-gate assessment #
# -------------------- #

# Density of log10 CD4 post-gating
graphics::plot(stats::density(test_gate$obs[test_gate$obs$g1 == 1, 3]),
               main = "log10 CD4",
               lty = 2)
graphics::lines(stats::density(test_gate$obs[test_gate$obs$g1 == 2, 3]),
                lty = 3)
graphics::legend("topright",
                 legend = c("Sample 1", "Sample 2"),
                 lty = c(2, 3),
                 bty = "n")

# --------------------------------------------- #
# Perform a single gate without data extraction #
# --------------------------------------------- #

# Single condition
## A p-value uncorrected for multiple testing
test_rrs <- gateR::rrs(dat = obs_dat,
                       p_cor = "none")

# Two conditions
## Create a second condition (randomly split the data)
## In practice, use data with a measured second condition
g2 <- stats::rbinom(nrow(obs_dat), 1, 0.5)
obs_dat$g2 <- as.factor(g2)
obs_dat <- obs_dat[ , c(1:2,7,3:6)]

## A p-value uncorrected for multiple testing
test_lotrrs <- gateR::lotrrs(dat = obs_dat,
                             p_cor = "none")
```
![](man/figures/gate1.png)

![](man/figures/gate2.png)

![](man/figures/postgate.png)