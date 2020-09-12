install.packages(c("devtools", "roxygen2", "testthat", "knitr"))
library(devtools)
devtools::has_devel()
library(roxygen2)
library(testthat)

devtools::load_all()

getOption("gateR")

# Convert roxygen components to .Rd files
devtools::document()
?gateR

# Create Vignette
install()
build()

# Testing
use_testthat()
use_test()
test()

# NAMESPACE
document()
install()

# Check
check()

# Ignore .R files from /build directory
usethis::use_build_ignore(c("build"))
usethis::use_build_ignore(c("figures"))

# Example in README
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
## Two gates (Four markers) "CD4", "CD38", "CD8", and "CD3"
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
                             doplot = FALSE,
                             p_cor = "none")


# Example in lotrrs()
  library(flowWorkspaceData)
  library(ncdfFlow)
  library(stats)

# Use 'extdata' from the {flowWorkspaceData} package
  flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
  fcsFiles <- list.files(pattern = "CytoTrol", flowDataPath, full = TRUE)
  ncfs  <- ncdfFlow::read.ncdfFlowSet(fcsFiles)
  fr1 <- ncfs[[1]]
  fr2 <- ncfs[[2]]

## Comparison of two samples at two time points (two conditions) "g1" and "g2"
## (Create a random binary variable for "g2")
## One gates (Two markers) "CD4", "CD38"
## Log10 Transformation for both markers
## Remove cells with NA and Inf values

# First sample
  obs_dat1 <- data.frame("id" = seq(1, nrow(fr1@exprs), 1),
                         "g1" = rep(1, nrow(fr1@exprs)),
                         "g2" = stats::rbinom(nrow(fr1@exprs), 1, 0.5),
                         "log10_CD4" = log(fr1@exprs[ , 5], 10),
                         "log10_CD38" = log(fr1@exprs[ , 6], 10))
# Second sample
  obs_dat2 <- data.frame("id" = seq(1, nrow(fr2@exprs), 1),
                         "g1" = rep(2, nrow(fr2@exprs)),
                         "g2" = stats::rbinom(nrow(fr2@exprs), 1, 0.5),
                         "log10_CD4" = log(fr2@exprs[ , 5], 10),
                         "log10_CD38" = log(fr2@exprs[ , 6], 10))
# Full set
  obs_dat <- rbind(obs_dat1, obs_dat2)
  obs_dat <- obs_dat[complete.cases(obs_dat), ] # remove NAs
  obs_dat <- obs_dat[is.finite(rowSums(obs_dat)),] # remove Infs
  obs_dat$g1 <- as.factor(obs_dat$g1) # set "g1" as binary factor
  obs_dat$g2 <- as.factor(obs_dat$g2) # set "g2" as binary factor

# Run lotrrs() function
  test_lotrrs <- lotrrs(dat = obs_dat, p_cor = "none")



plot(out_gate$rrs[[2]])
