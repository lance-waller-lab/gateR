# code to prepare `randCyto` dataset

set.seed(1234)

# ------------------ #
# Necessary packages #
# ------------------ #

library(dplyr)
library(flowWorkspaceData)
library(ncdfFlow)
library(stats)
library(usethis)

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
## Arcsinh Transformation for all markers
## Remove cells with NA and Inf values

# First sample
obs_dat1 <- data.frame("id" = seq(1, nrow(fr1@exprs), 1),
                       "g1" = rep(1, nrow(fr1@exprs)),
                       "arcsinh_CD4" = asinh(fr1@exprs[ , 5]),
                       "arcsinh_CD38" = asinh(fr1@exprs[ , 6]),
                       "arcsinh_CD8" = asinh(fr1@exprs[ , 7]),
                       "arcsinh_CD3" = asinh(fr1@exprs[ , 8]))
# Second sample
obs_dat2 <- data.frame("id" = seq(1, nrow(fr2@exprs), 1),
                       "g1" = rep(2, nrow(fr2@exprs)),
                       "arcsinh_CD4" = asinh(fr2@exprs[ , 5]),
                       "arcsinh_CD38" = asinh(fr2@exprs[ , 6]),
                       "arcsinh_CD8" = asinh(fr2@exprs[ , 7]),
                       "arcsinh_CD3" = asinh(fr2@exprs[ , 8]))

# Full set
obs_dat <- rbind(obs_dat1, obs_dat2)
obs_dat <- obs_dat[complete.cases(obs_dat), ] # remove NAs
obs_dat <- obs_dat[is.finite(rowSums(obs_dat)), ] # remove Infs
obs_dat$g1 <- as.factor(obs_dat$g1) # set "g1" as binary factor

## Create a second condition (randomly split the data)
## In practice, use data with a measured second condition
g2 <- stats::rbinom(nrow(obs_dat), 1, 0.5)
obs_dat$g2 <- as.factor(g2)
obs_dat <- obs_dat[ , c(1:2,7,3:6)]

# Reduce size for CRAN exportation
randCyto <- dplyr::sample_frac(obs_dat, size = 0.1) # saved as `randCyto` data

# Export
usethis::use_data(randCyto, overwrite = TRUE)
