context("gating")

###################
# gating testthat #
###################

# Generate testing data
## Use 'extdata' from the {flowWorkspaceData} package
flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
fcsFiles <- list.files(pattern = "CytoTrol", flowDataPath, full = TRUE)
ncfs  <- ncdfFlow::read.ncdfFlowSet(fcsFiles)
fr1 <- ncfs[[1]]
fr2 <- ncfs[[2]]

## Comparison of two samples (single condition) "g1"
## Two gates (Four markers) "CD4", "CD38", "CD8", and "CD3"
## Log10 Transformation for all markers
## Remove cells with NA and Inf values

## First sample
suppressMessages(suppressWarnings(
  obs_dat1 <- data.frame("id" = seq(1, nrow(fr1@exprs), 1),
                         "g1" = rep(1, nrow(fr1@exprs)),
                         "log10_CD4" = log(fr1@exprs[ , 5], 10),
                         "log10_CD38" = log(fr1@exprs[ , 6], 10),
                         "log10_CD8" = log(fr1@exprs[ , 7], 10),
                         "log10_CD3" = log(fr1@exprs[ , 8], 10))))
## Second sample
suppressMessages(suppressWarnings(
  obs_dat2 <- data.frame("id" = seq(1, nrow(fr2@exprs), 1),
                         "g1" = rep(2, nrow(fr2@exprs)),
                         "log10_CD4" = log(fr2@exprs[ , 5], 10),
                         "log10_CD38" = log(fr2@exprs[ , 6], 10),
                         "log10_CD8" = log(fr2@exprs[ , 7], 10),
                         "log10_CD3" = log(fr2@exprs[ , 8], 10))))

## Full set
obs_dat <- rbind(obs_dat1, obs_dat2)
obs_dat <- obs_dat[complete.cases(obs_dat), ] # remove NAs
obs_dat <- obs_dat[is.finite(rowSums(obs_dat)),] # remove Infs
obs_dat$g1 <- as.factor(obs_dat$g1) # set "g1" as binary factor

## Create a second condition (randomly split the data)
## In practice, use data with a measured second condition
g2 <- stats::rbinom(nrow(obs_dat), 1, 0.5)
obs_dat$g2 <- as.factor(g2)
obs_dat <- obs_dat[ , c(1:2,7,3:6)]

# Faulty datasets
fubar <- obs_dat
fubar$g1 <- as.numeric(fubar$g1)

# Tests

test_that("gating throws error with invalid arguments", {

  # Odd-numbered amount of vars
  expect_error(
    gating(dat = obs_dat,
           vars = c("log10_CD4"),
           n_condition = 1,
           p_correct = "none")
  )

  # A var not available
  expect_error(
    gating(dat = obs_dat,
           vars = c("log10_CD4", "log10_CD38",
                    "log10_CD8", "x"),
           n_condition = 1,
           p_correct = "none")
  )

  # Non-binary second feature
  expect_error(
    gating(dat = fubar,
           vars = c("log10_CD4", "log10_CD38",
                    "log10_CD8", "log10_CD3"),
           n_condition = 1,
           p_correct = "none")
  )

  # expect_error(
  #   gating(dat = fubar,
  #          vars = c("log10_CD4", "log10_CD38",
  #                   "log10_CD8", "log10_CD3"),
  #          n_condition = 2,
  #          p_correct = "none")
  # )

  # Non-binary third feature
  expect_error(
    gating(dat = fubar1,
           vars = c("log10_CD4", "log10_CD38",
                    "log10_CD8", "log10_CD3"),
           n_condition = 2,
           p_correct = "none")
  )

  # Incorrectly specified alpha
  expect_error(
    gating(dat = obs_dat,
           vars = c("log10_CD4", "log10_CD38",
                    "log10_CD8", "log10_CD3"),
           n_condition = 1,
           alpha = 0,
           p_correct = "none")
  )

  # Incorrectly specified p_correct
  expect_error(
    gating(dat = obs_dat,
           vars = c("log10_CD4", "log10_CD38",
                    "log10_CD8", "log10_CD3"),
           n_condition = 1,
           p_correct = "x")
  )

  # Incorrectly specified n_condition
  expect_error(
    gating(dat = obs_dat,
           vars = c("log10_CD4", "log10_CD38",
                    "log10_CD8", "log10_CD3"),
           n_condition = 3,
           p_correct = "none")
  )

  # Missing specified n_condition
  expect_error(
    gating(dat = obs_dat,
           vars = c("log10_CD4", "log10_CD38",
                    "log10_CD8", "log10_CD3"),
           p_correct = "none")
  )

}
)

test_that("gating works", {

  # expect_named(
  #     gating(dat = obs_dat,
  #            vars = c("log10_CD4", "log10_CD38"),
  #            n_condition = 1,
  #            p_correct = "none")
  # )
  # 
  # expect_named(
  #   gating(dat = obs_dat,
  #          vars = c("log10_CD4", "log10_CD38"),
  #          n_condition = 2,
  #          p_correct = "none")
  # )
  # 
  # expect_named(
  #   gating(dat = obs_dat,
  #          vars = c("log10_CD4", "log10_CD38"),
  #          n_condition = 1,
  #          p_correct = "uncorrelated")
  # )

  expect_named(
    gating(dat = obs_dat,
           vars = c("log10_CD4", "log10_CD38",
                    "log10_CD8", "log10_CD3"),
           n_condition = 1,
           p_correct = "none")
  )

  # expect_named(
  #   gating(dat = obs_dat,
  #          vars = c("log10_CD4", "log10_CD38",
  #                   "log10_CD8", "log10_CD3"),
  #          n_condition = 2,
  #          p_correct = "none")
  # )
  # 
  # expect_named(
  #   gating(dat = obs_dat,
  #          vars = c("log10_CD4", "log10_CD38"),
  #          n_condition = 1,
  #          alpha = 0.1,
  #          p_correct = "none")
  # )
  # 
  # expect_named(
  #   gating(dat = obs_dat,
  #          vars = c("log10_CD4", "log10_CD38"),
  #          n_condition = 1,
  #          alpha = 0.01,
  #          p_correct = "none")
  # )
  # 
  # expect_named(
  #   gating(dat = obs_dat,
  #          vars = c("log10_CD4", "log10_CD38"),
  #          n_condition = 1,
  #          numerator = FALSE,
  #          p_correct = "none")
  # )
  # 
  # expect_named(
  #   gating(dat = obs_dat,
  #          vars = c("log10_CD4", "log10_CD38"),
  #          n_condition = 2,
  #          doplot = TRUE,
  #          rcols = c("green", "yellow", "purple"),
  #          p_correct = "none")
  # )
  # 
  # expect_named(
  #   gating(dat = obs_dat,
  #          vars = c("log10_CD4", "log10_CD38"),
  #          n_condition = 1,
  #          resolution = 40,
  #          p_correct = "correlated")
  # )
}
)
