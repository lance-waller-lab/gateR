context("gating")

###################
# gating testthat #
###################

# Faulty datasets
fubar <- randCyto
fubar$g1 <- as.numeric(fubar$g1)

# Tests

test_that("gating throws error with invalid arguments", {

  # Odd-numbered amount of vars
  expect_error(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4"),
           n_condition = 1)
  )

  # A var not available
  expect_error(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "x"),
           n_condition = 1)
  )

  # Non-binary second feature
  expect_error(
    gating(dat = fubar,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "arcsinh_CD3"),
           n_condition = 1)
  )

  expect_error(
    gating(dat = fubar,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "arcsinh_CD3"),
           n_condition = 2)
  )

  # Non-binary third feature
  expect_error(
    gating(dat = fubar1,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "arcsinh_CD3"),
           n_condition = 2)
  )

  # Incorrectly specified alpha
  expect_error(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "arcsinh_CD3"),
           n_condition = 1,
           alpha = 0)
  )

  # Incorrectly specified p_correct
  expect_error(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "arcsinh_CD3"),
           n_condition = 1,
           p_correct = "x")
  )

  # Incorrectly specified n_condition
  expect_error(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "arcsinh_CD3"),
           n_condition = 3)
  )

  # Missing specified n_condition
  expect_error(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "arcsinh_CD3"))
  )
  
  # no results in Gate 1
  expect_error(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38"),
           n_condition = 1,
           p_correct = "uncorrelated Bonferroni")
  )

}
)

test_that("gating works", {

  expect_named(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38"),
           n_condition = 1)
  )

  expect_named(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38"),
           n_condition = 2,
           p_correct = "none")
  )
  
  expect_named(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38"),
           numerator = FALSE,
           n_condition = 1)
  )
  
  expect_named(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38"),
           n_condition = 2,
           resolution = 50,
           p_correct = "correlated Bonferroni")
  )
  
  expect_named(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38"),
           n_condition = 2,
           p_correct = "uncorrelated Bonferroni")
  )
  
  expect_named(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38"),
           n_condition = 2,
           bandw = 1,
           p_correct = "Friston")
  )

  expect_named(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "arcsinh_CD3"),
           n_condition = 1)
  )

  expect_named(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "arcsinh_CD3"),
           n_condition = 2)
  )

  expect_named(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38"),
           n_condition = 1,
           alpha = 0.1)
  )

  expect_named(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38"),
           n_condition = 2,
           plot_gate = TRUE,
           rcols = c("green", "yellow", "purple"))
  )

}
)
