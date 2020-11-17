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
           n_condition = 1,
           p_correct = "none")
  )

  # A var not available
  expect_error(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "x"),
           n_condition = 1,
           p_correct = "none")
  )

  # Non-binary second feature
  expect_error(
    gating(dat = fubar,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "arcsinh_CD3"),
           n_condition = 1,
           p_correct = "none")
  )

  expect_error(
    gating(dat = fubar,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "arcsinh_CD3"),
           n_condition = 2,
           p_correct = "none")
  )

  # Non-binary third feature
  expect_error(
    gating(dat = fubar1,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "arcsinh_CD3"),
           n_condition = 2,
           p_correct = "none")
  )

  # Incorrectly specified alpha
  expect_error(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "arcsinh_CD3"),
           n_condition = 1,
           alpha = 0,
           p_correct = "none")
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
           n_condition = 3,
           p_correct = "none")
  )

  # Missing specified n_condition
  expect_error(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "arcsinh_CD3"),
           p_correct = "none")
  )

}
)

test_that("gating works", {

  expect_named(
      gating(dat = randCyto,
             vars = c("arcsinh_CD4", "arcsinh_CD38"),
             n_condition = 1,
             p_correct = "none")
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
           n_condition = 1,
           p_correct = "uncorrelated")
  )

  expect_named(
    gating(dat = randCyto,
           vars = c("arcsinh_CD4", "arcsinh_CD38",
                    "arcsinh_CD8", "arcsinh_CD3"),
           n_condition = 1,
           p_correct = "none")
  )

  # expect_named(
  #   gating(dat = randCyto,
  #          vars = c("arcsinh_CD4", "arcsinh_CD38",
  #                   "arcsinh_CD8", "arcsinh_CD3"),
  #          n_condition = 2,
  #          p_correct = "none")
  # )
  # 
  # expect_named(
  #   gating(dat = randCyto,
  #          vars = c("arcsinh_CD4", "arcsinh_CD38"),
  #          n_condition = 1,
  #          alpha = 0.1,
  #          p_correct = "none")
  # )
  # 
  # expect_named(
  #   gating(dat = randCyto,
  #          vars = c("arcsinh_CD4", "arcsinh_CD38"),
  #          n_condition = 1,
  #          alpha = 0.01,
  #          p_correct = "none")
  # )
  # 
  # expect_named(
  #   gating(dat = randCyto,
  #          vars = c("arcsinh_CD4", "arcsinh_CD38"),
  #          n_condition = 1,
  #          numerator = FALSE,
  #          p_correct = "none")
  # )
  # 
  # expect_named(
  #   gating(dat = randCyto,
  #          vars = c("arcsinh_CD4", "arcsinh_CD38"),
  #          n_condition = 2,
  #          doplot = TRUE,
  #          rcols = c("green", "yellow", "purple"),
  #          p_correct = "none")
  # )
  # 
  # expect_named(
  #   gating(dat = randCyto,
  #          vars = c("arcsinh_CD4", "arcsinh_CD38"),
  #          n_condition = 1,
  #          resolution = 40,
  #          p_correct = "correlated")
  # )
}
)
