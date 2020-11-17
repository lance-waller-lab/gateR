context("rrs")

################
# rrs testthat #
################

# Faulty datasets
fubar <- randCyto
fubar$g1 <- as.numeric(fubar$g1)

# Tests

test_that("rrs throws error with invalid arguments", {

  # Odd-numbered amount of vars
  expect_error(
    rrs(dat = randCyto[, 1:4],
           p_correct = "none")
  )

  # Non-binary second feature
  expect_error(
    rrs(dat = fubar,
           p_correct = "none")
  )

  # Incorrectly specified alpha
  expect_error(
    rrs(dat = randCyto,
           alpha = 0,
           p_correct = "none")
  )

  # Incorrectly specified p_correct
  expect_error(
    rrs(dat = randCyto,
           p_correct = "x")
  )

}
)

test_that("rrs works", {

  expect_named(
    rrs(dat = randCyto,
        p_correct = "none")
  )

  expect_named(
    rrs(dat = randCyto,
        p_correct = "uncorrelated")
  )

  # expect_named(
  #   rrs(dat = randCyto,
  #       resolution = 40,
  #       p_correct = "correlated")
  # )
  # 
  # expect_named(
  #   rrs(dat = randCyto,
  #       alpha = 0.1,
  #       p_correct = "none")
  # )
  # 
  # expect_named(
  #   rrs(dat = randCyto,
  #       alpha = 0.01,
  #       p_correct = "none")
  # )
  # 
  # expect_named(
  #   rrs(dat = randCyto,
  #       doplot = TRUE,
  #       p_correct = "none")
  # )
  # 
  # expect_named(
  #   rrs(dat = randCyto,
  #       doplot = TRUE,
  #       rcols = c("green", "yellow", "purple"),
  #       p_correct = "none")
  # )
}
)
