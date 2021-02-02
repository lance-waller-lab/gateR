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
    rrs(dat = randCyto[, 1:4])
  )

  # Non-binary second feature
  expect_error(
    rrs(dat = fubar)
  )

  # Incorrectly specified alpha
  expect_error(
    rrs(dat = randCyto,
           alpha = 0)
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
    rrs(dat = randCyto)
  )
  
  expect_named(
    rrs(dat = randCyto,
        p_correct = "FDR")
  )

}
)
