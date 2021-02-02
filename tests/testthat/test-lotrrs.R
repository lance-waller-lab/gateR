context("lotrrs")

###################
# lotrrs testthat #
###################

# Faulty datasets
fubar <- fubar1 <- randCyto
fubar$g1 <- as.numeric(fubar$g1)
fubar1$g2 <- as.numeric(fubar1$g2)


# Tests

test_that("lotrrs throws error with invalid arguments", {
  
  # Odd-numbered amount of vars
  expect_error(
    lotrrs(dat = randCyto[ , 1:4])
  )
  
  # Non-binary second feature
  expect_error(
    lotrrs(dat = fubar)
  )
  
  # Non-binary third feature
  expect_error(
    lotrrs(dat = fubar1)
  )
  
  # Incorrectly specified alpha
  expect_error(
    lotrrs(dat = randCyto,
        alpha = 0)
  )
  
  # Incorrectly specified p_correct
  expect_error(
    lotrrs(dat = randCyto,
        p_correct = "x")
  )
  
}
)

test_that("lotrrs works", {

  expect_named(
    lotrrs(dat = randCyto)
  )

  expect_named(
    lotrrs(dat = randCyto,
           p_correct = "FDR")
  )
  
  expect_named(
    lotrrs(dat = randCyto,
           resolution = 50,
           p_correct = "correlated Sidak")
  )
  
  expect_named(
    lotrrs(dat = randCyto,
        p_correct = "uncorrelated Sidak")
  )
  
  expect_named(
    lotrrs(dat = randCyto,
           p_correct = "Adler and Hasofer")
  )
  
}
)
