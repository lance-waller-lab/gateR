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
    lotrrs(dat = randCyto[ , 1:4],
        p_correct = "none")
  )
  
  # Non-binary second feature
  expect_error(
    lotrrs(dat = fubar,
        p_correct = "none")
  )
  
  # Non-binary third feature
  expect_error(
    lotrrs(dat = fubar1,
           p_correct = "none")
  )
  
  # Incorrectly specified alpha
  expect_error(
    lotrrs(dat = randCyto,
        alpha = 0,
        p_correct = "none")
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
    lotrrs(dat = randCyto,
        p_correct = "none")
  )

  expect_named(
    lotrrs(dat = randCyto,
        p_correct = "uncorrelated")
  )

  # expect_named(
  #   lotrrs(dat = randCyto,
  #       resolution = 40,
  #       p_correct = "correlated")
  # )
  # 
  # expect_named(
  #   lotrrs(dat = randCyto,
  #       alpha = 0.1,
  #       p_correct = "none")
  # )
  # 
  # expect_named(
  #   lotrrs(dat = randCyto,
  #       alpha = 0.01,
  #       p_correct = "none")
  # )
  # 
  # expect_named(
  #   lotrrs(dat = randCyto,
  #       doplot = TRUE,
  #       p_correct = "none")
  # )
  # 
  # expect_named(
  #   lotrrs(dat = randCyto,
  #       doplot = TRUE,
  #       rcols = c("green", "yellow", "purple"),
  #       p_correct = "none")
  # )
}
)
