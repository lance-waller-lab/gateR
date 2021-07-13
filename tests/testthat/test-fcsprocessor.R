context("fcsprocessor")

#########################
# fcsprocessor testthat #
#########################

# Faulty datasets
# fubar <- randCyto

# Tests

test_that("fcsprocessor throws error with invalid arguments", {

  # LABEL DESCRIBING ERROR
  expect_error(
    fcsprocessor(...) # creates error for now, develop more informative test(s)
  )

}
)

test_that("fcsprocessor works", {

  expect_named(
    #fcsprocessor(...) # uncomment and create test
    print(list("x" = "THIS IS A PLACEHOLDER", "y" = "DELETE THIS LINE")) # Delete after developing more informative test(s) 
  )

}
)
