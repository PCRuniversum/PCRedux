library(PCRedux)

context("autocorrelation_test")

test_that("autocorrelation_test gives the correct dimensions and properties", {
  library(qpcR)
  res_ac_positive <- autocorrelation_test(testdat[, 2])
 
  expect_that(res_ac_positive, is_a("numeric"))
})
