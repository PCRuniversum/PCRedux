library(PCRedux)

context("autocorrelation_test")

test_that("autocorrelation_test gives the correct dimensions and properties", {
  library(qpcR)
  res_ac_positive <- autocorrelation_test(testdat[, 2])
  res_ac_negative <- autocorrelation_test(testdat[, 4])

  expect_that(res_ac_positive, is_a("numeric"))
  expect_that(res_ac_negative, is_a("numeric"))
#   expect_that(res_ac_negative, is_a("character"))
  expect_that(signif(res_ac_positive[[1]]), equals(0.581476))
#   expect_that(res_ac_negative[[1]], equals("n.s."))
  expect_that(res_ac_negative[[1]], equals(0))
})
