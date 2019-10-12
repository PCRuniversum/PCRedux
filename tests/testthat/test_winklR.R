library(PCRedux)

context("winklR")

test_that("winklR gives the correct angle", {
  data("RAS002")
  res <- winklR(x = RAS002[, 1], y = RAS002[, 2], normalize = FALSE)

  expect_that(res[[1]], is_a("numeric"))
})
