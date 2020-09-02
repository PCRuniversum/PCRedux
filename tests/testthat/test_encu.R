library(PCRedux)

context("encu")

test_that("encu is a function to calculate numerous features from amplification curve data from a quantitative PCR experiment", {
  library(qpcR)
  res <- encu(testdat[, 1:2])

  expect_is(res$f.tdp, "numeric")
  expect_is(res, "data.frame")
  expect_that(round(res$cp_e.agglo, 4) == 0.0408, is_true())
  expect_that(res$hookreg_hook == 1, is_true())
  expect_that(res$top == 10, is_true())
  expect_length(res, 60)
  expect_true(res$amptester_shapiro == FALSE)
  expect_true(res$amptester_rgt != FALSE)
})

