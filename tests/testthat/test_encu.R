library(PCRedux)

context("encu")

test_that("encu is a function to calculate numerous features from amplification curve data from a quantitative PCR experiment", {
  library(qpcR)
  res <- encu(RAS002[, 1:2])
  
  expect_is(res$f.tdp, "numeric")
  expect_is(res, "data.frame")
  expect_true(res$hookreg_hook == 0)
  expect_true(res$top == 25)
  expect_length(res, 92)
  expect_true(res$amptester_shapiro == FALSE)
  expect_true(res$amptester_rgt == TRUE)
})

