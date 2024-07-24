library(PCRedux)

context("hookregNL")

test_that("hookregNL gives the correct dimensions and properties", {
  library(qpcR)
  res_hookregNL_positve <- hookregNL(x = boggy[, 1], y = boggy[, 2])
  res_hookregNL_positve_complex <- hookregNL(x = boggy[, 1], y = boggy[, 2], simple = FALSE)
  res_hookregNL_negative <- hookregNL(x = boggy[, 1], y = boggy[, 6])

  expect_that(res_hookregNL_positve, is_a("data.frame"))
  expect_true(round(res_hookregNL_positve[["slope"]], 8) == -0.01596799)
  expect_that(res_hookregNL_positve_complex, is_a("list"))
  expect_that(res_hookregNL_negative, is_a("data.frame"))
  expect_true(res_hookregNL_positve[["hook"]] == 1)
  expect_true(res_hookregNL_positve_complex[["hook"]] == 1)
  expect_true(res_hookregNL_positve_complex[["fit"]][["weights"]][1] == 1)
  expect_true(round(res_hookregNL_positve_complex[["slope"]], 8) == -0.01596799)
  expect_false(res_hookregNL_negative[["hook"]] == 1)
})
