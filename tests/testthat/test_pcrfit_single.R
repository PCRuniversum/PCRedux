library(PCRedux)

context("pcrfit_single")

test_that("pcrfit_single gives the correct dimensions and properties", {
  res_pcrfit_single <- pcrfit_single(RAS002[, 2])

  expect_that(res_pcrfit_single, is_a("data.frame"))
  expect_that(round(res_pcrfit_single[["cpDdiff"]], 2) == 2.67, expect_true())
  expect_equal(res_pcrfit_single[["top"]], 25)
  expect_equal(res_pcrfit_single[["bg.stop"]], 15)
  expect_equal(res_pcrfit_single[["amp.stop"]], 40)
  expect_that(res_pcrfit_single[["hookreg_hook"]] == 0, expect_true())
  expect_that(res_pcrfit_single[["amptester_shapiro"]] == FALSE, expect_true())
  expect_that(res_pcrfit_single[["amptester_lrt"]] == TRUE, expect_true())
  expect_that(res_pcrfit_single[["amptester_rgt"]] == TRUE, expect_true())
  expect_that(res_pcrfit_single[["amptester_tht"]] == TRUE, expect_true())
  expect_that(res_pcrfit_single[["amptester_slt"]] == TRUE, expect_true())
  expect_length(res_pcrfit_single, 89)
})
