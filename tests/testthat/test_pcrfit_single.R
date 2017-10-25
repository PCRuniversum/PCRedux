library(PCRedux)

context("pcrfit_single")

test_that("pcrfit_single gives the correct dimensions and properties", {

    library(qpcR)
    res_pcrfit_single <- pcrfit_single(boggy[, 2])

    expect_that(res_pcrfit_single, is_a("data.frame"))
    expect_that(res_pcrfit_single[["changepoint.e.agglo"]] == 3, is_true())
})