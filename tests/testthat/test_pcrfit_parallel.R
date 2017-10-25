library(PCRedux)

context("pcrfit_parallel")

test_that("pcrfit_parallel gives the correct dimensions and properties", {

    library(qpcR)
    res_pcrfit_parallel <- pcrfit_parallel(boggy[, 1:7])

    expect_that(res_pcrfit_parallel, is_a("data.frame") )
    expect_that(nrow(res_pcrfit_parallel) == 6, is_true())
    expect_that(ncol(res_pcrfit_parallel) == 44, is_true())
})