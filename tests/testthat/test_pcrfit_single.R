library(PCRedux)

context("pcrfit_single")

test_that("pcrfit_single gives the correct dimensions and properties", {

    library(qpcR)
    res_pcrfit_single <- pcrfit_single(boggy[, 2])

    expect_that(res_pcrfit_single, is_a("data.frame"))
    expect_that(res_pcrfit_single[["changepoint.e.agglo"]] == 3, is_true())
    expect_that(res_pcrfit_single[["changepoint.e.agglo"]] == 3, is_true())
    expect_that(round(res_pcrfit_single[["cpDdiff"]], 2) == 2.54, is_true())
    expect_that(res_pcrfit_single[["top"]] == 5, is_true())
    expect_that(res_pcrfit_single[["bg.start_normalized"]] == 0.05, is_true())
    expect_that(res_pcrfit_single[["bg.stop_normalized"]] == 0.125, is_true())
    expect_that(res_pcrfit_single[["amp.stop_normalized"]] == 0.625, is_true())
    expect_that(res_pcrfit_single[["hookreg_hook"]] == 1, is_true())
    expect_that(res_pcrfit_single[["amptester_shap.noisy"]] == FALSE, is_true())
    expect_that(res_pcrfit_single[["amptester_lrt.test"]] == TRUE, is_true())
    expect_that(res_pcrfit_single[["amptester_rgt.dec"]] == TRUE, is_true())
    expect_that(res_pcrfit_single[["amptester_tht.dec"]] == TRUE, is_true())
    expect_that(res_pcrfit_single[["amptester_slt.dec"]] == TRUE, is_true())
})
