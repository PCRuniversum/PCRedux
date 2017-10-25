library(PCRedux)

context("visdat_pcrfit_parallel")

test_that("visdat_pcrfit_parallel uses all types", {

    library(qpcR)
    res_pcrfit_parallel <- pcrfit_parallel(boggy[, 1:7])

    plot_visdat_pcrfit_parallel <- visdat_pcrfit_parallel(res_pcrfit_parallel, type="all", interactive=FALSE)

    expect_that(plot_visdat_pcrfit_parallel, is_a("gg"))
    expect_that(nrow(plot_visdat_pcrfit_parallel[["data"]]), equals(264))
})

test_that("visdat_pcrfit_parallel uses the qpcR package types only", {

    library(qpcR)
    res_pcrfit_parallel <- pcrfit_parallel(boggy[, 1:7])

    plot_visdat_pcrfit_parallel <- visdat_pcrfit_parallel(res_pcrfit_parallel, type="qpcR", interactive=FALSE)

    expect_that(plot_visdat_pcrfit_parallel, is_a("gg"))
    expect_that(nrow(plot_visdat_pcrfit_parallel[["data"]]), equals(78))
})