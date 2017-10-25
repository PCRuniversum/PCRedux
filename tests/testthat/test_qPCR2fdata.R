library(PCRedux)

context("qPCR2fdata")

test_that("qPCR2fdata gives the correct dimensions and properties", {
    library(qpcR)
    res_fdata <- qPCR2fdata(testdat)

    expect_that(res_fdata, is_a("fdata"))
    expect_that(length(res_fdata$rangeval) == 2 && res_fdata$rangeval[2] == 49, is_true())
})