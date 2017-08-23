#' A function to the autocorrelation of amplification curve data 
#' from a quantitative PCR experiment
#' 
#' \code{autocorrelation_test} is a function an autocorrelation analysis
#' from a quantitative PCR experiment.
#' 
#' @param y is the cycle dependent fluorescence amplitude (y-axis).
#' @param n is  
#' @param sig.level is the significance level for the correlation test.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords autocorrelation
#' @examples
#' 
#' # calculate head to tail ratio on noise (negative) amplification curve data
#'
#' autocorrelation_test(y=rnorm(35))
#' 
#' @export autocorrelation_test

autocorrelation_test <- function(y, n=3, sig.level=0.01) {
            cycle_RFU <- try(as.zoo(y), silent=TRUE)
            if(class(cycle_RFU) != "try-error") {
                cycle_RFU_n <- lag(cycle_RFU, k=-n, na.pad=TRUE)
                res_autocorrelation <- cor.test(cycle_RFU[!is.na(cycle_RFU_n)], cycle_RFU_n[!is.na(cycle_RFU_n)])
                if(res_autocorrelation$p.value <= sig.level) {
                    res_autocorrelation <- res_autocorrelation$estimate
                } else {
                    res_autocorrelation <- "n.s."
                }
            } else {
                res_autocorrelation <- NA
            }
            res_autocorrelation
}
