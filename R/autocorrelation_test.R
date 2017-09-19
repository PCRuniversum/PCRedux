#' A function to the autocorrelation of amplification curve data 
#' from a quantitative PCR experiment
#' 
#' \code{autocorrelation_test} is a function an autocorrelation analysis
#' from a quantitative PCR experiment. The result of the funciton is either
#' a correlation coefficient in case the result is significant at a given
#' significance level, or a "n.s." (non-significant) of no correlation could
#' be determined.
#' 
#' @param y is the cycle dependent fluorescence amplitude (y-axis).
#' @param n is  
#' @param sig.level is the significance level for the correlation test., Default: 0.01
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords autocorrelation
#' @rdname autocorrelation_test
#' @export autocorrelation_test
#' @seealso 
#'  \code{\link[zoo]{as.zoo}}, \code{\link[stats]{lag}}, \code{\link[stats]{cor.test}}
#' @examples
#' # Test for autocorrelation in amplification curve data
#' # Load the libraries magrittr for pipes and qpcR for the data
#' library(magrittr)
#' library(qpcR)
#' # Test for autocorrelation in the testdat dataset
#' res_ac <- sapply(2:ncol(testdat), function(i) {
#'                     autocorrelation_test(testdat[, i])
#'                 }
#'          )
#' 
#' # Plot curve data as overview
#' # Define the colors for the amplification curves
#' colors <- rainbow(ncol(testdat)-1, alpha=0.3)
#' # Names of samples
#' samples <- colnames(testdat)[-1]
#v layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE), respect = TRUE)
#' matplot(testdat[, 1], testdat[, -1], xlab="Cycle", ylab="RFU", 
#'         main="testdat dataset", type="l", lty=1, col=colors, lwd=2)
#' legend("topleft", samples, pch=19, col=colors, ncol=2, bty="n")
#'     
#' # Curves rated by a human after analysis of the overview. 1 = positive, 
#' # 0 = negative
#' human_rating <- c(1,1,0,0,1,1,0,0,
#'                   1,1,0,0,1,1,0,0,
#'                   1,1,0,0,1,1,0,0)
#' 
#' # Convert the n.s. (not significant) in 0 and others to 1. 
#' # Combine the results of the atomatic autocorrelation_test as variable "ac", 
#' # the human rated values as variable "hr" in a new data frame (res_ac_hr).
#' res_ac_hr <- data.frame(ac=ifelse(res_ac=="n.s.", 0, 1), 
#'                         hr=human_rating) %>% as.matrix
#' 
#' # Add ratings by humna and autocorrelation_test to the plot
#' par(las=2)
#' barplot(res_ac_hr[, "ac"], main="Rated by autocorrelation_test", col=colors, 
#'         xlab="Sample", ylab="rating", names.arg=samples)
#' barplot(res_ac_hr[, "hr"], main="Rated by a human", col=colors, xlab="Sample", 
#'         ylab="rating", names.arg=samples)

autocorrelation_test <- function (y, n = 3, sig.level = 0.01) 
{
    # Coercing object to class "zoo".
    cycle_RFU <- try(zoo::as.zoo(y), silent = TRUE)
    
    if(class(cycle_RFU) == "zoo") {
        # Compute a lagged version of the cycle, shifting the cycle (time) base 
        # back by a given number of observations
        cycle_RFU_n <- stats::lag(cycle_RFU, k = -n, na.pad = TRUE)
        # Test for correlation between paired samples (cycle & lagged cycle)
        res_autocorrelation <- stats::cor.test(cycle_RFU[!is.na(cycle_RFU_n)], 
            cycle_RFU_n[!is.na(cycle_RFU_n)])
        # Logical analysis of the correlation test and output
        if(res_autocorrelation$p.value <= sig.level) {
            res_autocorrelation <- res_autocorrelation$estimate
        } else{
                res_autocorrelation <- "n.s."
                }
    } else{
           res_autocorrelation <- NA
           }
    res_autocorrelation
}
