#' A function to calculate the slope and intercept of an amplification curve data 
#' from a quantitative PCR experiment.
#' 
#' \code{earlyreg} is a function to calculate the slope and intercept of an 
#' amplification curve data from a quantitative PCR experiment. The number
#' of cycles to be analyzed is defined by the user (default 6 cycles).
#' The output contains the Maximal Information Coefficient (MIC), which 
#' can be interpreted as a correlation measure with a range of [0,1]. 
#' A value of 0 mean statistically independent data and 1 approaches 
#' in "probability for noiseless functional relationships" (see 
#' original study Reshef, D. N. et al. Detecting novel associations in 
#' large data sets. Science, 334, 1518–1524 (2011)).
#' 
#' @param x is the cycle numbers (x-axis).
#' @param y is the cycle dependent fluorescence amplitude (y-axis).
#' @param range is the number of cycles to be used for the regression.
#' @param normalize is a logical parameter which indicates if the amplification curve
#' data should be normalized to the 99 percent percentile of the amplification curve.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords slope intercept
#' @seealso 
#'  \code{\link[robustbase]{lmrob}}
#'  \code{\link[stats]{coefficients}}
#'  \code{\link[minerva]{mine}}
#' @examples
#' 
#' # Calculate slope and intercept on noise (negative) amplification curve data
#' # for the cycles 2 to 7 for the C316.amp data set
#' library(chipPCR)
#' data(C316.amp)
#' 
#' # Plot the data
#' plot(C316.amp[, 2], y=C316.amp[, 3], xlab="Cycle", ylab="RFU", 
#'      main="C316.amp data set", lty=1, type="l")
#' res <- earlyreg(x=C316.amp[, 2], y=C316.amp[, 3], range=6)
#' res
#' @export earlyreg

earlyreg <- function(x, y, range=6, normalize=FALSE) {
    data <- na.omit(cbind(x = x, y = y))
    x <- data[, "x"]
    y <- data[, "y"]
    if (normalize) 
        y <- y/quantile(y, 0.999)
    range <- head(x[-1], range)
    suppressWarnings
    res_lm_fit <- try(suppressWarnings(coefficients(lmrob(y[range] ~ 
        x[range]))), silent = TRUE)
    if (class(res_lm_fit) == "try-error") {
        res_lm_fit <- c(NA, NA)
    }
    
    res_mine <- try(suppressWarnings(mine(x[range], y[range])), silent = TRUE)
    if (class(res_mine) == "try-error") {
       res_mine <- c(NA)
    }
    
    res_lm_fit <- c(res_lm_fit, res_mine[["MIC"]])
    names(res_lm_fit) <- c("intercept", "slope", "earlyreg.MIC")
    res_lm_fit
}
