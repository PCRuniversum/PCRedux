#' A function to calculate the slope and intercept of an amplific curve data 
#' from a quantitative PCR experiment.
#' 
#' \code{earlyreg} is a function to calculate the slope and intercept of an 
#' amplific curve data from a quantitative PCR experiment.
#' 
#' @param y is the cycle dependent fluorescence amplitude (y-axis).
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords slope intercept
#' @examples
#' 
#' # Calculate slope and intercept on noise (negative) amplification curve data
#' # for the cycles 2 to 10.
#'
#' earlyreg(1:35, y=rnorm(35), range=3:11)
#' 
#' @export earlyreg

earlyreg <- function(x, y, range=10) {
            
            range <- head(x[-1], range)
            
            res_lm_fit <- try(coefficients(lm(y[range] ~ x[range])), silent=TRUE)
            
            if(class(res_lm_fit) == "try-error") {
                res_lm_fit <- c(NA, NA)
                }
            names(res_lm_fit) <- c("intercept", "slope")
            res_lm_fit
}
