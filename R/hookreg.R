#' A function to calculate the slope and intercept of an amplification curve data 
#' from a quantitative PCR experiment at the end of the data stream.
#' 
#' \code{hookreg} is a function to calculate the slope and intercept of an 
#' amplification curve data from a quantitative PCR experiment. The idea is that
#' a strong negative slope at the end of an amplification curve is indicative for
#' a hook effect (see Barratt and Mackay 2002).
#'
#' @param x is the cycle numbers (x-axis).
#' @param y is the cycle dependent fluorescence amplitude (y-axis).
#' @author Stefan Roediger, Michal Burdukiewcz
#' @references K. Barratt, J.F. Mackay, \emph{Improving Real-Time PCR Genotyping 
#' Assays by Asymmetric Amplification}, J. Clin. Microbiol. 40 (2002) 1571–1572. 
#' doi:10.1128/JCM.40.4.1571-1572.2002.
#' @keywords slope intercept hook
#' @examples
#' 
#' # Calculate slope and intercept on noise (negative) amplification curve data
#' # for the last eight cycles.
#'
#' hookreg(x=1:35, y=rnorm(35))
#' 
#' @export hookreg

hookreg <- function(x, y) {
    # narrow the range of the potential hook region by a 75% quantile 
    # filter
    hook_quantile_range <- which(y >= quantile(y, c(0.75)))[1]:length(y)
    
    # narrow the range of the potential hook region by a 75% quantile 
    # filter
    hook_max_range <- which(y == max(y[hook_quantile_range]))
    
    if(hook_max_range < length(x) && length(hook_max_range:length(x) >= 4)) {
    # Determine putative hook range
    range <- hook_max_range:length(x)         
    # Regression for putative hook range   
    res_lm_fit <- try(coefficients(lmrob(y[range] ~ x[range])), silent=TRUE)
    
    if(class(res_lm_fit) == "try-error") {res_lm_fit <- c(0, 0)}
    names(res_lm_fit) <- c("intercept", "slope")
    res_lm_fit
    } else {
        res_lm_fit <- c(0, 0)
        names(res_lm_fit) <- c("intercept", "slope")
        res_lm_fit
    }
}
