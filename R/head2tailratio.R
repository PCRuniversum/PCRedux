#' A function to calculate to head to tail ratio of amplification curve data 
#' from a quantitative PCR experiment
#' 
#' \code{head2tailratio} is a function to calculate the ratio of the head and the
#' tail of a quantitative PCR amplification curve.
#' 
#' @param y is the cycle dependent fluorescence amplitude (y-axis).
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords ratio head tail
#' @examples
#' 
#' # calculate head to tail ratio on noise (negative) amplification curve data
#'
#' head2tailratio(y=rnorm(35))
#' 
#' @export head2tailratio
head2tailratio <- function(y) {
    res_head_tail_ratio <- try(
        median(head(y))/median(tail(y))
    )
    if(class(res_head_tail_ratio) != "numeric") {
        res_head_tail_ratio <- NA
        }
    res_head_tail_ratio
}
