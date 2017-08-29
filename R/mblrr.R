#' A function perform the Median based Local Robust Regression (mblrr) 
#' from a quantitative PCR experiment
#' 
#' \code{mblrr} is a function to perform the Median based Local Robust Regression (mblrr) 
#' from a quantitative PCR experiment. In detail, this function attempts to break the 
#' amplification curve in two parts (head (~background) and tail (~plateau)). 
#' Subsequent, a robust linear regression analysis (\code{\link{lmrob}}) is 
#' preformed individually on both parts. The rational behind this analysis is 
#' that the slope and intercept of an amplification curve differ in the 
#' background and plateau region.
#' @details 
#' \emph{mblrr_intercept_less} is the intercept of the head region,
#' \emph{mblrr_slope_less} is the slope of the head region,
#' \emph{mblrr_cor_less} is the coefficient of correlation of the head region,
#' \emph{mblrr_intercept_more} is the intercept of the tail region,
#' \emph{mblrr_intercept_more} is the slope of the tail region,
#' \emph{mblrr_cor_more} is the coefficient of correlation of the tail region
#' 
#' @param x is the cycle numbers (x-axis).
#' @param y is the cycle dependent fluorescence amplitude (y-axis).
#' @param sig.level is the significance level for the correlation test.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords segmented regression
#' @examples
#' 
#' # Perform an mblrr analysis on noise (negative) amplification data of qPCR data
#' # with 35 cycles.
#'
#' mblrr(x=1:35, y=rnorm(35))
#' 
#' @export mblrr
mblrr <-
function(x, y, sig.level=0.01) {
            res_less_than_median <- y < median(y)
            res_more_than_median <- y > median(y)
            
            if(class(res_less_than_median) == "logical" && class(res_less_than_median) == "logical") {
                res_less_than_median_lm <- try(lmrob(y[res_less_than_median] ~ x[res_less_than_median]), silent=TRUE)
                res_more_than_median_lm <- try(lmrob(y[res_more_than_median] ~ x[res_more_than_median]), silent=TRUE)
                
                if(class(res_less_than_median_lm) != "try-error" & class(res_more_than_median_lm) != "try-error"){
                    res_less_than_median_cor.test <- try(cor.test(x[res_less_than_median], y[res_less_than_median]), silent=TRUE)
                    res_more_than_median_cor.test <- try(cor.test(x[res_more_than_median], y[res_more_than_median]), silent=TRUE)
                    
                    if(class(res_less_than_median_cor.test) != "try-error" & class(res_more_than_median_cor.test) != "try-error"){
                        ifelse(res_less_than_median_cor.test$p.value <= sig.level, 
                                        res_less_than_median_cor.test_estimate <- res_less_than_median_cor.test$estimate, 
                                        res_less_than_median_cor.test_estimate <- NA
                              )
                        ifelse(res_more_than_median_cor.test$p.value <= sig.level, 
                                        res_more_than_median_cor.test_estimate <- res_more_than_median_cor.test$estimate, 
                                        res_more_than_median_cor.test_estimate <- NA
                              )
                    }
                res_mblrr <- c(res_less_than_median_lm$coefficients[1],
                               res_less_than_median_lm$coefficients[2],
                               res_less_than_median_cor.test_estimate,
                               res_more_than_median_lm$coefficients[1],
                               res_more_than_median_lm$coefficients[2],
                               res_more_than_median_cor.test_estimate
                               )
                                
                }
            } else {
                res_mblrr <- c(NA, NA, NA, NA, NA, NA)
            }
            names(res_mblrr) <- c("mblrr_intercept_less",
                                      "mblrr_slope_less",
                                      "mblrr_cor_less",
                                      "mblrr_intercept_more",
                                      "mblrr_slope_more",
                                      "mblrr_cor_more")
            res_mblrr                          
}
