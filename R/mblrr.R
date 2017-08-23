#' A function perform the Median based Local Robust Regression (mblrr) 
#' from a quantitative PCR experiment
#' 
#' \code{mblrr} is a function to perform the Median based Local Robust Regression (mblrr) 
#' from a quantitative PCR experiment.
#' 
#' @param x is the cycles.
#' @param y is the cycle dependent fluorescence amplitude (y-axis).
#' @param sig.level is the significance level for the correlation test.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords segmented regression
#' @examples
#' 
#' # Perfrom an mblrr analysis on noise (negative) amplification data of
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
                res_less_than_median_lm <- try(robustbase::lmrob(y[res_less_than_median] ~ x[res_less_than_median]), silent=TRUE)
                res_more_than_median_lm <- try(robustbase::lmrob(y[res_more_than_median] ~ x[res_more_than_median]), silent=TRUE)
                
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
