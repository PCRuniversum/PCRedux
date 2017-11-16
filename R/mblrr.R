#' A function to perform a Qunantile-filter based Local Robust Regression
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
#' @param normalize is a logical parameter, which indicates if the amplification curve
#' data should be normalized to the 99 percent quantile of the amplification curve.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords segmented regression
#' @examples
#' 
#' # Perform an mblrr analysis on noise (negative) amplification data of qPCR data
#' # with 35 cycles.
#' library(qpcR)
#' mblrr(x=boggy[, 1], y=boggy[, 2], normalize=TRUE)
#' 
#' @export mblrr
mblrr <-
function(x, y, sig.level=0.01, normalize=FALSE) {

            res_mblrr <- c(NA, NA, NA, NA, NA, NA)

            if(normalize) y <- y/quantile(y, 0.999, na.rm=TRUE)

            res_q25 <- y < quantile(y, 0.25, na.rm=TRUE)
            res_q75 <- y > quantile(y, 0.75, na.rm=TRUE)

            res_q25_lm <- try(suppressWarnings(lmrob(y[res_q25] ~ x[res_q25])), silent=TRUE)
            res_q75_lm <- try(suppressWarnings(lmrob(y[res_q75] ~ x[res_q75])), silent=TRUE)
            
            res_q25_cor.test <- try(cor.test(x[res_q25], y[res_q25]), silent=TRUE)
            res_q75_cor.test <- try(cor.test(x[res_q75], y[res_q75]), silent=TRUE)

            if(class(res_q25_lm) == "try-error" || is.na(res_q25_lm[[1]])[1] == TRUE) {
                                    mblrr_intercept_less <- NA
                                    mblrr_slope_less <- NA
                                    } else {
                                            mblrr_intercept_less <- res_q25_lm$coefficients[[1]]
                                            mblrr_slope_less <- res_q25_lm$coefficients[[2]]
                                    }
            
            if(class(res_q75_lm) == "try-error" || is.na(res_q75_lm[[1]])[1] == TRUE) {
                        mblrr_intercept_more <- NA
                        mblrr_slope_more <- NA
                        } else {
                                mblrr_intercept_more <- res_q75_lm$coefficients[[1]]
                                mblrr_slope_more <- res_q75_lm$coefficients[[2]]
                        }
            
           

            if(class(res_q25_cor.test) == "try-error" || is.na(res_q25_cor.test$p.value)) {
                                            res_q25_cor.test_estimate <- NA
                                    } else {
                                            ifelse(res_q25_cor.test$p.value < sig.level,  
                                                   res_q25_cor.test_estimate <- res_q25_cor.test$estimate, 
                                                   res_q25_cor.test_estimate <- NA)
                                            }

            if(class(res_q75_cor.test) == "try-error" || is.na(res_q75_cor.test$p.value)) {
                                            res_q75_cor.test_estimate <- NA
                                    } else {
                                            ifelse(res_q75_cor.test$p.value < sig.level,  
                                                   res_q75_cor.test_estimate <- res_q75_cor.test$estimate, 
                                                   res_q75_cor.test_estimate <- NA)
                                            }
                                            

            res_mblrr <- c(mblrr_intercept_less,
                           mblrr_slope_less,
                           res_q25_cor.test_estimate,
                           mblrr_intercept_more,
                           mblrr_slope_more,
                           res_q75_cor.test_estimate
                        )
            
            names(res_mblrr) <- c("mblrr_intercept_less",
                                  "mblrr_slope_less",
                                  "mblrr_cor_less",
                                  "mblrr_intercept_more",
                                  "mblrr_slope_more",
                                  "mblrr_cor_more")
            res_mblrr
}
