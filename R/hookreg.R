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
#' @param normalize is a logical parameter indicating if the data should be 
#' normalized to the 0.999 quantile
#' @param sig.level defines the significance level to test for a significant 
#' regression
#' @param CI.level confidence level required for the slope
#' @param robust is a logical parameter indicating if the data should be 
#' analyzed be a robust linear regression (\code{lmrob}).
#' @author Stefan Roediger, Michal Burdukiewcz
#' @references K. Barratt, J.F. Mackay, \emph{Improving Real-Time PCR Genotyping 
#' Assays by Asymmetric Amplification}, J. Clin. Microbiol. 40 (2002) 1571--1572. 
#' doi:10.1128/JCM.40.4.1571-1572.2002.
#' @keywords slope intercept hook
#' @examples
#' 
#' # Calculate slope and intercept on noise (negative) amplification curve data
#' # for the last eight cycles.
#' 
#' library(qpcR)
#' library(magrittr)
#' 
#'  res_hook <- sapply(2:ncol(boggy), function(i) {
#'      hookreg(x=boggy[, 1], y=boggy[, i])}) %>% t %>% 
#'      data.frame(sample=colnames(boggy)[-1],.)
#' res_hook
#' 
#' data_colors <- rainbow(ncol(boggy[, -1]), alpha=0.5)
#' cl <- kmeans(na.omit(res_hook[, 2:3]), 2)$cluster
#'
#' par(mfrow=c(1,2))
#' matplot(x=boggy[, 1], y=boggy[, -1], xlab="Cycle", ylab="RFU", 
#'  main="boggy Data Set", type="l", lty=1, lwd=2, col=data_colors)
#'  legend("topleft", as.character(res_hook$sample), pch=19, 
#'          col=data_colors, bty="n")
#'
#' plot(res_hook$intercept, res_hook$slope, pch=19, cex=2, col=data_colors,
#'  xlab="intercept", ylab="Slope", 
#'  main="Clusters of Amplification Curves with an Hook Effect-like Curvature\nboggy Data Set")
#'  points(res_hook$intercept, res_hook$slope, col=cl, pch=cl, cex=cl)
#'  legend("topright", c("Strong Hook effect", " Weak Hook effect"), pch=c(1,2), col=c(1,2), bty="n")
#'  text(res_hook$intercept, res_hook$slope, res_hook$sample)
#'  
#' @export hookreg

hookreg <- function(x, y, normalize=TRUE, sig.level=0.005, CI.level=0.995, robust=FALSE) {
    # Remove missing values from data
    data <- na.omit(data.frame(x=x, y=y))
    x <- data[, "x"]
    y <- data[, "y"]
    
    # Quantile Normaliation of amplification data
    if(normalize) {y <- y / quantile(y, 0.999)}
    
    # narrow the range of the potential hook region by a 75% quantile 
    # filter
    hook_quantile_range <- which(y >= quantile(y, c(0.75)))[1]:length(y)

    # narrow the range of the potential hook region by a 75% quantile 
    # filter
    hook_max_range <- which(y == max(y[hook_quantile_range]))[1]

    # Determine putative hook range
    range <- hook_max_range:length(x)
    
    hook_delta <- length(hook_max_range:length(x))

    if(hook_max_range < length(x) && hook_delta >= 5) {
            # Regression for putative hook range
            if(robust) {
                        res_lm_fit <- try(lmrob(y[range] ~ x[range]), silent=TRUE)
                        } else {
                        res_lm_fit <- try(lm(y[range] ~ x[range]), silent=TRUE)
                        }

            # Statistics for regression
            res_lm_fit_summary <- try(summary(res_lm_fit))$coefficients[2, 4]
            res_lm_fit_coefficients <- coefficients(res_lm_fit)
            res_lm_fit_confint <- confint(res_lm_fit, level=CI.level)
            res_hook_significance <- ifelse(res_lm_fit_summary < sig.level, TRUE, FALSE)
            res_lm_fit_confint_decision <- ifelse(res_lm_fit_confint[2, 1] < 0 && 
                                                   res_lm_fit_confint[2, 2] < 0, 
                                                   TRUE, FALSE)
            dec_hook <- ifelse(res_hook_significance == TRUE || res_lm_fit_confint_decision == TRUE, TRUE, FALSE)
            
            res_hookreg <- c(res_lm_fit_coefficients[[1]],
                            res_lm_fit_coefficients[[2]],
                            hook_max_range,
                            hook_delta,
                            res_lm_fit_summary,
                            res_lm_fit_confint[1, 2],
                            res_lm_fit_confint[2, 2],
                            res_hook_significance,
                            res_lm_fit_confint_decision,
                            dec_hook)
            
    } else {
            res_hookreg <- c(NA, NA, NA, NA, NA, NA, NA, FALSE, FALSE, FALSE)
            }
    names(res_hookreg) <- c("intercept", "slope", "hook.start", "hook.delta", "p.value", 
                            "CI.low", "CI.up", "hook.fit", "hook.CI", "hook")
    res_hookreg
}
