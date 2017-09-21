#' @title hookregNL
#' @description \code{hookregNL} is just a nice thing
#' @param x PARAM_DESCRIPTION
#' @param y PARAM_DESCRIPTION
#' @param plot PARAM_DESCRIPTION, Default: FALSE
#' @param level PARAM_DESCRIPTION, Default: 0.99
#' @param simple is a logical parameter. If TRUE (default) only the slope, 
#' confidence interval and decisions are shown as output 
#' (\code{\link[base]{data.frame}}). If FALSE, a \code{\link[base]{list}} 
#' including the 6-parameter model is the output.
#' @author Andrej-Nikolai Spiess, Stefan Roediger, Michal Burdukiewcz
#' @return OUTPUT_DESCRIPTION 
#' @details DETAILS
#' @examples 
#' # Analyze data from the boggy data set for potential hook effect like 
#' # curvature
#' library(qpcR)
#' # has hook
#' res <- hookregNL(boggy[, 1], boggy[, 12])
#' res
#' 
#' # has no hook
#' res <- hookregNL(reps[, 1], reps[, 2])
#' res
#' @seealso 
#'  \code{\link[qpcR]{pcrfit}}
#'  \code{\link[stats]{confint}}
#' @rdname hookregNL
#' @export hookregNL

hookregNL <- function(x, y, plot=FALSE, level=0.99, simple=TRUE) {
  # Create data and remove first 5 cycles to 
  # Avoid fitting baseline slopes
  data <- cbind(cycles = x, fluo = y)
  data <- data[-(1:1), ]
  
  # fit a 6-parameter log-logistic model
  fit <- try(qpcR::pcrfit(data, 1, 2, model=l6), silent = TRUE)
  l6 <- NULL
  if (inherits(fit, "try-error")) {
    print("fitting failed.")
    return(NA)
  }
  if (plot) plot(fit)
  
  # Confidence interval for slope parameter 'k'
  slope <- coefficients(fit)[6]
  confslope <- try(stats::confint(fit, level = level)[6, ], silent = TRUE)
  if (inherits(confslope, "try-error")) {
    print("Could not calculate confidence interval.")
    confslope <- NA
  }
  
  # Decision
  hook <- ifelse(!is.na(confslope) && confslope[1] < 0 && confslope[2] < 0, TRUE, FALSE)
  
  confslope_simple <- if(!is.na(confslope)) {
                            data.frame(CI.low=confslope[1], CI.up=confslope[2])  
                      } else{
                            data.frame(CI.low=NA, CI.up=NA)
                      }
  
  # Output
  if(simple){
        return(data.frame(slope=slope, conf=confslope_simple, hook=hook))
  } else {
        return(list(fit=fit, slope=slope, conf=confslope, hook=hook))
        }
}
