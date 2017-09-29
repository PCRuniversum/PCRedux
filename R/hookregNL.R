#' @title hookregNL
#' @description \code{hookregNL} is just a nice thing
#' @param x PARAM_DESCRIPTION
#' @param y PARAM_DESCRIPTION
#' @param plot PARAM_DESCRIPTION, Default: FALSE
#' @param level PARAM_DESCRIPTION, Default: 0.99
#' @param simple is a logical parameter. If TRUE (default) only the slope, 
#' confidence interval and decisions are shown as output 
#' @param manualtrim is the number of cycles that should be reomoved from the 
#' background.
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
#' res <- hookregNL(boggy[, 1], boggy[, 2])
#' res
#' 
#' # has no hook
#' res <- hookregNL(boggy[, 1], boggy[, 12])
#' res
#' @seealso 
#'  \code{\link[qpcR]{pcrfit}}
#'  \code{\link[stats]{confint}}
#' @rdname hookregNL
#' @export hookregNL

hookregNL <- function(x, y, plot=FALSE, level=0.99, simple=TRUE, manualtrim=5) {
  # Create data, remove missing values (manualtrim) and remove first 5 cycles to 
  # avoid fitting baseline slopes.
  
  data <- na.omit(data.frame(cycles = x, fluo = y))
  data <- data[-(1:manualtrim), ]
  
  # fit a 6-parameter log-logistic model
  fit <- try(pcrfit(data, 1, 2, l6), silent = TRUE)
  l6 <- NULL
  if(inherits(fit, "try-error")) {message("fitting failed.")}
  
  if(plot && !inherits(fit, "try-error")) plot(fit)
  
  # Confidence interval for slope parameter 'k'
  if(inherits(fit, "try-error")) {
        slope <- NA
        confslope <- c(NA,NA)
        message("Could not calculate confidence interval.")
        } else {
            slope <- coefficients(fit)[6]
            confslope <- try(stats::confint(fit, level = level)[6, ], silent = TRUE)
            if (inherits(confslope, "try-error")) {
                confslope <- c(NA,NA)
            }
        }
  
  # Decision
  hook <- ifelse(!is.na(confslope[1]) && 
                       confslope[1] < 0 && 
                       confslope[2] < 0  && 
                       class(fit) != "try-error", 
                       1, 0)
  
  confslope_simple <- if(!is.na(confslope)) {
                            data.frame(CI.low=confslope[[1]], CI.up=confslope[[2]])  
                      } else{
                            data.frame(CI.low=NA, CI.up=NA)
                      }
  
  # Output
  if(simple){
        res <- data.frame(slope=slope[[1]], CI.low=confslope_simple[[1]], CI.up=confslope_simple[[2]], hook=hook)
  } else {
        res <- list(fit=fit, slope=slope, conf=confslope, hook=hook)
        }
   res
}
