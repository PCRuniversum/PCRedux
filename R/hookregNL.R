#' @title hookregNL - A function to calculate the slope of amplification curves in the tail region
#' @description \code{hookregNL} is a function to calculate the slope and intercept of an
#' amplification curve from a quantitative PCR experiment. The idea is that a
#' strong negative slope at the end of an amplification curve is indicative for a
#' hook effect (see Barratt and Mackay 2002). In contrast to
#' \code{\link[PCRedux]{hookreg}} fits this function a sex-parameter model to the
#' amplification curve and extracts the coefficient, which describes the slope.
#' @return gives a \code{numeric} (S3 class, type of \code{double}) as output 
#' for the detection of a hook
#' @param x is the cycle numbers (x-axis).
#' @param y is the cycle dependent fluorescence amplitude (y-axis).
#' @param plot is a logical parameter indicating if the data should be plotted, Default: FALSE.
#' @param level the confidence level required, Default: 0.99.
#' @param simple is a logical parameter. If TRUE (default) only the slope,
#' confidence interval and decisions are shown as output
#' @param manualtrim is the number of cycles that should be removed from the
#' background.
#' (\code{\link[base]{data.frame}}). If FALSE, a \code{\link[base]{list}}
#' including the 6-parameter model is the output.
#' @author Andrej-Nikolai Spiess, Stefan Roediger, Michal Burdukiewcz
#' @references K. Barratt, J.F. Mackay, \emph{Improving Real-Time PCR Genotyping
#' Assays by Asymmetric Amplification}, J. Clin. Microbiol. 40 (2002) 1571--1572.
#' doi:10.1128/JCM.40.4.1571-1572.2002.
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

hookregNL <- function(x, y, plot=FALSE, level=0.995, simple=TRUE, manualtrim=5) {
  # Create data, remove missing values (manualtrim) and remove first 5 cycles to
  # avoid fitting baseline slopes.
  
  data <- na.omit(data.frame(cycles = x, fluo = y))
  data <- data[-(1:manualtrim), ]
  
  # fit a 6-parameter log-logistic model
  fit <- try(pcrfit(data, 1, 2, l6), silent = TRUE)
  l6 <- NULL
  if (inherits_error(fit)) {
    message("fitting failed.")
  }
  
  if (plot && !inherits_error(fit)) plot(fit)
  
  # Confidence interval for slope parameter 'k'
  if (inherits_error(fit)) {
    slope <- NA
    confslope <- c(NA, NA)
    message("Could not calculate confidence interval.")
  } else {
    slope <- coefficients(fit)[6]
    confslope <- try(stats::confint(fit, level = level)[6, ], silent = TRUE)
    if (inherits_error(confslope)) {
      confslope <- c(NA, NA)
    }
  }
  
  # Decision
  hook <- if(!is.na(confslope[1]) &&
             confslope[1] < 0 &&
             confslope[2] < 0 &&
             !inherits_error(fit)) {
    1
  }  else {
    0 
  }
  
  confslope_simple <- if (!is.na(confslope[1])) {
    data.frame(CI.low = confslope[[1]], CI.up = confslope[[2]])
  } else {
    data.frame(CI.low = NA, CI.up = NA)
  }
  
  # Output
  if (simple) {
    res <- data.frame(slope = slope[[1]], CI.low = confslope_simple[[1]], CI.up = confslope_simple[[2]], hook = hook)
  } else {
    res <- list(fit = fit, slope = slope, conf = confslope, hook = hook)
  }
  res
}
