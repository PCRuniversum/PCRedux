#' A function to calculate the slope and intercept of an amplification curve data
#' from a quantitative PCR experiment.
#'
#' \code{earlyreg} is a function to calculate the slope and intercept of an
#' amplification curve data from a quantitative PCR experiment. The number
#' of cycles to be analyzed is defined by the user (default 6 cycles).
#' The output contains the Maximal Information Coefficient (MIC), which
#' can be interpreted as a correlation measure with a range of [0,1].
#' A value of 0 mean statistically independent data and 1 approaches
#' in "probability for noiseless functional relationships" (see
#' original study by Reshef, D. N. et al. Detecting novel associations in
#' large data sets. Science, 334, 1518-1524 (2011)).
#'
#' @param x is the cycle numbers (x-axis).
#' @param y is the cycle dependent fluorescence amplitude (y-axis).
#' @param range is the number of cycles to be used for the regression.
#' @param normalize is a logical parameter which indicates if the amplification curve
#' data should be normalized to the 99 percent percentile of the amplification curve.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords slope intercept
#' @seealso
#'  \code{\link[robustbase]{lmrob}}
#'  \code{\link[stats]{coefficients}}
#' @examples
#'
#' # Calculate slope and intercept on noise (negative) amplification curve data
#' # for the cycles 2 to 7 for the C316.amp data set
#' library(chipPCR)
#' data(C316.amp)
#'
#' # Plot the data
#' plot(C316.amp[, 2], y=C316.amp[, 3], xlab="Cycle", ylab="RFU",
#'      main="C316.amp data set", lty=1, type="l")
#' res <- earlyreg(x=C316.amp[, 2], y=C316.amp[, 3], range=5)
#' res
#' @export earlyreg

earlyreg <- function(x, y, range = 5, normalize = FALSE) {
  data <- na.omit(cbind(x = x, y = y))

  x <- data[, "x"]
  if (is.integer(x) == FALSE) {
    x <- 1L:length(x)
  }
  if (range < 2) {
    range <- 2
  }
  if (range > length(x)) {
    range <- length(x)
  }


  y <- data[, "y"]

  if (normalize) {
    y <- y / quantile(y, 0.99)
  }

  range_ht <- head(x, range)

#   model <- try(suppressWarnings(lmrob(y[range_ht] ~ x[range_ht])), silent = TRUE)
  model <- try(suppressWarnings(lm(y[range_ht] ~ x[range_ht])), silent = TRUE)
  
  res_lm_fit <- try(suppressWarnings(
    coefficients(model)
  ), silent = TRUE)
  
  res_lm_sigma <- try(suppressWarnings(summary(model)$sigma), silent = TRUE)

  if (inherits_error(res_lm_fit)) {
    res_lm_fit <- c(0, 0)
  }
  
  if (inherits_error(res_lm_sigma)) {
    res_lm_sigma <- c(0)
  }

  output <- c(res_lm_fit, res_lm_sigma)
  
  names(output) <- c("intercept", "slope", "sigma")
  output
}
