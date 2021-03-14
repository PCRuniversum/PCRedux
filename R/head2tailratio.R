#' A function to calculate to head to tail ratio of amplification curve data
#' from a quantitative PCR experiment
#'
#' \code{head2tailratio} is a function to calculate the ratio of the head and the
#' tail of a quantitative PCR amplification curve. In this test, only the head
#' (first six cycles) and the tail (last six cycles) form the region of interest
#' (ROI).
#' @return gives a \code{numeric} (S3 class, type of \code{double}) as output for the head to tail ratio
#'
#' @param y is the cycle dependent fluorescence amplitude (y-axis).
#' @param normalize is a logical parameter, which indicates if the amplification curve.
#' @param slope_normalizer is a logical parameter, which indicates if the
#' head2tailratio should be normalized to the slope of the ROI.
#' @param verbose is a logical parameter, which indicates if all the values,
#' parameters and coefficients of the analysis should be shown.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords ratio head tail
#' @examples
#'
#' # calculate head to tail ratio on amplification curve data
#'
#' library(qpcR)
#'
#' res_head2tailratio <- sapply(2:ncol(competimer), function(i) {
#'    head2tailratio(y=competimer[, i], normalize=TRUE, slope_normalizer=TRUE)
#' })
#'
#' res_head2tailratio_cluster <- kmeans(res_head2tailratio, 3)$cluster
#'
#' matplot(competimer[, 1], competimer[, -1], xlab="Cycle", ylab="RFU",
#'         main="competimer data set", type="l", lty=1, col=res_head2tailratio_cluster, lwd=2)
#'
#' @export head2tailratio
head2tailratio <- function(y, normalize=FALSE, slope_normalizer=FALSE, verbose=FALSE) {

  # Remove missing values
  y <- na.omit(y)
  # Normalize data if needed
  if (normalize) y <- y / quantile(y, 0.999)
  # Create denovo the abscissa values
  y_length <- length(y)


  # Determine the head and tail values
  y_head <- head(y, 6L)
  y_tail <- tail(y, 6L)
  y_reg <- c(y_head, y_tail)

  x_reg <- c(1:6, (y_length - 5L):y_length)

  # Calculate the ratio of the head and the tail
  res_head_tail_ratio <- try(median(y_head) / median(y_tail))

  # Regression analysis of the ROI

  head_tail_fit <- try(lm(y_reg ~ x_reg))

  # Normalize the the head/tail ratio to the slope of the data set,
  # provided that the slope is significant.
  if (slope_normalizer) {
    res_lm_fit_summary <- try(summary(head_tail_fit))$coefficients[2, 4]

    res_hook_significance <- ifelse(res_lm_fit_summary < 0.01, TRUE, FALSE)

    head_tail_slope <- try(coefficients(lm(y_reg ~ x_reg))[[2]])

    if (!inherits_error(head_tail_slope) || res_hook_significance == TRUE) {
      res_head_tail_ratio <- res_head_tail_ratio / head_tail_slope
    }
  }

  res_head_tail_ratio_verbose <- list(x_roi = x_reg, y_roi = y_reg, fit_lm = head_tail_fit, head_tail_ratio = res_head_tail_ratio)

  if (!inherits(res_head_tail_ratio, "numeric")) {
    res_head_tail_ratio <- NA
  }
  if (verbose) {
    res_head_tail_ratio_verbose
  } else {
    res_head_tail_ratio
  }
}
