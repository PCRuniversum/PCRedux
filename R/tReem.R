#' A function to Group Amplification Curves According to their Shape
#'
#' \code{tReem} is a function to group amplification curves
#' from a quantitative PCR experiment according to their shape.
#' Either the Pearson correlation coefficient or the Hausdorff
#' distance is used as measure. In most cases the grouping based
#' on the Pearson correlation coefficient is sufficient. The
#' grouping based on the Hausdorff distance can be very slow
#' for large data sets.
#'
#' @param data is the cycle dependent fluorescence amplitude (y-axis).
#' @param cor is a logical parameter. If set true, the Pearson
#' correlation is used as distance measure. If set FALSE the
#' Hausdorff distance will be used.
#' @param k an integer scalar or vector with the desired number of groups.
#' @author Stefan Roediger, Andrej-Nikolai Spiess
#' @keywords autocorrelation
#' @rdname tReem
#' @export tReem
#' @seealso
#'  \code{\link[fda.usc]{metric.hausdorff}}, \code{\link[stats]{cutree}},
#'  \code{\link[PCRedux]{qPCR2fdata}}, \code{\link[stats]{hclust}},
#'  \code{\link[stats]{cor}}
#' @examples
#' # Classify amplification curve data by Hausdorff distance
#' \dontrun{
#' library(qpcR)
#' tReem(testdat[, 1:5])
#' }

tReem <- function(data, cor = TRUE, k = 2) {
  if (is.integer(k) == FALSE) {
    k <- round(k)
  }
  if (k < 2) {
    k <- 2
  }
  y_lim <- range(data[, -1])
  y_quantile <- quantile(unlist(data[, -1]), c(0.25, 0.5, .75))

  if (cor) {
    data_hclust <- hclust(dist(cor(data[, -1], method = "pearson")))
  } else {
    data_hclust <- hclust(dist(metric.hausdorff(qPCR2fdata(data))))
  }
  res_cutree <- cutree(data_hclust, k = k)
  dec <- vector()

  for (i in 1L:k) {
    matplot(
      data[, 1], data[, which(res_cutree == i) + 1],
      type = "l", lty = 1, lwd = 2, main = i, xlab = "Cycle",
      ylab = "RFU", ylim = y_lim
    )
    abline(h = y_quantile, lwd = c(1, 1.5, 1), col = "grey", lty = 2)
    r_tmp <- readline(prompt = "Press [enter] to continue")
    dec <- c(dec, r_tmp)
  }

  df <- data.frame(names(res_cutree), res_cutree, dec = NA)

  for (i in 1L:length(dec)) {
    df[which(res_cutree == i), "dec"] <- dec[i]
  }
  df
}
