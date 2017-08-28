#' A helper function to convert amplification curve data to the fdata format.
#' 
#' \code{qPCR2fdata} is a helper function 
#' @param data is a data set containing the amplification cycles (1. column) 
#' and the fluorescence (subsequent columns).
#' @author Stefan Roediger, Michal Burdukiewcz
#' @references M. Febrero-Bande, M.O. de la Fuente, others, \emph{Statistical 
#' computing in functional data analysis: The R package fda.usc}, Journal of 
#' Statistical Software. 51 (2012) 1–28. http://www.jstatsoft.org/v51/i04/
#' @keywords fdata
#' @examples
#' 
#' # Calculate slope and intercept on noise (negative) amplification curve data
#' # for the last eight cycles.
#' library(qpcR)
#' library(fda.usc)
#' library(magrittr)
#' 
#' res_fdata <- qPCR2fdata(testdat)
#' res_fdata_colnames <- testdat[-1] %>% colnames()
#' data_colors <- rainbow(length(res_fdata_colnames), alpha=0.5)
#' 
#' par(mfrow=c(1,2))
#' res_fdata %>% plot(., xlab="cycles", ylab="RFU", main="testdat", type="l", 
#'                    lty=1, lwd=2, col=data_colors)
#' legend("topleft", as.character(res_fdata_colnames), pch=19, 
#'          col=data_colors, bty="n", ncol=2)
#' res_fdata_hclust <- metric.hausdorff(res_fdata)
#' plot(hclust(as.dist(res_fdata_hclust))) 
#'  
#' @export qPCR2fdata

qPCR2fdata <- function(data) {
            fdata(t(data[, -1]), argvals=data[, 1], rangeval=range(data[, 1]))
}
