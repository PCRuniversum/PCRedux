#' PCRedux - quantitative PCR machine learning helper tool
#'
#' @description \code{PCRedux} package is a toolbox for the analysis of sigmoid curve (qPCR) data.
#' 
#' @section ML:
#' T.B.D.
#' 
#' @importFrom bcp bcp
#' @importFrom changepoint cpt.meanvar
#' @importFrom chipPCR amptester bg.max CPP smoother
#' @importFrom ecp e.agglo
#' @importFrom fda.usc fdata
#' @importFrom graphics matplot par
#' @importFrom grDevices rainbow
#' @importFrom stats coefficients confint cor.test lag lm median quantile
#' @importFrom utils head tail
#' @importFrom magrittr %>%
#' @importFrom MBmca diffQ diffQ2 mcaPeaks
#' @importFrom pracma polyarea
#' @importFrom qpcR efficiency LRE mselect pcrfit sliwin takeoff
#' @importFrom robustbase lmrob
#' @importFrom zoo as.zoo
#' @author Stefan Roediger, Michal Burdukiewcz, Andrej-Nikolai Spiess
#' @docType package
#' @name PCRedux-package
#' @aliases PCRedux
#' @examples
#' # Use the mblrr functopm to analyse 
#' mblrr(x=1:35, y=rnorm(35))

NULL
