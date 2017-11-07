#' PCRedux - quantitative PCR machine learning helper tool
#'
#' @description \code{PCRedux} package is a toolbox for the analysis of sigmoid curve (qPCR) data.
#' 
#' @section Machine learning:
#' In machine learning and statistics, classification is the task to identify a new observation which is assigned to a set of categories. A foundation are training data set, which containing observations with known memberships of the categories. In the context of sigmoid amplification curves this could be an assignment into "negative", "ambiguous" or "positive" classes. Basically, a set of descriptors (features of the curvature) is need to perform an assignment to a class. The PCRedux package contains function for feature extraction and human rated amplification curve with classes. 
#' 
#' @importFrom bcp bcp
#' @importFrom changepoint cpt.meanvar
#' @importFrom chipPCR amptester bg.max CPP smoother
#' @importFrom ecp e.agglo
#' @importFrom fda.usc fdata
#' @importFrom FFTrees FFTrees
#' @importFrom graphics matplot par
#' @importFrom grDevices rainbow
#' @importFrom magrittr %>%
#' @importFrom MBmca diffQ diffQ2 mcaPeaks
#' @importFrom plotly ggplotly
#' @importFrom pracma polyarea
#' @importFrom qpcR AICc efficiency LRE mselect pcrfit sliwin takeoff
#' @importFrom robustbase lmrob
#' @importFrom stats coefficients confint cor.test lag lm median na.omit quantile
#' @importFrom testthat context test_that
#' @importFrom utils head tail
#' @importFrom utils data
#' @importFrom visdat vis_dat
#' @importFrom zoo as.zoo
#' @author Stefan Roediger, Michal Burdukiewcz, Andrej-Nikolai Spiess, Konstantin A. Blagodatskikh
#' @docType package
#' @name PCRedux-package
#' @aliases PCRedux
#' @examples
#' # Use the mblrr function to analyse amplification curves
#' library(qpcR)
#' mblrr(x=boggy[, 1], y=boggy[, 2])

NULL
