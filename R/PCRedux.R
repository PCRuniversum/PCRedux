#' PCRedux - quantitative PCR machine learning helper tool
#'
#' @description \code{PCRedux} package is a toolbox for the analysis of sigmoid curve (qPCR) data.
#' 
#' @section Machine learning:
#' In machine learning and statistics, the classification should be used to identify a new unknown observation. This observation is assigned to a number of categories. One basis is training data sets containing observations with known classes. Using the example of sigmoid amplification curves, this could be an assignment to the class "negative","ambiguous" or "positive". Basically, a number of descriptors (e. g., characteristics of curvature) are required to be able to assign classes. This package contains functions for extracting characteristics. In addition, the package contains data sets of classified amplification curves.
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
