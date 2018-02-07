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
#' @importFrom xray anomalies
#' @importFrom zoo as.zoo
#' @author Stefan Roediger, Michal Burdukiewcz, Andrej-Nikolai Spiess, Konstantin A. Blagodatskikh
#' @docType package
#' @name PCRedux-package
#' @aliases PCRedux
#' @examples
#' # Use the mblrr function to analyse amplification curves
#' library(qpcR)
#' mblrr(x=boggy[, 1], y=boggy[, 2])

l7 <- list(expr = "Fluo ~ c + (k1 * Cycles) + (k2 * Cycles^2) + (d - c)/((1 + exp(b * (log(Cycles) - log(e))))^f)", 
    fct = function (x, parm) 
    {
        b <- parm[1]
        c <- parm[2]
        d <- parm[3]
        e <- parm[4]
        f <- parm[5]
        k1 <- parm[6]
        k2 <- parm[7]
        c + (k1 * x) + (k2 * x^2) + (d - c)/((1 + exp(b * (log(x) - 
            log(e))))^f)
    }, ssFct = function (x, y) 
    {
        d <- max(y) + 0.001
        c <- min(y) - 0.001
        x2 <- x[y > 0]
        y2 <- y[y > 0]
        logitTrans <- log((d - y2)/(y2 - c))
        lmFit <- lm(logitTrans ~ log(x2))
        coefVec <- coef(lmFit)
        b <- coefVec[2]
        e <- exp(-coefVec[1]/b)
        f <- 1
        lmFit2 <- lm(y2[1:10] ~ x2[1:10])
        k1 <- coef(lmFit2)[2]
        k2 <- -0.01 * k1
        ssVal <- as.numeric(c(b, c, d, e, f, k1, k2))
        names(ssVal) <- l7$parnames
        return(ssVal)
    }, d1 = function (x, parm) 
    {
        b <- parm[1]
        c <- parm[2]
        d <- parm[3]
        e <- parm[4]
        f <- parm[5]
        k1 <- parm[6]
        k2 <- parm[7]
        k1 + 2 * k2 * x + b * (c - d) * e^-b * f * x^(-1 + b) * 
            (1 + e^-b * x^b)^(-1 - f)
    }, d2 = function (x, parm) 
    {
        b <- parm[1]
        c <- parm[2]
        d <- parm[3]
        e <- parm[4]
        f <- parm[5]
        k1 <- parm[6]
        k2 <- parm[7]
        (e^(-2 * b) * (1 + e^-b * x^b)^(-2 - f) * (-b * (c - 
            d) * f * x^b * (e^b + x^b) + 2 * k2 * x^2 * (e^b + 
            x^b)^2 * (1 + e^-b * x^b)^f + b^2 * (c - d) * f * 
            x^b * (e^b - f * x^b)))/x^2
    }, inv = function (y, parm) 
    {
        x <- 1:100
        fn <- function(x, parm) l7$fct(x, parm) - y
        uniroot(fn, interval = c(1, 100), parm)$root
    }, expr.grad = expression(c + (k1 * Cycles) + (k2 * Cycles^2) + 
        (d - c)/((1 + exp(b * (log(Cycles) - log(e))))^f)), inv.grad = NULL, 
    parnames = c("b", "c", "d", "e", "f", "k1", "k2"), name = "l7", 
    type = "seven-parameter log-logistic")

NULL
