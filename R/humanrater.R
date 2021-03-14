#' Human Rater 2.0
#'
#' Launches graphical user interface for the manual annotation of
#' large amplification curve data sets, similarly to
#' the \code{\link[chipPCR]{humanrater}} function.
#' @return No return value, called for side effects
#'
#' @section Warning : Any ad-blocking software may cause malfunctions.
#' @export humanrater2
humanrater2 <- function()
  runApp(system.file("humanrater2", package = "PCRedux"))
