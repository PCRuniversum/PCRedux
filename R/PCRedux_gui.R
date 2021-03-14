#' PCRedux Graphical User Interface
#'
#' Launches graphical user interface that generates numerous
#' features of a large amplification curve data set, similarly to
#' the \code{\link[PCRedux]{encu}} function.
#' @return No return value, called for side effects
#'
#' @section Warning : Any ad-blocking software may cause malfunctions.
#' @export PCRedux_gui
PCRedux_gui <- function()
  runApp(system.file("PCRedux_gui", package = "PCRedux"))
