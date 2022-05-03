#' PCRedux app
#'
#' A graphical user interface for computing the properties of amplification curves. 
#' Take a look at the vignette to learn more about the different ways to start the app.
#' @return null.
#' @export
#' @return No return value, called for side effects
#' @seealso \code{\link{encu}}, \code{\link[shiny]{runApp}}.
#' @importFrom shiny runApp
#' @note
#' Any ad-blocking software may cause malfunctions.

run_PCRedux <- function() {
  runApp(system.file("PCRedux_gui", package = "PCRedux"))
}
