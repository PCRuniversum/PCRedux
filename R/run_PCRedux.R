#' PCRedux app
#'
#' A graphical user interface for compute the properties of amplification curves.
#' @return null.
#' @export
#' @seealso \code{\link{encu}}, \code{\link[shiny]{runApp}}.
#' @importFrom shiny runApp
#' @note
#' Any ad-blocking software may cause malfunctions.

run_PCRedux <- function() {
  runApp(system.file("PCRedux-app", package = "PCRedux"))
}
