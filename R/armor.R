#' armor: fetch errors gently
#'
#' \code{armor} is a helper function that catches errors and creates an output 
#' that can be used for further processing.
#' 
#' @return gives a \code{numeric} value (S3 class) as output for errors
#'
#' @param f is the function that needs armor.
#' @param n is the number of Zero repeats if a function fails.
#' @author Andrej Nikolai Spiess, Stefan Roediger
#' @keywords error
#' @seealso
#'  \code{\link[base:suppressMessages]{base::suppressMessages()}}
#'  \code{\link[base:inherits]{base::inherits()}}
#' @examples
#'
#' # Fetch the error from the diffQ function
#' require(MBmca)
#' # In the following the approximate derivative of the amplification curve data
#' # x <- RAS002[, 1] and y <- RAS002[, 2] is calculated by diffQ(). 
#' # This will not give an error.
#' x <- RAS002[, 1]
#' y <- RAS002[, 2]
#' armor_diffQ_passes <- armor(MBmca::diffQ(cbind(x, y), verbose = TRUE)$xy)
#' armor_diffQ_passes
#' #
#' # In the following the approximate derivative of the sequences x <- 1:40
#' # and y <- 1:40 is calculated by diffQ(). However, this will fail.
#' # This will give the "internal" error
#' # > 
#' # Error in list.res[[i]][[8]] : subscript out of bounds
#' # that is resolved to 0.
#' x <- 1:40
#' y <- 1:40
#' armor_diffQ_fails <- armor(MBmca::diffQ(cbind(x, y), verbose = TRUE)$xy)
#' armor_diffQ_fails
#' @export armor

armor <- function(f, n = 1)  {
  OUT <- try(suppressMessages(f), silent = TRUE)
  if (inherits(OUT, "try-error")) return(rep(0, n)) else return(OUT)
} 
