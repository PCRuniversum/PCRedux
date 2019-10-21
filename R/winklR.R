#' winklR: A function to calculate the angle based on the first and the second 
#' derivative of an amplification curve data from a quantitative PCR experiment.
#'
#' \code{winklR} is a function to calculate the in the trajectory of the first
#' and the second derivatives maxima and minima  of an amplification curve data
#' from a quantitative PCR experiment. For the determination of the angle 
#' (central angle), the origin is the maximum of the first derivative. On this 
#' basis, the vectors to the minimum and maximum of the second figure are 
#' determined. The vectors result from the relation of the maximum of the 
#' first derivative to the minimum of the second derivative and from the 
#' maximum of the first derivative to the maximum of the second derivative. 
#' In a simple trigonometric approach, the scalar product of the two vectors 
#' is formed first. Then the absolute values are calculated and multiplied by 
#' each other. Finally, the value is converted into an angle with the cosine. The
#' assumption is that flat (negative amplification curves) have a large angle
#' and sigmoid (positive amplification curves) have a smaller angle. Another
#' assumption is that this angle is independent of the rotation of the
#' amplification curve. This means that systematic off-sets, such as those
#' caused by incorrect background correction, are of no consequence.
#' The cycles to be analyzed is defined by the user.
#' The output contains the angle.
#'
#' @param x is the cycle numbers (x-axis). By default the first ten cycles are removed.
#' @param y is the cycle dependent fluorescence amplitude (y-axis).
#' @param normalize is a logical parameter, which indicates if the amplification curve
#' data should be normalized to the 99 percent percentile of the amplification curve.
#' @param preprocess is a logical parameter, which indicates if the amplification curve
#' data should be smoothed (moving average filter, useful for noisy, jagged data).
#' @author Stefan Roediger
#' @keywords angle derivative
#' @seealso
#'  \code{\link[base]{acos}}
#'  \code{\link[MBmca]{diffQ2}}
#' @examples
#'
#' # Calculate the angles for amplification curve data from the RAS002 data set
#' data(RAS002)
#'
#' # Plot the data
#' plot(RAS002[, 1],
#'   y = RAS002[, 2], xlab = "Cycle", ylab = "RFU",
#'   main = "RAS002 data set", lty = 1, type = "l"
#' )
#' res <- winklR(x = RAS002[, 1], y = RAS002[, 2])
#' res
#' plot(rbind(res$origin, res$p1, res$p2), col = c("black", "green", "blue"))
#'
#' plot(RAS002[, 1],
#'   y = RAS002[, 7], xlab = "Cycle", ylab = "RFU",
#'   main = "RAS002 data set", lty = 1, type = "l"
#' )
#' res <- winklR(x = RAS002[, 1], y = RAS002[, 7])
#' res
#' plot(rbind(res$origin, res$p1, res$p2), col = c("black", "green", "blue"))
#'
#' res_angles <- unlist(lapply(2:21, function(i) {
#'   winklR(RAS002[, 1], RAS002[, i])$angle
#' }))
#' cdplot(RAS002_decisions[1L:20] ~ res_angles, xlab = "angle", ylab = "decision")
#' @export winklR


winklR <- function(x, y, normalize = FALSE, preprocess = TRUE) {
  data <- na.omit(cbind(x = x, y = y))

  x <- data[, "x"]

  if (is.integer(x) == FALSE) {
    x <- 1L:length(x)
  }

  y <- data[, "y"]

  # Optinal normalize data
  if (normalize) {
    y <- y / quantile(y, 0.99)
  }

  # Smooth data with Savitzky-Golay smoothing filter for other data
  # analysis steps.
  if (preprocess) {
    y <- chipPCR::smoother(x, y, method = "mova")
  } else {
    y <- data[, "y"]
  }

  guess_direction <- ifelse(median(head(y, 5)) > (median(tail(y, 5)) + mad(tail(y, 5))), 'max', 'min')
  
  # Calculate the point of the first and the second derivatives
  res <- try(suppressMessages(MBmca::diffQ2(cbind(x[-c(1:10)], y[-c(1:10)]),
    inder = TRUE, verbose = TRUE,
    fct = get(guess_direction, pos = "package:base")
  )), silent = TRUE)

  origin <- data.frame(res[["TmD1"]][1], res[["TmD1"]][2])
  point_x1 <- data.frame(res[["xTm1.2.D2"]][1], res[["yTm1.2.D2"]][1])
  point_x2 <- data.frame(res[["xTm1.2.D2"]][2], res[["yTm1.2.D2"]][2])

  # Calculate the distance between the points
  #distance_a <- sqrt((origin[1] - res[["xTm1.2.D2"]][1])^2 + (origin[2] - res[["yTm1.2.D2"]][1])^2)
  #distance_b <- sqrt((origin[1] - res[["xTm1.2.D2"]][2])^2 + (origin[2] - res[["yTm1.2.D2"]][2])^2)

  # Calculation of vectors (origin is the starting point)
  u <- data.frame(
    u_x = origin[1] - point_x1[1],
    u_y = origin[2] - point_x1[2]
  )

  v <- data.frame(
    v_x = origin[1] - point_x2[1],
    V_y = origin[2] - point_x2[2]
  )


  # Determine the dot product and the absolute values
  dot_product <- sum(u * v)

  length_of_vectors <- sqrt(sum(u^2)) * sqrt(sum(v^2))

  # Calculate angle

  # angle <- acos(dot_product / length_of_vectors) * 180 / pi

  angle <- dot_product / length_of_vectors

  rownames(origin) <- "origin"
  colnames(origin) <- c("x", "y")

  rownames(point_x1) <- "p1"
  colnames(point_x1) <- c("x", "y")
  rownames(point_x2) <- "p2"
  colnames(point_x2) <- c("x", "y")

  rownames(u) <- "u"
  rownames(v) <- "v"

  colnames(u) <- c("x", "y")
  colnames(v) <- c("x", "y")



  output <- list(angle, origin, point_x1, point_x2)

  names(output) <- c("angle", "origin", "p1", "p2")
  output
}
