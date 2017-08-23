#' A helper function to plot the frequency of classes of amplification curves
#' 
#' \code{class_plotter} is a function to to plot the frequency of classes of 
#' amplification curves.
#' 
#' @param data is the data set containing the amplification curves.
#' @param classes is the data set containing the classes of the amplification curves.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords plot classes
#' @examples
#' 
#' # Use class_plotter to visualize the classified amplification curve data
#' library(qpcR)
#' library(data.table)
#' filename <- system.file("decission_res_htPCR.csv", package = "PCRedux")
#' class_plotter(data=htPCR, classes=as.data.frame(fread(filename)))
#' 
#' @export class_plotter


class_plotter <- function(data, classes) {
    
    range_ylim <- range(data[, -1])
    
    positive_classes    <- classes[classes[["test.result.3"]]=="y", ][, 1]
    ambiguous_classes   <- classes[classes[["test.result.3"]]=="a", ][, 1]
    negative_classes    <- classes[classes[["test.result.3"]]=="n", ][, 1]
    
    colors_positive_classes     <- rainbow(length(positive_classes), alpha=0.3)
    colors_ambiguous_classes    <- rainbow(length(ambiguous_classes), alpha=0.3) 
    colors_negative_classes     <- rainbow(length(negative_classes), alpha=0.3)
    
    par(mfrow=c(1,3))
    
    try(matplot(data[, 1], data[, colnames(data) %in% positive_classes], type="l", 
            col=colors_positive_classes, lty=1, main="Positive", xlab="Cycle", 
            lwd=1.5, ylab="RFU", ylim = c(range_ylim[1], range_ylim[2])), silent=TRUE)
    try(matplot(data[, 1], data[, colnames(data) %in% ambiguous_classes], type="l", 
            col=colors_ambiguous_classes, lty=1, main="Ambiguous", xlab="Cycle", 
            lwd=1.5, ylab="RFU", ylim = c(range_ylim[1], range_ylim[2])), silent=TRUE)
    try(matplot(data[, 1], data[, colnames(data) %in% negative_classes], type="l", 
            col=colors_negative_classes, lty=1, main="Negative", xlab="Cycle", 
            lwd=1.5, ylab="RFU", ylim = c(range_ylim[1], range_ylim[2])), silent=TRUE)
}
