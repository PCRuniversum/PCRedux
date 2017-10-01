#' @title Performance analysis for binary classification
#' @description This function performs an analysis sensitivity and specificity 
#' to asses the performance of a binary classification test. For further reading
#' the studies by Brenner and Gefeller 1997 and by Kuhn 2008 are a good starting
#' point. 
#' @param data is a vector with logical decisions (0, 1) of the test system.
#' @param reference  is a vector with logical decisions (0, 1) of the reference 
#' system.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords sensitivity specificity precision accuracy
#' @seealso 
#'  \code{\link[caret]{sensitivity}}
#'  \code{\link[caret]{specificity}}
#'  \code{\link[caret]{negPredValue}}
#' @references H. Brenner, O. Gefeller, others, Variation of sensitivity, 
#' specificity, likelihood ratios and predictive values with disease prevalence, 
#' \emph{Statistics in Medicine}. 16 (1997) 981--991.
#'
#' M. Kuhn, Building Predictive Models in R Using the caret Package, 
#' \emph{Journal of Statistical Software}. 28 (2008). doi:10.18637/jss.v028.i05.
#' @details TP, true positive; FP, false positive; TN, true negative; FN, false 
#' negative
#' 
#' Sensitivity - TPR, true positive rate
#' TPR = TP / (TP + FN)
#' 
#' Specificity - SPC, true negative rate
#' SPC = TN / (TN + FP)
#' 
#' Precision - PPV, positive predictive value
#'  PPV = TP / (TP +  FP)
#'
#' Negative predictive value - NPV
#'  NPV = TN / (TN + FN)
#'
#' Fall-out, FPR, false positive rate
#'  FPR = FP / (FP + TN) = 1 - SPC
#'
#' False negative rate - FNR
#'  FNR = FN / (TN + FN) = 1 - TPR
#'
#' False discovery rate - FDR
#' FDR = FP / (TP + FP) = 1 - PPV
#'
#' Accuracy - ACC
#' ACC = (TP + TN) / (TP + FP + FN + TN)
#'
#' F1 score
#' F1 = 2TP / (2TP + FP + FN)
#'
#' Likelihood ratio positive - LRp
#' LRp = TPR/(1-SPC)
#'
#' Matthews correlation coefficient (MCC)
#' MCC = (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
#' @examples 
#' # Produce some arbitrary binary decisions data
#' # test_data is the new test or method that should be analyzed
#' # reference_data is the reference data set that should be analyzed
#' test_data <- c(0,0,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1) 
#' reference_data <- c(0,0,0,0,1,1,1,1,0,1,0,1,0,1,0,1,0,1,0,1,1,1,1,1)
#' 
#' # Plot the data of the decisions
#' plot(1:length(test_data), test_data, xlab="Sample", ylab="Decisions", 
#'      yaxt="n", pch=19)
#' axis(2, at=c(0,1), labels=c("negative", "positive"), las=2)
#' points(1:length(reference_data), reference_data, pch=1, cex=2, col="blue")
#' legend("topleft", c("Test data", "Reference data"), pch=c(19,1), 
#'         cex=c(1.5,1.5), bty="n", col=c("black","blue"))
#' 
#' # Do the statistical analysis with the performance function
#' performance(data=test_data, reference=reference_data)
#' @rdname performance
#' @export performance

performance <- function(data, reference) {
    data <- data.frame(d=data,r=reference)
    
    counts <- nrow(data)
    
    # TP, true positive
    # FP, false positive
    # TN, true negative
    # FN, false negative
    
    # Sensitivity - TPR, true positive rate
    # TPR = TP / (TP + FN)
    
    data_tmp <- data[data[, "r"] == 1, ]
    TP <- length(which(data_tmp[["d"]] == data_tmp[["r"]]))
    FN <- length(which(data_tmp[["d"]] != data_tmp[["r"]]))
    TPR <- TP / (TP + FN)
        
    # Specificity - SPC, true negative rate
    # SPC = TN / (TN + FP)
    data_tmp <- data[data[, "r"] == 0, ]
    TN <- length(which(data_tmp[["d"]] == data_tmp[["r"]]))
    FP <- length(which(data_tmp[["d"]] != data_tmp[["r"]]))
    SPC <- TN / (TN + FP)
    
    # Precision - PPV, positive predictive value
    # PPV = TP / (TP +  FP)
    PPV <- TP / (TP +  FP)
    
    # Negative predictive value - NPV
    # NPV = TN / (TN + FN)
    
    NPV <- TN / (TN + FN)
    
    # Fall-out, FPR, false positive rate
    # FPR = FP / (FP + TN) = 1 - SPC
    FPR <- 1 - SPC
    
    # False negative rate - FNR
    # FNR = FN / (TN + FN) = 1 - TPR
    FNR <- 1 - TPR
    
    # False discovery rate - FDR
    # FDR = FP / (TP + FP) = 1 - PPV
    FDR <- 1 - PPV
    
    # Accuracy - ACC
    # ACC = (TP + TN) / (TP + FP + FN + TN)
    ACC <- (TP + TN) / (TP + FP + FN + TN)
    
    # F1 score
    # F1 = 2TP / (2TP + FP + FN)
    F1 <- 2*TP / (2*TP + FP + FN)
    
    # Matthews correlation coefficient (MCC)
    # MCC = (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    MCC <- (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    
    # Likelihood ratio positive
    # LRp = TPR/(1-SPC)
    LRp <- TPR/(1 - SPC)

    # Combination of all results
    res <- data.frame(
        TPR=TPR,
        SPC=SPC,
        PPV=PPV,
        NPV=NPV,
        FPR=FPR,
        FNR=FNR,
        FDR=FDR,
        ACC=ACC,
        F1=F1,
        MCC=MCC,
        LRp=LRp,
        TP=TP,
        TN=TN,
        FP=FP,
        FN=FN,
        counts=counts
    )
    res
}
