#' @title Performance analysis for binary classification
#' @description This function performs an analysis sensitivity and specificity
#' to asses the performance of a binary classification test. For further reading
#' the studies by Brenner and Gefeller 1997, James 2013 by Kuhn 2008 are a good
#' starting point.
#'
#' @return gives a \code{data.frame} (S3 class, type of \code{list}) as output 
#' for the performance
#'
#' @param sample is a vector with logical decisions (0, 1) of the test system.
#' @param reference is a vector with logical decisions (0, 1) of the reference
#' system.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords sensitivity specificity precision accuracy
#' @references H. Brenner, O. Gefeller, others, Variation of sensitivity,
#' specificity, likelihood ratios and predictive values with disease prevalence,
#' \emph{Statistics in Medicine}. 16 (1997) 981--991.
#'
#' M. Kuhn, Building Predictive Models in R Using the caret Package,
#' \emph{Journal of Statistical Software}. 28 (2008). doi:10.18637/jss.v028.i05.
#'
#' G. James, D. Witten, T. Hastie, R. Tibshirani, An Introduction to Statistical
#' Learning, \emph{Springer New York, New York, NY}, (2013).
#' doi:10.1007/978-1-4614-7138-7.
#'
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
#' MCC = (TP*TN - FP*FN) /
#        (sqrt(TP + FP) * sqrt(TP + FN) *
#'        sqrt(TN + FP) * sqrt(TN+FN)
#'       )
#'
#' Cohen's kappa (binary classification)
#' kappa=(p0-pc)/(1-p0)
#'
#' r (reference) is the trusted label and s (sample) is the predicted value
#'
#' \tabular{ccc}{
#'       \tab r=1 \tab r=0 \cr
#'   s=1 \tab a   \tab b   \cr
#'   s=0 \tab c   \tab d
#' }
#'
#' \deqn{n = a + b + c + d}
#'
#' pc=((a+b)/n)((a+c)/n)+((c+d)/n)((b+d)/n)
#'
#' po=(a+d)/n
#'
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
#' legend("topleft", c("Sample", "Reference"), pch=c(19,1),
#'         cex=c(1.5,1.5), bty="n", col=c("black","blue"))
#'
#' # Do the statistical analysis with the performeR function
#' performeR(sample=test_data, reference=reference_data)
#' @rdname performeR
#' @export performeR

performeR <- function(sample, reference) {
  data <- data.frame(s = sample, r = reference)


  # TP, true positive
  # FP, false positive
  # TN, true negative
  # FN, false negative

  # Sensitivity - TPR, true positive rate
  # TPR = TP / (TP + FN)

  data_tmp <- data[data[, "r"] == 1, ]
  TP <- length(which(data_tmp[["s"]] == data_tmp[["r"]]))
  FN <- length(which(data_tmp[["s"]] != data_tmp[["r"]]))
  TPR <- TP / (TP + FN)

  # Specificity - SPC, true negative rate
  # SPC = TN / (TN + FP)
  data_tmp <- data[data[, "r"] == 0, ]
  TN <- length(which(data_tmp[["s"]] == data_tmp[["r"]]))
  FP <- length(which(data_tmp[["s"]] != data_tmp[["r"]]))
  SPC <- TN / (TN + FP)

  # Precision - PPV, positive predictive value
  # PPV = TP / (TP +  FP)
  PPV <- TP / (TP + FP)

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
  F1 <- 2 * TP / (2 * TP + FP + FN)

  # Matthews correlation coefficient (MCC)
  # MCC = (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  MCC <- (TP * TN - FP * FN) / (sqrt(TP + FP) * sqrt(TP + FN) * sqrt(TN + FP) * sqrt(TN + FN))

  # Likelihood ratio positive
  # LRp = TPR/(1-SPC)
  LRp <- TPR / (1 - SPC)

  # Number of samples for analysis

  n <- nrow(data)
  if (TP + TN + FP + FN == n) {
    counts <- n
  } else {
    message("error")
  }

  # Cohen's kappa coefficient
  # r is the trusted label and s is the predicted value
  #
  #      | r= 1| r=0
  # -----|-----|-----
  # s=1  | a   | b
  # s=0  | c   | d
  #
  # n=a+b+c+d
  #
  # pc=((a+b)/n)((a+c)/n)+((c+d)/n)((b+d)/n)
  #
  # po=(a+d)/n
  #
  # k=(po-pc)/(1-pc)

  a <- TP
  b <- FP
  c <- FN
  d <- TN

  pc <- ((a + b) / n) * ((a + c) / n) + ((c + d) / n) * ((b + d) / n)
  po <- (a + d) / n

  kappa <- (po - pc) / (1 - pc)

  # Combination of all results
  res <- data.frame(
    TPR = TPR,
    SPC = SPC,
    PPV = PPV,
    NPV = NPV,
    FPR = FPR,
    FNR = FNR,
    FDR = FDR,
    ACC = ACC,
    F1 = F1,
    MCC = MCC,
    LRp = LRp,
    kappa = kappa,
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN,
    counts = counts
  )
  res
}
