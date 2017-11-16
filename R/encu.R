#' A function to calculate numerous features from amplification curve data 
#' from a quantitative PCR experiment.
#' 
#' encu (ENcode CUrves) is a function to calculate numerous 
#' features of a large amplification curve data set. 
#' The \code{\link[PCRedux]{pcrfit_single}} is 
#' performing the analysis for a single process.  
#'
#' 
#' @param data is the data set containing the cycles and fluorescence amplitudes.
#' @param detection_chemistry contains additional meta information about the 
#' detection chemistry (e.g., probes, intercalating dye) that was used.
#' @param device contains additional meta information about the qPCR system that 
#' was used.
#' @return The output of the encu function is 
#' identical to the \code{\link[PCRedux]{pcrfit_single}} function.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords slope intercept preprocessing normalization
#' @importFrom pbapply pblapply
#' @examples
#' 
#' # Calculate curve features of an amplification curve data. Note that not all 
#' # available CPU cores are used. If need set "all" to use all available cores.
#' # In this example the testdat data set from the qpcR package is used.
#' # The samples F1.1 and F1.2 are positive amplification curves. The samples 
#' # F1.3 and F1.4 are negative.
#' 
#' library(qpcR)
#' res_encu <- encu(testdat[, 1:4])
#' res_encu
#'  
#' @export

encu <- function(data, detection_chemistry = NA, device = NA) {
  # Determine the number of available cores and register them
  
  # Prepare the data for further processing
  # Normalize RFU values to the alpha quantiles (0.999)
  cycles <- data.frame(cycles=data[, 1])
  data_RFU <- data.frame(data[, -1])
  data_RFU_colnames <- colnames(data_RFU)
  data_RFU <- sapply(1L:ncol(data_RFU), function(i) {
    data_RFU[, i] / quantile(data_RFU[, i], 0.999, na.rm=TRUE)
  })
  colnames(data_RFU) <- data_RFU_colnames
  
  # just to shut RCHeck for NSE we define ith_cycle
  ith_cycle <- 1
  
  run_res <- do.call(cbind,
                     pblapply(1L:ncol(data_RFU), function(ith_run) pcrfit_single(data_RFU[, ith_run])))
  
  res <- cbind(runs = colnames(data_RFU), run_res, 
               detection_chemistry = detection_chemistry,
               device = device)
  
  rownames(res) <- NULL
  
  res
  
}
