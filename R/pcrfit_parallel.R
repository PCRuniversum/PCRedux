#' A function to calculate numerous features from amplification curve data 
#' from a quantitative PCR experiment.
#' 
#' \code{pcrfit_parallel} is a function to calculate numerous features of a large 
#' amplification curve data set. pcrfit_parallel makes use of parallelized code
#' to make use of multi-core architectures.
#' 
#' @param data is the data set containing the cycles and fluorescence amplitudes.
#' @param n_cores defines the numbers of cores that should be left unused by this 
#' function. By default \code{pcrfit_parallel} is using only two cores. Use
#' \code{"all"} to use all available cores.
#' @param detection_chemistry contains additional meta information about the 
#' detection chemistry (e.g., probes, intercalating dye) that was used.
#' @param device contains additional meta information about the qPCR system that 
#' was used.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords slope intercept
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach
#' @examples
#' 
#' # Calculate curve features of an amplification curve data set by using all 
#' # available CPU cores.
#' library(qpcR)
#' res_pcrfit_parallel <- pcrfit_parallel(boggy[, 1:7])
#' res_pcrfit_parallel
#'  
#' @export pcrfit_parallel

pcrfit_parallel <- function(data, n_cores = 1, detection_chemistry = NA, device = NA) {
    # Determine the number of available cores and registrate them    
    if(n_cores == "all")
      n_cores <- detectCores() 
    
    registerDoParallel(n_cores)
    
    # Prepare the data for further processing
    # Normalize RFU values to the alpha quantiles (0.999)
    cycles <- data.frame(cycles=data[, 1])
    data_RFU <- data.frame(data[, -1])
    data_RFU_colnames <- colnames(data_RFU)
    data_RFU <- sapply(1L:ncol(data_RFU), function(i) {
        data_RFU[, i] / quantile(data_RFU[, i], 0.999)
    })
    colnames(data_RFU) <- data_RFU_colnames

    # just to shut RCHeck for NSE we define ith_cycle
    ith_cycle <- 1
    run_res <- foreach(ith_cycle = 1L:ncol(data_RFU), 
            .packages=c("bcp", "changepoint", "chipPCR", "ecp", "MBmca",
                        "PCRedux", "pracma", "qpcR", "robustbase", 
                        "zoo"),
            .combine = rbind) %dopar% {
                          prcfit_single(data_RFU[, ith_cycle])
                        }

    res <- cbind(runs = colnames(data_RFU), run_res, 
                 detection_chemistry = detection_chemistry,
                 device = device)
    
    res
    
}
