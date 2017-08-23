#' A function to calculate numerous features from amplification curve data 
#' from a quantitative PCR experiment.
#' 
#' \code{pcrfit_parallel} is a function to calculate numerous features of a large 
#' amplification curve data set. pcrfit_parallel makes use of parallelized code
#' to make use of multi-core architectures.
#' 
#' @param data is the data set containing the cycles and fluorescence amplitudes.
#' @param less_cores defines the numbers of cores that should be left unused by this function. By default is pcrfit_parallel using all cores.
#' @param detection_chemistry contains additional meta information about the detection chemistry (e.g., probes, intercalating dye) that was used.
#' @param device contains additional meta information about the qPCR system that was used.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords slope intercept
#' @examples
#' 
#' # Calculate curve features of an amplification curve data set by using all 
#' # available CPU cores.
#'
#' # T.B.D.
#' 
#' @export pcrfit_parallel

pcrfit_parallel <- function(data, less_cores=0, detection_chemistry=NA, device=NA) {
    # Determine the number of available cores and registrate them    
    n_cores <- detectCores() - less_cores
    if(n_cores < 2) {n_cores <- 2} 
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
    
    # Cut the data.frame into smaller blocks for parallel analysis
    if(ncol(data_RFU) <= 3) {
        registerDoParallel(1)
        cuts <- rep(1, ncol(data_RFU))
    } else {
        cuts <- cut(1L:ncol(data_RFU), breaks = n_cores, label=FALSE)
    }
    
    # Define names of feature labels used for the ML analysis
    features <- c("sample",
                  "eff",
                  "cpD1",
                  "cpD2",
                  "fluo",
                  "init1",
                  "init2",
                  "top", 
                  "f.top",
                  "resLRE",
                  "ressliwin",
                  "cpDdiff",
                  "slope_background",
                  "intercept_background",
                  "polyarea",
                  "changepoint.e.agglo",
                  "changepoint.cpt.mean",
                  "changepoint.bcp",
                  "qPCRmodel",
                  "amptester_shap.noisy",
                  "amptester_lrt.test",
                  "amptester_rgt.dec",
                  "amptester_tht.dec",
                  "amptester_slt.dec",
                  "amptester_polygon",
                  "amptester_slope.ratio",
                  "minRFU", 
                  "maxRFU",
                  "head_to_tail_ratio",
                  "bg.start_normalized",
                  "bg.stop_normalized",
                  "amp.stop_normalized",
                  "autocorellation",
                  "detection_chemistry",
                  "device",
                  "mblrr_intercept_less",
                  "mblrr_slope_less",
                  "mblrr_cor_less",
                  "mblrr_intercept_more",
                  "mblrr_slope_more",
                  "mblrr_cor_more")
    
    # Do the parallel analysis
    res_tmp <- foreach(block=unique(cuts), 
                       .packages=c("bcp", "changepoint", "chipPCR", "ecp", 
                                   "pracma", "qpcR", "robustbase", "zoo"), 
                       .combine=cbind, .export=c("simple_fn")) %dopar% {
        
        # Combine the cycle values and a specific data block that was previously 
        # cutted from the input data
        # Select the appropriate data for a further analysis
        dat <- data.frame(cycles, data_RFU[, which(cuts == block)])
        column_names <- colnames(dat)
        
        res_block <- sapply(2L:ncol(dat), function(bc) {
            
            # Determine highest and lowest amplification curve values
            fluo_range <- quantile(dat[, bc], c(0.01, 0.99))
            
            res_bg.max_tmp <- bg.max(unlist(cycles), dat[, bc])
            length_cycle <- nrow(cycles)
            res_bg.max <- c(bg.start=res_bg.max_tmp@bg.start/length_cycle,
                                bg.stop=res_bg.max_tmp@bg.stop/length_cycle,
                                amp.stop=res_bg.max_tmp@amp.stop/length_cycle)
 
            
            # Determine the head to tail ratio
            res_head_tail_ratio <- head2tailratio(dat[, bc])
            
            # Determine the slope from the cycles 2 to 11
            res_lm_fit <- earlyreg(x= dat[, 1], dat[, bc])
            
            # Calculates the area of the amplification curve
            res_polyarea <- try(polyarea(dat[, 1], dat[, bc]), silent=TRUE)
            if(class(res_polyarea) == "try-error") {res_polyarea <- NA}
            
            # Calculate change points
            # Agglomerative hierarchical estimation algorithm for multiple change point analysis
            res_changepoint_e.agglo <- try(length(e.agglo(as.matrix(dat[, bc]))$estimates)/length_cycle, silent=TRUE)
            if(class(res_changepoint_e.agglo) == "try-error") {res_changepoint_e.agglo <- NA}
            
            # Binary Segmentation
            res_changepoint_cpt.mean <- try(length(cpt.meanvar(as.matrix(dat[, bc]), penalty="CROPS",method="PELT")@param.est)/length_cycle, silent=TRUE)
            if(class(res_changepoint_cpt.mean) == "try-error") {res_changepoint_cpt.mean <- NA}
            
            # Bayesian analysis of change points
            res_bcp_tmp <- bcp(dat[, bc])            
            res_bcp_tmp <- res_bcp_tmp$posterior.prob > 0.45
            res_bcp <- try((which(as.factor(res_bcp_tmp) == TRUE) %>% length)/length_cycle)
            if(class(res_bcp) == "try-error") {res_bcp <- NA}
            
            # Median based local robust regression (mblrr)
            res_mblrr <- mblrr(dat[, 1], dat[, bc])

            # Calculate amptester results
            res_amptester <- try(amptester(dat[, bc]))
            
            # Perform an autocorrelation analysis           
            res_autocorrelation <- autocorrelation_test(y=dat[, bc])
            
            # Fit sigmoidal models to curve data
            res_fit <- try(mselect(pcrfit(dat, cyc=1, fluo=bc), 
                                   verbose=FALSE, do.all=TRUE), silent=TRUE)

            # Determine the model suggested by the mselect function based on the AICc
            res_fit_model <- try(names(which(res_fit[["retMat"]][, "AICc"] == min(res_fit[["retMat"]][, "AICc"]))), silent=TRUE)
            if(class(res_fit_model) == "try-error") {res_fit_model <- NA}
            
            if(class(res_fit)[1] != "try-error") {
                # TakeOff Point
                # Calculates the first significant cycle of the exponential region 
                # (takeoff point) using externally studentized residuals as described 
                # in Tichopad et al. (2003).
                res_takeoff <- try(takeoff(res_fit), silent=TRUE)
                if(class(res_takeoff) == "try-error") {res_takeoff <- list(NA, NA)}
                
                # LRE qPCR efficiency
                # Calculation of qPCR efficiency by the 'linear regression of 
                # efficiency' method
                
                res_LRE <- try(LRE(res_fit, wsize=5, border=c(0,3), base=0, 
                                   plot=FALSE, verbose=FALSE)$eff, silent=TRUE)
                if(class(res_LRE) == "try-error") {res_LRE <- NA}
                
                # sliwin qPCR efficiency
                # Calculation of the qPCR efficiency by the 'window-of-linearity' method
                res_sliwin <- try(sliwin(res_fit, wsize=5, plot=FALSE, verbose=FALSE)$eff, 
                                  silent=TRUE)
                if(class(res_sliwin) == "try-error") {res_sliwin <- NA}
                
                # Cq of the amplification curve
                # Determine the Cq and other parameters
                res_efficiency_tmp <- try(
                    efficiency(res_fit, plot=FALSE)[c("eff",
                                                      "cpD1", "cpD2",
                                                      "fluo",
                                                      "init1", "init2")], 
                                          silent=TRUE)
                res_cpDdiff <- try(res_efficiency_tmp[["cpD1"]] - res_efficiency_tmp[["cpD2"]])
                
            } else {
                res_efficiency_tmp <- list(NA, NA, NA, NA, NA, NA)
                res_takeoff <- list(NA, NA)
                res_LRE <- NA
                res_sliwin <- NA
                res_cpDdiff <- NA
            }
  
            res_efficiency <- list(
                    sample=column_names[bc],
                    eff=res_efficiency_tmp[["eff"]],
                    cpD1=res_efficiency_tmp[["cpD1"]],
                    cpD2=res_efficiency_tmp[["cpD2"]],
                    fluo=res_efficiency_tmp[["fluo"]],
                    init1=res_efficiency_tmp[["init1"]],
                    init2=res_efficiency_tmp[["init2"]],
                    top=res_takeoff[[1]], 
                    f.top=res_takeoff[[2]],
                    resLRE=res_LRE[1],
                    ressliwin=res_sliwin[[1]],
                    cpDdiff=res_cpDdiff,
                    slope_background=res_lm_fit[["slope"]],
                    intercept_background=res_lm_fit[["intercept"]],
                    polyarea=res_polyarea,
                    changepoint.e.agglo=res_changepoint_e.agglo,
                    changepoint.cpt.mean=res_changepoint_cpt.mean,
                    changepoint.bcp=res_bcp,
                    qPCRmodel=res_fit_model[[1]],
                    amptester_shap.noisy=res_amptester@decisions["shap.noisy"][[1]],
                    amptester_lrt.test=res_amptester@decisions["lrt.test"][[1]],
                    amptester_rgt.dec=res_amptester@decisions["rgt.dec"][[1]],
                    amptester_tht.dec=res_amptester@decisions["tht.dec"][[1]],
                    amptester_slt.dec=res_amptester@decisions["slt.dec"][[1]],
                    amptester_polygon=res_amptester@"polygon",
                    amptester_slope.ratio=res_amptester@"slope.ratio",
                    minRFU=fluo_range[[1]], 
                    maxRFU=fluo_range[[2]],
                    bg.start_normalized=res_bg.max[1],
                    bg.stop_normalized=res_bg.max[2],
                    amp.stop_normalized=res_bg.max[3],
                    head_to_tail_ratio=res_head_tail_ratio,
                    autocorellation=res_autocorrelation,
                    detection_chemistry=detection_chemistry,
                    device=device,
                    mblrr_intercept_less=res_mblrr[1],
                    mblrr_slope_less=res_mblrr[2],
                    mblrr_cor_less=res_mblrr[3],
                    mblrr_intercept_more=res_mblrr[4],
                    mblrr_slope_more=res_mblrr[5],
                    mblrr_cor_more=res_mblrr[6]                    
                )
        })
    }
}
