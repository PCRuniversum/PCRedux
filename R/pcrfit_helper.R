#' @title prcfit_single - A function to extract features from an amplification curve
#' @description The prcfit_single is responsible for the 
#' extraction of features from amplification curve data. This function is the basis
#' of the \code{\link[PCRedux]{pcrfit_parallel}} function. The later performs the 
#' parallelized analysis of amplification curve data
#' @param x is the data set containing the fluorescence amplitudes.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' # Load the chipPCR package and analyse from the C126EG685 the first qPCR run
#' # "A01" (column 2).
#' library(chipPCR)
#' res <- prcfit_single(C126EG685[, 2])
#' @seealso 
#'  \code{\link[bcp]{bcp}}
#'  \code{\link[chipPCR]{bg.max}},\code{\link[chipPCR]{amptester}},\code{\link[chipPCR]{smoother}}
#'  \code{\link[ecp]{e.agglo}}
#'  \code{\link[MBmca]{diffQ}},\code{\link[MBmca]{mcaPeaks}},\code{\link[MBmca]{diffQ2}}
#'  \code{\link{head2tailratio}},\code{\link{earlyreg}},\code{\link{hookreg}},\code{\link{mblrr}},\code{\link{autocorrelation_test}}
#'  \code{\link[pracma]{polyarea}}
#'  \code{\link[qpcR]{pcrfit}},\code{\link[qpcR]{takeoff}},\code{\link[qpcR]{LRE}},\code{\link[qpcR]{sliwin}},\code{\link[qpcR]{efficiency}}
#'  \code{\link[base]{diff}}
#'  \code{\link[stats]{quantile}}
#' @rdname prcfit_single
#' @export prcfit_single

prcfit_single <- function(x) {

  # Normalize RFU values to the alpha quantiles (0.999)
  x <- x/quantile(x, 0.999)
  length_cycle <- length(x)
  cycles <- 1L:length_cycle
  # Determine highest and lowest amplification curve values
  fluo_range <- stats::quantile(x, c(0.01, 0.99), na.rm=TRUE)

  # for qpcR
  dat <- cbind(cyc=cycles, fluo=x)

  res_bg.max_tmp <- chipPCR::bg.max(cycles, x)
  res_bg.max <- c(bg.start=res_bg.max_tmp@bg.start/length_cycle,
                  bg.stop=res_bg.max_tmp@bg.stop/length_cycle,
                  amp.stop=res_bg.max_tmp@amp.stop/length_cycle)

  # Determine the head to tail ratio
  res_head_tail_ratio <- PCRedux::head2tailratio(x)

  # Determine the slope from the cycles 2 to 11
  res_lm_fit <- PCRedux::earlyreg(x=cycles, x)

  # Try to estimate the slope and the intercept of the tail region,
  # which might be indicative of a hook effect (strong negative slope)
  res_hookreg <- PCRedux::hookreg(x=cycles, x)

  # Calculates the area of the amplification curve
  res_polyarea <- try(pracma::polyarea(cycles, x), silent=TRUE)
  if(class(res_polyarea) == "try-error") {res_polyarea <- NA}
  
  # Calculate change points
  # Agglomerative hierarchical estimation algorithm for multiple change point analysis
  res_changepoint_e.agglo <- try(length(ecp::e.agglo(as.matrix(x))$estimates), silent=TRUE)
  if(class(res_changepoint_e.agglo) == "try-error") {res_changepoint_e.agglo <- NA}

  # Bayesian analysis of change points
  res_bcp_tmp <- bcp::bcp(x)            
  res_bcp_tmp <- res_bcp_tmp$posterior.prob > 0.45
  res_bcp <- try((which(as.factor(res_bcp_tmp) == TRUE) %>% length))
  if(class(res_bcp) == "try-error") {res_bcp <- NA}

  # Median based local robust regression (mblrr)
  res_mblrr <- PCRedux::mblrr(cycles, x)

  # Calculate amptester results
  res_amptester <- try(chipPCR::amptester(x))
  
  # Estimate the spread of the approximate local minima and maxima of the curve data
  dat_smoothed <- chipPCR::smoother(cycles, x)

  res_diffQ <- suppressMessages(MBmca::diffQ(cbind(cycles, dat_smoothed), verbose = TRUE)$xy)
  res_mcaPeaks <- MBmca::mcaPeaks(res_diffQ[, 1], res_diffQ[, 2])
  mcaPeaks_minima_maxima_ratio <- base::diff(range(res_mcaPeaks$p.max[, 2])) / diff(range(res_mcaPeaks$p.min[, 2]))

  # Estimate the slope between the minimum and the maximum of the second derivative
  res_diffQ2 <- suppressMessages(MBmca::diffQ2(cbind(cycles, dat_smoothed), verbose=FALSE, fct=min))
  range_Cq <- diff(res_diffQ2[[3]])
  if(res_diffQ2[[3]][1] < res_diffQ2[1] && res_diffQ2[1] < res_diffQ2[[3]][2]) {range_Cq} else {range_Cq <- 0}
  if(res_diffQ2[[3]][1] < res_diffQ2[1] && res_diffQ2[1] < res_diffQ2[[3]][2] && range_Cq > 1 && range_Cq < 9) {
    res_diffQ2_slope <- coefficients(lm(unlist(c(res_diffQ2[[4]])) ~ unlist(c(res_diffQ2[[3]]))))[2]
  } else {res_diffQ2_slope <- 0}

  # Perform an autocorrelation analysis
  res_autocorrelation <- PCRedux::autocorrelation_test(y=x)

  browser()
  # Fit sigmoidal models to curve data
  res_fit <- try(mselect(qpcR::pcrfit(dat, 1, 2, model = qpcR::l4), 
                         verbose=FALSE, do.all=TRUE), silent=TRUE)

  # Determine the model suggested by the mselect function based on the AICc
  res_fit_model <- try(names(which(res_fit[["retMat"]][, "AICc"] == min(res_fit[["retMat"]][, "AICc"]))), silent=TRUE)
  if(class(res_fit_model) == "try-error") {res_fit_model <- NA}

  if(class(res_fit)[1] != "try-error") {
    # TakeOff Point
    # Calculates the first significant cycle of the exponential region 
    # (takeoff point) using externally studentized residuals as described 
    # in Tichopad et al. (2003).
    res_takeoff <- try(qpcR::takeoff(res_fit), silent=TRUE)
    if(class(res_takeoff) == "try-error") {res_takeoff <- list(NA, NA)}

    # LRE qPCR efficiency
    # Calculation of qPCR efficiency by the 'linear regression of 
    # efficiency' method

    res_LRE <- try(qpcR::LRE(res_fit, plot=FALSE, verbose=FALSE)$eff, silent=TRUE)
    if(class(res_LRE) == "try-error") {res_LRE <- NA}

    # sliwin qPCR efficiency
    # Calculation of the qPCR efficiency by the 'window-of-linearity' method
    res_sliwin <- try(qpcR::sliwin(res_fit, plot=FALSE, verbose=FALSE)$eff, 
                      silent=TRUE)
    if(class(res_sliwin) == "try-error") {res_sliwin <- NA}

    # Cq of the amplification curve
    # Determine the Cq and other parameters
    res_efficiency_tmp <- try(
      qpcR::efficiency(res_fit, plot=FALSE)[c("eff",
                                        "cpD1", "cpD2",
                                        "fluo",
                                        "init1", "init2")],
      silent=TRUE)
    if(class(res_efficiency_tmp) != "try-error") {
      res_cpDdiff <- try(res_efficiency_tmp[["cpD1"]] - res_efficiency_tmp[["cpD2"]])
    } else {
      res_efficiency_tmp <- list(eff = NA,
                                 cpD1 = NA,
                                 cpD2 = NA,
                                 fluo = NA,
                                 init1 = NA,
                                 init2 = NA)
      res_cpDdiff <- NA
    }
  } else {
    res_efficiency_tmp <- list(eff = NA,
                               cpD1 = NA,
                               cpD2 = NA,
                               fluo = NA,
                               init1 = NA,
                               init2 = NA)
    res_takeoff <- list(NA, NA)
    res_LRE <- NA
    res_sliwin <- NA
    res_cpDdiff <- NA
  }
  
  res_efficiency <- data.frame(
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
    mblrr_intercept_less=res_mblrr[1],
    mblrr_slope_less=res_mblrr[2],
    mblrr_cor_less=res_mblrr[3],
    mblrr_intercept_more=res_mblrr[4],
    mblrr_slope_more=res_mblrr[5],
    mblrr_cor_more=res_mblrr[6],
    hookreg_hook=res_hookreg[["hook"]],
    mcaPeaks_minima_maxima_ratio=mcaPeaks_minima_maxima_ratio,
    diffQ2_slope=res_diffQ2_slope,
    diffQ2_Cq_range=range_Cq
  )
}
