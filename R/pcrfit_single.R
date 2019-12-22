#' @title pcrfit_single - A function to extract features from an amplification curve
#' @description The pcrfit_single is responsible for the
#' extraction of features from amplification curve data. The function can be used
#' for custom functions for a paralleled analysis of amplification curve data.
#' An example is given in the vignette.
#' @param x is the data set containing the fluorescence amplitudes.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @references M. Febrero-Bande, M.O. de la Fuente, others, \emph{Statistical
#' computing in functional data analysis: The R package fda.usc}, Journal of
#' Statistical Software. 51 (2012) 1--28. http://www.jstatsoft.org/v51/i04/
#'
#' A.-N. Spiess, C. Deutschmann, M. Burdukiewicz, R. Himmelreich, K. Klat, P.
#' Schierack, S. Roediger, Impact of Smoothing on Parameter Estimation in
#' Quantitative DNA Amplification Experiments, Clinical Chemistry. 61 (2015)
#' 379--388. doi:10.1373/clinchem.2014.230656.
#'
#' S. Roediger, A. Boehm, I. Schimke, Surface Melting Curve Analysis with R,
#' \emph{The R Journal}. 5 (2013) 37--53.
#' http://journal.r-project.org/archive/2013-2/roediger-bohm-schimke.pdf.
#'
#' S. Roediger, M. Burdukiewicz, K.A. Blagodatskikh, P. Schierack, R as an
#' Environment for the Reproducible Analysis of DNA Amplification Experiments,
#' \emph{The R Journal}. 7 (2015) 127--150.
#' http://journal.r-project.org/archive/2015-1/RJ-2015-1.pdf.
#'
#' S. Pabinger, S. Roediger, A. Kriegner, K. Vierlinger, A. Weinhauusel, A
#' survey of tools for the analysis of quantitative PCR (qPCR) data, \emph{Biomolecular
#' Detection and Quantification}. 1 (2014) 23--33. doi:10.1016/j.bdq.2014.08.002.
#'
#' S. Roediger, M. Burdukiewicz, P. Schierack, \emph{chipPCR: an R package to
#' pre-process raw data of amplification curves}, \emph{Bioinformatics}. 31 (2015)
#' 2900--2902. doi:10.1093/bioinformatics/btv205.
#'
#' @return Output Description
#' \tabular{llr}{
#'   "eff" \tab qPCR amplification efficiency \tab numeric \cr
#'   "cpD1" \tab maximum of the first derivative curve \tab numeric \cr
#'   "cpD2" \tab maximum of the second derivative curve \tab numeric \cr
#'   "cpD2_approx" \tab maximum of the second derivative curve calculated by the approximate derivative \tab numeric \cr
#'   "cpD2_ratio" \tab a value calculated from the ratio between cpD2 and cpD2_approx \tab numeric \cr
#'   "fluo" \tab raw fluorescence value at the point defined by cpD2 \tab  numeric \cr
#'   "init2" \tab initial template fluorescence from an exponential model \tab numeric \cr
#'   "top" \tab takeoff point. When no top can be determined, the tob value is set to the first cycle number. \tab numeric \cr
#'   "f.top" \tab fluorescence at takeoff point. When no f.tdp can be determined, the f.tdp value is set to the RFU value at the first cycle number. \tab  numeric \cr
#'   "tdp" \tab takes the maximum x fluorescence subtracted by reverse values of the fluorescence and calculates then the fake takeoff point. It is so to speak the take down point (tdp). When no tdp can be determined, the tdb value is set to the last cycle number. \tab numeric \cr
#'   "f.tdp" \tab fluorescence at tdp point. When no f.tdp can be determined, the f.tdp value is set to the RFU value at the last cycle number. \tab  numeric \cr
#'   "sliwin" \tab PCR efficiency by the 'window-of-linearity' method \tab numeric \cr
#'   "b_slope" \tab Is the slope of the seven parameter model \tab  numeric \cr
#'   "b_model_param" \tab Is the b model parameter of the five parameter model \tab  numeric \cr
#'   "c_model_param" \tab Is the c model parameter of the five parameter model \tab  numeric \cr
#'   "d_model_param" \tab Is the d model parameter of the five parameter model \tab  numeric \cr
#'   "e_model_param" \tab Is the e model parameter of the five parameter model \tab  numeric \cr
#'   "f_model_param" \tab Is the f model parameter of the five parameter model \tab  numeric \cr
#'   "f_intercept" \tab Is the intercept of the seven parameter model \tab  numeric \cr
#'   "k1_model_param" \tab Is the k1 model parameter of the seven parameter model \tab  numeric \cr
#'   "k2_model_param" \tab Is the k2 model parameter of the seven parameter model \tab  numeric \cr
#'   "convInfo_iteratons" \tab Number of iterations needed to fit the model \tab numeric \cr
#'   "cpDdiff" \tab absolute difference between cpD1 and cpD2 \tab numeric \cr
#'   "slope_bg" \tab slope of the first cycles \tab numeric \cr
#'   "intercept_bg" \tab intercept of the first cycles \tab numeric \cr
#'   "polyarea" \tab area of a polygon given by the vertices in the vectors cycles and fluorescence \tab numeric \cr
#'   "cp_e.agglo" \tab agglomerative hierarchical estimate for multiple change points \tab numeric \cr
#'   "cp_bcp" \tab change point by Bayesian analysis methods \tab numeric \cr
#'   "qPCRmodel" \tab non-linear model determined for the analysis \tab factor \cr
#'   "amptester_shapiro" \tab tests based on the Shapiro-Wilk normality test if the amplification curve is just noise \tab binary \cr
#'   "amptester_lrt" \tab performs a cycle dependent linear regression and determines if the coefficients of determination deviates from a threshold \tab binary \cr
#'   "amptester_rgt" \tab Resids growth test (RGt) tests if fluorescence values in a linear phase are stable \tab binary \cr
#'   "amptester_tht" \tab  Threshold test (THt) takes the first 20 percent and the last 15 percent of any input data set and performs a Wilcoxon rank sum tests. \tab binary \cr
#'   "amptester_slt" \tab Signal level test compares 1. the signals by a robust "sigma" rule by median + 2 * mad and 2. by comparison of the signal/noise ratio \tab binary \cr
#'   "amptester_polygon" \tab pco test (pco) determines if the points in an amplification curve (like a polygon, in particular non-convex polygons) are in a "clockwise" order. \tab binary \cr
#'   "amptester_slope.ratio" \tab SlR uses the inder function to find the approximated first derivative maximum, second derivative minimum and the second derivative maximum. These are used for a regression analysis with the corresponding fluorescence amplitude data. \tab numeric \cr
#'   "minRFU" \tab minimum of fluorescence amplitude (percentile 0.01) \tab numeric \cr
#'   "maxRFU" \tab maximum of fluorescence amplitude (percentile 0.99) \tab numeric \cr
#'   "bg.start_norm" \tab takes the start (cycle) the amplification curve background based on the bg.max function and normalizes it to the total cycle number \tab numeric \cr
#'   "bg.stop" \tab estimates the end (cycle) the amplification curve background based on the bg.max function and normalizes it to the total cycle number \tab numeric \cr
#'   "amp.stop" \tab estimates the end (cycle) of the amplification curve based in the bg.max function and normalizes it to the total cycle number \tab numeric \cr
#'   "head2tail_ratio" \tab \tab numeric \cr
#'   "autocorrelation" \tab  \tab numeric \cr
#'   "mblrr_intercept_bg" \tab  \tab numeric \cr
#'   "mblrr_slope_bg" \tab \tab numeric \cr
#'   "mblrr_cor_bg" \tab \tab numeric \cr
#'   "mblrr_intercept_pt" \tab \tab numeric \cr
#'   "mblrr_slope_pt" \tab \tab numeric \cr
#'   "mblrr_cor_pt" \tab \tab numeric \cr
#'   "amp_cor_MIC" \tab \tab numeric \cr
#'   "hookreg_hook" \tab estimate of hook effect like curvature \tab binary \cr
#'   "hookreg_hook_slope" \tab estimate of slope of the hook effect like curvature \tab numeric \cr
#'   "peaks_ratio" \tab Takes the estimate approximate local minimums and maximums \tab \cr
#'   "loglin_slope" \tab slope determined by a linear model of the data points from the minimum and maximum of the second derivative \tab numeric \cr
#'   "cpD2_range" \tab cycle difference between the maximum and the minimum of the second derivative curve \tab numeric \cr
#'   "sd_bg" \tab shows the standard deviation of the fluorescence in the ground phase \tab numeric \cr
#'   "central_angle" \tab shows the central angle calculated from the maxima and minima of the derivatives \tab numeric \cr
#' }
#' @details Details can be found in the vignette.
#' @importFrom qpcR pcrfit
#' @examples
#' # Load the chipPCR package and analyze from the C126EG685 the first qPCR run
#' # "A01" (column 2).
#' library(chipPCR)
#' res <- pcrfit_single(C126EG685[, 2])
#' @seealso
#'  \code{\link[bcp]{bcp}}
#'  \code{\link[chipPCR]{bg.max}},\code{\link[chipPCR]{amptester}},\code{\link[chipPCR]{smoother}}
#'  \code{\link[ecp]{e.agglo}}
#'  \code{\link[MBmca]{diffQ}},\code{\link[MBmca]{mcaPeaks}},\code{\link[MBmca]{diffQ2}}
#'  \code{\link{head2tailratio}},\code{\link{earlyreg}},\code{\link{hookreg}},\code{\link{hookregNL}},\code{\link{mblrr}},\code{\link{autocorrelation_test}}
#'  \code{\link[pracma]{polyarea}}
#'  \code{\link[qpcR]{pcrfit}},\code{\link[qpcR]{takeoff}},\code{\link[qpcR]{sliwin}},\code{\link[qpcR]{efficiency}}
#'  \code{\link[base]{diff}}
#'  \code{\link[stats]{quantile}}
#'
#' @rdname pcrfit_single
#' @export pcrfit_single

pcrfit_single <- function(x) {

  # Normalize RFU values to the alpha percentile (0.99)
  x <- x / quantile(x, 0.99, na.rm = TRUE)
  length_cycle <- length(x)
  cycles <- 1L:length_cycle
  
  # Smooth data with moving average smoothing filter for further data
  # analysis steps.
  dat_smoothed <- chipPCR::smoother(cycles, x, method = "mova")
  
  guess_direction <- ifelse(median(head(dat_smoothed, 5)) > (median(tail(dat_smoothed, 5)) + mad(tail(dat_smoothed, 5))), 
           'max', 'min')
  
  # Calculate the first derivative
  res_diffQ <- try(suppressMessages(MBmca::diffQ(cbind(cycles[-c(1:3)], dat_smoothed[-c(1:3)]), verbose = TRUE)$xy), silent = TRUE)
  
  # Determine highest and lowest amplification curve values
  fluo_range <- stats::quantile(x, c(0.01, 0.99), na.rm = TRUE)

  # for qpcR
  dat <- cbind(cyc = cycles, fluo = x)
  dat_reverse <- cbind(cyc = cycles, fluo = (max(dat_smoothed) - rev(dat_smoothed)))
  # inefficient kludge to find l4 model
  data("sysdata", package = "qpcR", envir = parent.frame())

  res_bg.max_tmp <- try(chipPCR::bg.max(cycles, x), silent = TRUE)
  if (inherits_error(res_bg.max_tmp)) {
    res_bg.max <- c(
      bg.start = 0,
      bg.stop = 0,
      amp.stop = 0
    )
  } else {
    res_bg.max <- c(
      bg.start = res_bg.max_tmp@bg.start,
      bg.stop = res_bg.max_tmp@bg.stop,
      amp.stop = res_bg.max_tmp@amp.stop
    )
  }

  # Determine the head to tail ratio
  res_head2tail_ratio <- PCRedux::head2tailratio(x)
  if(res_head2tail_ratio > 50) {res_head2tail_ratio <- 0}

  # Determine the slope from the cycles 2 to 11
  res_lm_fit <- PCRedux::earlyreg(x = cycles, x)

  # Try to estimate the slope and the intercept of the tail region,
  # which might be indicative of a hook effect (strong negative slope)
  res_hookreg_simple <- PCRedux::hookreg(x = cycles, x)
  res_hookregNL <- suppressMessages(PCRedux::hookregNL(x = cycles, x))

  res_hookreg <- ifelse(res_hookreg_simple["hook"] == 1 || res_hookregNL["hook"] == 1, 1, 0)

  # Calculates the area of the amplification curve
  res_polyarea <- try(pracma::polyarea(cycles, x), silent = TRUE) / length_cycle
  if (inherits_error(res_polyarea)) {
    res_polyarea <- 0
  }

  # Calculate change points
  res_cp_e.agglo <- try(length(ecp::e.agglo(as.matrix(res_diffQ[["d(F) / dT"]]))$estimates), silent = TRUE)
  if (inherits_error(res_cp_e.agglo)) {
    res_cp_e.agglo <- length_cycle
  }

  # Bayesian analysis of change points
  res_bcp_tmp <- try(bcp::bcp(res_diffQ[["d(F) / dT"]]), silent = TRUE)
  res_bcp_tmp <- try(res_bcp_tmp$posterior.prob >= 0.6, silent = TRUE)
  res_bcp <- try(length(which(as.factor(res_bcp_tmp) == TRUE)))
  if (inherits_error(res_bcp)) {
    res_bcp <- 0
  }

  # Median based local robust regression (mblrr)
  res_mblrr <- PCRedux::mblrr(cycles, x)


  # Calculate amptester results
  res_amptester <- suppressMessages(try(chipPCR::amptester(x)))

  # Estimate the spread of the approximate local minima and maxima of the curve data

  res_mcaPeaks <- try(MBmca::mcaPeaks(res_diffQ[, 1], res_diffQ[, 2]), silent = TRUE)
  res_peaks_ratio <- try(base::diff(range(res_mcaPeaks$p.max[, 2])) / base::diff(range(res_mcaPeaks$p.min[, 2])), silent = TRUE)
  if (is.infinite(res_peaks_ratio) || inherits_error(res_peaks_ratio) || is.na(res_peaks_ratio)) {
    res_peaks_ratio <- 0
  }

  # Estimate the slope between the minimum and the maximum of the second derivative
  res_diffQ2 <- suppressMessages(MBmca::diffQ2(cbind(cycles[-c(1:3)], dat_smoothed[-c(1:3)]), 
                                               verbose = FALSE, 
                                               fct = get(guess_direction, pos = "package:base"), 
                                               inder = TRUE))
  
  # difference between the minimum and the maximum of the approximate second derivative.
  cpD2_range <- diff(res_diffQ2[[3]])
  if(cpD2_range > 200) cpD2_range <- length_cycle
  if (res_diffQ2[[3]][1] < res_diffQ2[1] && res_diffQ2[1] < res_diffQ2[[3]][2]) {
    cpD2_range
  } else {
    cpD2_range <- length_cycle
  }
  if (res_diffQ2[[3]][1] < res_diffQ2[1] && res_diffQ2[1] < res_diffQ2[[3]][2] && cpD2_range > 1 && cpD2_range < 9) {
    ROI <- round(c(res_diffQ2[[1]], res_diffQ2[[3]]))
    res_loglin_slope <- try(coefficients(lm(x[ROI] ~ ROI))[["ROI"]], silent = FALSE)
    if(res_loglin_slope > 0.45) {res_loglin_slope <- 0}
  } else {
    res_loglin_slope <- 0
  }

  # Perform an autocorrelation analysis
  res_autocorrelation <- PCRedux::autocorrelation_test(y = dat_smoothed)
  if(res_autocorrelation == "n.s." || is.na(res_autocorrelation)) res_autocorrelation <- 0
  
  # Fitting of non-linear multiparameter models
  # Fit sigmoidal models to curve data
  
  b_val <- -10000
  c_val <- -10000
  d_val <- 10000
  e_val <- 10000
  f_val <- -10000
  k_val <- -10000
  k1_val <- -10000
  k2_val <- -10000
  
  # Get the parameters from the four-parameter model
  pcrfit_model_l4 <- try(qpcR::pcrfit(dat, 1, 2, model = l4), silent = TRUE)
  res_coef_model_l4 <- try(coefficients(pcrfit_model_l4), silent = TRUE)
  if (inherits_error(res_coef_model_l4)) {
    res_coef_model_l4 <- c(b = b_val, c = c_val, d = d_val, e = e_val)
  }
  
  # Get the parameters from the five-parameter model
  pcrfit_model_l5 <- try(qpcR::pcrfit(dat, 1, 2, model = l5), silent = TRUE)
  res_coef_model_l5 <- try(coefficients(pcrfit_model_l5), silent = TRUE)
  if (inherits_error(res_coef_model_l5)) {
    res_coef_model_l5 <- c(b = b_val, c = c_val, d = d_val, e = e_val, f = f_val)
  }
  
  # Get the parameters from the five-parameter model
  pcrfit_model_l6 <- try(qpcR::pcrfit(dat, 1, 2, model = l6), silent = TRUE)
  res_coef_model_l6 <- try(coefficients(pcrfit_model_l6), silent = TRUE)
  if (inherits_error(res_coef_model_l6)) {
    res_coef_model_l6 <- c(b = b_val, c = c_val, d = d_val, e = e_val, f = f_val, k = k_val)
  }
  
  # Get the parameters from the seven-parameter model
  pcrfit_model_l7 <- try(qpcR::pcrfit(dat, 1, 2, model = l7), silent = TRUE)
  res_coef <- try(coefficients(pcrfit_model_l7), silent = TRUE)
  if (inherits_error(res_coef)) {
    res_coef <- c(b = b_val, c = c_val, d = d_val, e = e_val, f = f_val, k1 = k1_val, k2 = k2_val)
  }

  res_convInfo_iteratons <- try(pcrfit_model_l7[["convInfo"]][["finIter"]], silent = TRUE)
  if (inherits_error(res_convInfo_iteratons)) {
    res_convInfo_iteratons <- 5000
  }

  # Determine the optimal model based on the AICc

  res_AICc <- list(l4 = if(inherits_error(pcrfit_model_l4)) {
    NA
  } else {
    try(AICc(pcrfit_model_l4), silent = TRUE)
  }, 
  l5 = if(inherits_error(pcrfit_model_l5)) {
    NA
  } else {
    try(AICc(pcrfit_model_l5), silent = TRUE)
  },
  l6 = if(inherits_error(pcrfit_model_l6)) {
    NA
  } else {
    try(AICc(pcrfit_model_l6), silent = TRUE)
  },
  l7 = if(inherits_error(pcrfit_model_l7)) {
    NA
  } else {
    try(AICc(pcrfit_model_l7), silent = TRUE)
  }
  )
  
  res_fit_model <- as.factor(which(res_AICc == min(unlist(res_AICc), na.rm = TRUE)))

  if (length(res_fit_model) == 0) {
      res_fit_model <- as.factor("l0")
      names(res_fit_model) <- "l0"
  }
  
  if(names(res_fit_model) == "l4") res_fit <- pcrfit_model_l4
  if(names(res_fit_model) == "l5") res_fit <- pcrfit_model_l5
  if(names(res_fit_model) == "l6") res_fit <- pcrfit_model_l6
  if(names(res_fit_model) == "l7") res_fit <- pcrfit_model_l7
  if(names(res_fit_model) == "l0") res_fit <- "try-error"
  
  # Fit the "reverse" model
  pcrfit_startmodel_reverse <- try(qpcR::pcrfit(dat_reverse, 1, 2), silent = TRUE)
  
  res_fit_reverse <- try(qpcR::mselect(
    pcrfit_startmodel_reverse,
    verbose = FALSE, fctList = list(l4, l5, l6, l7)
  ), silent = TRUE)

  # Determine the model for the reverse data suggested by the
  # mselect function based on the AICc
  res_fit_model_reverse <- try(names(which(res_fit_reverse[["retMat"]][, "AICc"] == min(res_fit_reverse[["retMat"]][, "AICc"]))), silent = TRUE)
  if (inherits_error(res_fit_model_reverse)) {
    res_fit_model_reverse <- as.factor("l0")
  }

  if (!inherits_error(res_fit_reverse)) {
    # TakeOff Point from the reverse data
    # Calculates the first significant cycle of the exponential region
    #
    #   Take Down Point tdp
    #
    # (takeoff point) using externally studentized residuals as described
    # in Tichopad et al. (2003).
    res_takeoff_reverse <- try(qpcR::takeoff(res_fit_reverse, nsig = 5), silent = TRUE)
    res_takeoff_reverse[[1]] <- length_cycle - res_takeoff_reverse[[1]]
    res_takeoff_reverse[[2]] <- x[res_takeoff_reverse[[1]]] -
      res_takeoff_reverse[[2]] + min(x)
    if (is.na(res_takeoff_reverse[[1]])) {
      res_takeoff_reverse <- list(1, x[length_cycle])
    }
    # Calculate the standard deviation of the fluorescence starting from
    # cylce 2 to the takeoff point
    if (!is.na(res_takeoff_reverse[[1]])) {
      res_sd_bg <- try(sd(x[2L:res_takeoff_reverse[[1]]]), silent = TRUE)
    } else {
      res_sd_bg <- sd(x[1L:8])
    }
    if (inherits_error(res_takeoff_reverse)) {
      res_takeoff_reverse <- list(0, 1)
    }
  } else {
    res_takeoff_reverse <- list(length_cycle, 1)
    # Calculate the standard deviation of the fluorescence starting from
    # cylce 2 to cycle 8 if the the takeoff point cannot be
    # determined
    res_sd_bg <- sd(x[1L:8])
  }
    names(res_takeoff_reverse) <- c("tdp", "f.tdp")

  if (!inherits_error(res_fit)) {
    # TakeOff Point
    # Calculates the first significant cycle of the exponential region
    # (takeoff point) using externally studentized residuals as described
    # in Tichopad et al. (2003).
    res_takeoff <- try(qpcR::takeoff(res_fit), silent = TRUE)
    if (inherits_error(res_takeoff)) {
      res_takeoff <- list(0, 1)
    }
    if (is.na(res_takeoff[[1]])) {
      res_takeoff <- list(0, 1)
    }
    names(res_takeoff) <- c("top", "f.top")

    # sliwin qPCR efficiency
    # Calculation of the qPCR efficiency by the 'window-of-linearity' method
    res_sliwin <- try(
      qpcR::sliwin(res_fit, plot = FALSE, verbose = FALSE)$eff,
      silent = TRUE
    )
    if (inherits_error(res_sliwin)) {
      res_sliwin <- 0
    }

    # Cq of the amplification curve
    # Determine the Cq and other parameters
    res_efficiency_tmp <- try(
      qpcR::efficiency(res_fit, plot = FALSE)[c(
        "eff",
        "cpD1", "cpD2",
        "fluo",
        "init2"
      )],
      silent = TRUE
    )
    
    
    if (!inherits_error(res_efficiency_tmp)) {
      res_cpDdiff <- try(abs(res_efficiency_tmp[["cpD1"]] - res_efficiency_tmp[["cpD2"]]))
      
      if(res_efficiency_tmp$init2 < -200 || as.character(res_efficiency_tmp$init2) == "NaN") {res_efficiency_tmp$init2 <- 0}
      if(res_efficiency_tmp$eff > 2.2 || as.character(res_efficiency_tmp$eff) == "NaN") {res_efficiency_tmp$eff <- 0}
      
    } else {
      res_efficiency_tmp <- list(
        eff = 0,
        cpD1 = 0,
        cpD2 = 0,
        fluo = 1,
        init2 = 1
      )
      res_cpDdiff <- 1
    }
  } else {
    res_efficiency_tmp <- list(
      eff = 0,
      cpD1 = 0,
      cpD2 = 0,
      fluo = 1,
      init2 = 1
    )
    res_takeoff <- list(0, x[1])
    res_takeoff_reverse <- list(0, 1)
    res_sliwin <- 0
    res_cpDdiff <- 0
  }
  
  # Calculate the ratio between the the approximate second derivative maximum 
  # cpD2_approx and the second derivative maximum cpD2
  
  cpD2_ratio_tmp <- res_efficiency_tmp[["cpD2"]] / res_diffQ2[[3]][1]
  if(cpD2_ratio_tmp == Inf) cpD2_ratio <- 0  
  cpD2_ratio <- cpD2_ratio_tmp

  
  # The cpD2_ratio is a binary value, which is based on a range empirically
  # determined on the data_sample data set.
  # 
  # if(cpD2_ratio_tmp > 0.8 && cpD2_ratio_tmp < 1.1) {
  #     cpD2_ratio <- 1
  # } else {
  #     cpD2_ratio <- 0
  # }

  
# # Calculate angle
#     # Get x coordinates for the vectors be the inder 
#     # function. The first derivative maximum (FDM) is the
#     # center. The second derivative maximum (SDM) is the left
#     # and second derivative minimum (SDm) is the right point of
#     # the vectors.

    origin <- data.frame(res_diffQ2[["Tm"]][1], res_diffQ2[["fluoTm"]][1])
    point_x1 <- data.frame(res_diffQ2[["xTm1.2.D2"]][1], res_diffQ2[["yTm1.2.D2"]][1])
    point_x2 <- data.frame(res_diffQ2[["xTm1.2.D2"]][2], res_diffQ2[["yTm1.2.D2"]][2])

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
    
    # Calculate central angle
    res_angle <- dot_product / length_of_vectors

  all_results <- data.frame(
    # Quantification points, derivatives, efficiencies,
    # curve fitting
    cpD1 = res_efficiency_tmp[["cpD1"]],
    cpD2 = res_efficiency_tmp[["cpD2"]],
    cpD2_approx = res_diffQ2[[3]][1],
    cpD2_ratio = cpD2_ratio,
    eff = res_efficiency_tmp[["eff"]],
    sliwin = res_sliwin[[1]],
    cpDdiff = res_cpDdiff,
    loglin_slope = res_loglin_slope,
    cpD2_range = cpD2_range,
    top = res_takeoff[[1]],
    f.top = res_takeoff[[2]],
    tdp = res_takeoff_reverse[[1]],
    f.tdp = res_takeoff_reverse[[2]],
    bg.stop = res_bg.max[2],
    amp.stop = res_bg.max[3],
    b_slope = res_coef[["b"]],
    b_model_param = res_coef_model_l5[["b"]],
    c_model_param = res_coef_model_l5[["c"]],
    d_model_param = res_coef_model_l5[["d"]],
    e_model_param = res_coef_model_l5[["e"]],
    f_model_param = res_coef_model_l5[["f"]],
    f_intercept = res_coef[["f"]],
#    k1_model_param = res_coef[["k1"]],
#    k2_model_param = res_coef[["k2"]],
    convInfo_iteratons = res_convInfo_iteratons,
    qPCRmodel = factor(names(res_fit_model)), # res_fit_model[[1]],
    qPCRmodelRF = res_fit_model_reverse[[1]],
    # Signal levels
    minRFU = fluo_range[[1]],
    maxRFU = fluo_range[[2]],
    init2 = res_efficiency_tmp[["init2"]],
    fluo = res_efficiency_tmp[["fluo"]],
    slope_bg = res_lm_fit[["slope"]],
    intercept_bg = res_lm_fit[["intercept"]],
    sigma_bg = res_lm_fit[["sigma"]],
    sd_bg = res_sd_bg,
    head2tail_ratio = res_head2tail_ratio,
    mblrr_slope_pt = res_mblrr[5],
    mblrr_intercept_bg = res_mblrr[1],
    mblrr_slope_bg = res_mblrr[2],
    mblrr_cor_bg = res_mblrr[3],
    mblrr_intercept_pt = res_mblrr[4],
    mblrr_cor_pt = res_mblrr[6],
    # Areas
    polyarea = res_polyarea,
    peaks_ratio = res_peaks_ratio,
    autocorrelation = res_autocorrelation,
    # Change points
    cp_e.agglo = res_cp_e.agglo,
    cp_bcp = res_bcp,
    # Amptester
    amptester_shapiro = res_amptester@decisions["shap.noisy"][[1]],
    amptester_lrt = res_amptester@decisions["lrt.test"][[1]],
    amptester_rgt = res_amptester@decisions["rgt.dec"][[1]],
    amptester_tht = res_amptester@decisions["tht.dec"][[1]],
    amptester_slt = res_amptester@decisions["slt.dec"][[1]],
    amptester_polygon = res_amptester@"polygon" / length_cycle,
    amptester_slope.ratio = ifelse(is.na(res_amptester@"slope.ratio"), 0, res_amptester@"slope.ratio"),
    # Curvature
    hookreg_hook = res_hookreg,
    hookreg_hook_slope = res_hookreg_simple[["slope"]],
    hookreg_hook_delta = res_hookreg_simple[["hook.delta"]],
    # Angle
    central_angle = res_angle,
    # Number of Cycles
    number_of_cycles = length_cycle,
    # Identifier
    row.names = "results"
  )
}
