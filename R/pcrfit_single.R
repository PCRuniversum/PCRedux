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
#'   "fluo" \tab raw fluorescence value at the point defined by cpD2 \tab  numeric \cr
#'   "init2" \tab initial template fluorescence from an exponential model \tab numeric \cr
#'   "top" \tab takeoff point. When no top can be determined, the tob value is set to the first cycle number. \tab numeric \cr
#'   "f.top" \tab fluorescence at takeoff point. When no f.tdp can be determined, the f.tdp value is set to the RFU value at the first cycle number. \tab  numeric \cr
#'   "tdp" \tab takes the maximum x fluorescence subtracted by reverse values of the fluorescence and calculates then the fake takeoff point. It is so to speak the take down point (tdp). When no tdp can be determined, the tdb value is set to the last cycle number. \tab numeric \cr
#'   "f.tdp" \tab fluorescence at tdp point. When no f.tdp can be determined, the f.tdp value is set to the RFU value at the last cycle number. \tab  numeric \cr
#'   "ressliwin" \tab PCR efficiency by the 'window-of-linearity' method \tab numeric \cr
#'   "b_slope" \tab Is the slope of seven parameter model \cr
#'   "f_intercept" \tab Is the intercept of the seven parameter model \cr
#'   "convInfo_iteratons" \tab Number of iterations needed to fit the model \tab numeric \cr
#'   "cpDdiff" \tab difference between cpD1 and cpD2 \tab numeric \cr
#'   "slope_bg" \tab slope of the first cycles \tab numeric \cr
#'   "intercept_bg" \tab intercept of the first cycles \tab numeric \cr
#'   "polyarea" \tab area of a polygon given by the vertices in the vectors cycles and fluorescence \tab numeric \cr
#'   "changepoint_e.agglo" \tab agglomerative hierarchical estimate for multiple change points \tab numeric \cr
#'   "changepoint_bcp" \tab change point by Bayesian analysis methods \tab numeric \cr
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
#'   "bg.stop_norm" \tab estimates the end (cycle) the amplification curve background based on the bg.max function and normalizes it to the total cycle number \tab numeric \cr
#'   "amp.stop_norm" \tab estimates the end (cycle) of the amplification curve based in the bg.max function and normalizes it to the total cycle number \tab numeric \cr
#'   "head_to_tail_ratio" \tab \tab numeric \cr
#'   "autocorellation" \tab  \tab numeric \cr
#'   "mblrr_intercept_bg" \tab  \tab numeric \cr
#'   "mblrr_slope_bg" \tab \tab numeric \cr
#'   "mblrr_cor_bg" \tab \tab numeric \cr
#'   "mblrr_intercept_pt" \tab \tab numeric \cr
#'   "mblrr_slope_pt" \tab \tab numeric \cr
#'   "mblrr_cor_pt" \tab \tab numeric \cr
#'   "amp_cor_MIC" \tab \tab numeric \cr
#'   "hookreg_hook" \tab estimate of hook effect like curvature \tab binary \cr
#'   "hookreg_hook_slope" \tab estimate of slope of the hook effect like curvature \tab numeric \cr
#'   "peaks_min_max_ratio" \tab Takes the estimate approximate local minimums and maximums \tab \cr
#'   "diffQ2_slope" \tab slope determined by a linear model of the data points from the minimum and maximum of the second derivative \tab numeric \cr
#'   "diffQ2_Cq_range" \tab cycle difference between the maximum and the minimum of the second derivative curve \tab numeric \cr
#'   "sd_bg \tab shows the standard deviation of the fluorescence in the ground phase \tab numeric \cr
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
  
  # Smooth data with moving average for other data
  # analysis steps.
  dat_smoothed <- chipPCR::smoother(cycles, x)
  
  # Calculate the first derivative
  res_diffQ <- suppressMessages(MBmca::diffQ(cbind(cycles, dat_smoothed), verbose = TRUE)$xy)
  
  # Determine highest and lowest amplification curve values
  fluo_range <- stats::quantile(x, c(0.01, 0.99), na.rm = TRUE)

  # for qpcR
  dat <- cbind(cyc = cycles, fluo = x)
  dat_reverse <- cbind(cyc = cycles, fluo = (max(x) - rev(x)))
  # inefficient kludge to find l4 model
  data("sysdata", package = "qpcR", envir = parent.frame())

  res_bg.max_tmp <- try(chipPCR::bg.max(cycles, x), silent = TRUE)
  if (class(res_bg.max_tmp) == "try-error") {
    res_bg.max <- c(
      bg.start = 0,
      bg.stop = 0,
      amp.stop = 0
    )
  } else {
    res_bg.max <- c(
      bg.start = res_bg.max_tmp@bg.start / length_cycle,
      bg.stop = res_bg.max_tmp@bg.stop / length_cycle,
      amp.stop = res_bg.max_tmp@amp.stop / length_cycle
    )
  }

  # Determine the head to tail ratio
  res_head_tail_ratio <- PCRedux::head2tailratio(x)

  # Determine the slope from the cycles 2 to 11
  res_lm_fit <- PCRedux::earlyreg(x = cycles, x)

  # Try to estimate the slope and the intercept of the tail region,
  # which might be indicative of a hook effect (strong negative slope)
  res_hookreg_simple <- PCRedux::hookreg(x = cycles, x)
  res_hookregNL <- suppressMessages(PCRedux::hookregNL(x = cycles, x))

  res_hookreg <- ifelse(res_hookreg_simple["hook"] == 1 || res_hookregNL["hook"] == 1, 1, 0)

  # Calculates the area of the amplification curve
  res_polyarea <- try(pracma::polyarea(cycles, x), silent = TRUE) / length_cycle
  if (class(res_polyarea) == "try-error") {
#     res_polyarea <- NA
    res_polyarea <- 0
  }

  # Calculate change points
  res_changepoint_e.agglo <- try(length(ecp::e.agglo(as.matrix(x))$estimates), silent = TRUE)
  res_changepoint_e.agglo <- try(length(ecp::e.agglo(as.matrix(res_diffQ[["d(F) / dT"]]))$estimates), silent = TRUE)
  if (class(res_changepoint_e.agglo) == "try-error") {
    res_changepoint_e.agglo <- length_cycle
  }

  # Bayesian analysis of change points
  res_bcp_tmp <- bcp::bcp(res_diffQ[["d(F) / dT"]])
  res_bcp_tmp <- res_bcp_tmp$posterior.prob >= 0.6
  res_bcp <- try((which(as.factor(res_bcp_tmp) == TRUE) %>% length()))
  if (class(res_bcp) == "try-error") {
    res_bcp <- 0
  }

  # Median based local robust regression (mblrr)
  res_mblrr <- PCRedux::mblrr(cycles, x)


  # Calculate amptester results
  res_amptester <- suppressMessages(try(chipPCR::amptester(x)))

  # Estimate the spread of the approximate local minima and maxima of the curve data

  res_mcaPeaks <- MBmca::mcaPeaks(res_diffQ[, 1], res_diffQ[, 2])
  peaks_min_max_ratio <- base::diff(range(res_mcaPeaks$p.max[, 2])) / base::diff(range(res_mcaPeaks$p.min[, 2]))
  if (is.infinite(peaks_min_max_ratio)) {
    peaks_min_max_ratio <- NA
  }

  # Estimate the slope between the minimum and the maximum of the second derivative
  res_diffQ2 <- suppressMessages(MBmca::diffQ2(cbind(cycles, dat_smoothed), verbose = FALSE, fct = min))
  range_Cq <- diff(res_diffQ2[[3]])
  if (res_diffQ2[[3]][1] < res_diffQ2[1] && res_diffQ2[1] < res_diffQ2[[3]][2]) {
    range_Cq
  } else {
    range_Cq <- 0
  }
  if (res_diffQ2[[3]][1] < res_diffQ2[1] && res_diffQ2[1] < res_diffQ2[[3]][2] && range_Cq > 1 && range_Cq < 9) {
    res_diffQ2_slope <- coefficients(lm(unlist(c(res_diffQ2[[4]])) ~ unlist(c(res_diffQ2[[3]]))))[2]
  } else {
    res_diffQ2_slope <- 0
  }

  # Perform an autocorrelation analysis
  res_autocorrelation <- PCRedux::autocorrelation_test(y = x)

  # Fit sigmoidal models to curve data

  pcrfit_startmodel <- try(qpcR::pcrfit(dat, 1, 2, model = l7), silent = TRUE)

  res_coef <- try(coefficients(pcrfit_startmodel), silent = TRUE)
  if (class(res_coef) == "try-error") {
    res_coef <- c(b = 0, f = 0)
  }


  res_convInfo_iteratons <- try(pcrfit_startmodel[["convInfo"]][["finIter"]], silent = TRUE)
  if (class(res_convInfo_iteratons) == "try-error") {
    res_convInfo_iteratons <- 5000
  }

  pcrfit_startmodel_reverse <- try(qpcR::pcrfit(dat_reverse, 1, 2), silent = TRUE)

  res_fit <- try(qpcR::mselect(
    pcrfit_startmodel,
    verbose = FALSE, fctList = list(l4, l5, l6, l7)
  ), silent = TRUE)

  res_fit_reverse <- try(qpcR::mselect(
    pcrfit_startmodel_reverse,
    verbose = FALSE, fctList = list(l4, l5, l6, l7)
  ), silent = TRUE)

  # Determine the model suggested by the mselect function based on the AICc
  res_fit_model <- try(names(which(res_fit[["retMat"]][, "AICc"] == min(res_fit[["retMat"]][, "AICc"]))), silent = TRUE)
  if (class(res_fit_model) == "try-error") {
    res_fit_model <- as.factor("l0")
  }

  # Determine the model for the reverse data suggested by the
  # mselect function based on the AICc
  res_fit_model_reverse <- try(names(which(res_fit_reverse[["retMat"]][, "AICc"] == min(res_fit_reverse[["retMat"]][, "AICc"]))), silent = TRUE)
  if (class(res_fit_model_reverse) == "try-error") {
    res_fit_model_reverse <- as.factor("l0")
  }

  if (class(res_fit_reverse)[1] != "try-error") {
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
      res_takeoff_reverse <- list(length_cycle, x[length_cycle])
    }
    # Calculate the standard deviation of the fluorescence starting from
    # cylce 2 to the takeoff point
    if (!is.na(res_takeoff_reverse[[1]])) {
      sd_bg <- try(sd(x[2L:res_takeoff_reverse[[1]]]), silent = TRUE)
    } else {
      sd_bg <- sd(x[2L:8])
    }
    if (class(res_takeoff_reverse) == "try-error") {
      res_takeoff_reverse <- list(length_cycle, 1)
    }
  } else {
    res_takeoff_reverse <- list(length_cycle, 1)
    # Calculate the standard deviation of the fluorescence starting from
    # cylce 2 to cycle 8 if the the takeoff point cannot be
    # determined
    sd_bg <- sd(x[2L:8])
  }
    names(res_takeoff_reverse) <- c("tdp", "f.tdp")

  if (class(res_fit)[1] != "try-error") {
    # TakeOff Point
    # Calculates the first significant cycle of the exponential region
    # (takeoff point) using externally studentized residuals as described
    # in Tichopad et al. (2003).
    res_takeoff <- try(qpcR::takeoff(res_fit), silent = TRUE)
    if (class(res_takeoff) == "try-error") {
      res_takeoff <- list(length_cycle, 1)
    }
    if (is.na(res_takeoff[[1]])) {
      res_takeoff <- list(length_cycle, 1)
    }
    names(res_takeoff) <- c("top", "f.top")

    # sliwin qPCR efficiency
    # Calculation of the qPCR efficiency by the 'window-of-linearity' method
    res_sliwin <- try(
      qpcR::sliwin(res_fit, plot = FALSE, verbose = FALSE)$eff,
      silent = TRUE
    )
    if (class(res_sliwin) == "try-error") {
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
    if (class(res_efficiency_tmp) != "try-error") {
      res_cpDdiff <- try(abs(res_efficiency_tmp[["cpD1"]] - res_efficiency_tmp[["cpD2"]]))
    } else {
      res_efficiency_tmp <- list(
        eff = 0,
        cpD1 = 0,
        cpD2 = 0,
        fluo = 1,
        init2 = 1
      )
      res_cpDdiff <- length_cycle
    }
  } else {
    res_efficiency_tmp <- list(
      eff = 0,
      cpD1 = 0,
      cpD2 = 0,
      fluo = 1,
      init2 = 1
    )
    res_takeoff <- list(1, x[1])
    res_takeoff_reverse <- list(length_cycle, 1)
    res_sliwin <- 0
    res_cpDdiff <- length_cycle
  }

  all_results <- data.frame(
    # Quantification points, derivatives, efficiencies,
    # curve fitting
    cpD1 = res_efficiency_tmp[["cpD1"]],
    cpD2 = res_efficiency_tmp[["cpD2"]],
    eff = res_efficiency_tmp[["eff"]],
    ressliwin = res_sliwin[[1]],
    cpDdiff = res_cpDdiff,
    diffQ2_slope = res_diffQ2_slope,
    diffQ2_Cq_range = range_Cq,
    top = res_takeoff[[1]],
    f.top = res_takeoff[[2]],
    tdp = res_takeoff_reverse[[1]],
    f.tdp = res_takeoff_reverse[[2]],
    bg.stop_norm = res_bg.max[2],
    amp.stop_norm = res_bg.max[3],
    b_slope = res_coef[["b"]],
    f_intercept = res_coef[["f"]],
    convInfo_iteratons = res_convInfo_iteratons,
    qPCRmodel = res_fit_model[[1]],
    qPCRmodel_reverse = res_fit_model_reverse[[1]],
    # Signal levels
    minRFU = fluo_range[[1]],
    maxRFU = fluo_range[[2]],
    init2 = res_efficiency_tmp[["init2"]],
    fluo = res_efficiency_tmp[["fluo"]],
    slope_bg = res_lm_fit[["slope"]],
    intercept_bg = res_lm_fit[["intercept"]],
    sd_bg = sd_bg,
    head_to_tail_ratio = res_head_tail_ratio,
    mblrr_slope_pt = res_mblrr[5],
    mblrr_intercept_bg = res_mblrr[1],
    mblrr_slope_bg = res_mblrr[2],
    mblrr_cor_bg = res_mblrr[3],
    mblrr_intercept_pt = res_mblrr[4],
    mblrr_cor_pt = res_mblrr[6],
    # Areas
    polyarea = res_polyarea,
    peaks_min_max_ratio = peaks_min_max_ratio,
    autocorellation = res_autocorrelation,
    # Change points
    changepoint_e.agglo = res_changepoint_e.agglo,
    changepoint_bcp = res_bcp,
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
    # Identifier
    row.names = "results"
  )
}
