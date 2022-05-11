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
#'   "cpD1" \tab maximum of the first derivative curve \tab numeric \cr
#'   "cpD2" \tab maximum of the second derivative curve \tab numeric \cr
#'   "cpD2_approx" \tab maximum of the second derivative curve calculated by the approximate derivative \tab numeric \cr
#'   "cpD2_ratio" \tab a value calculated from the ratio between cpD2 and cpD2_approx \tab numeric \cr
#'   "eff" \tab qPCR amplification efficiency \tab numeric \cr
#'   "sliwin" \tab qPCR amplification efficiency according the the 'window-of-linearity' method by Ruijter et al. (2009) \tab numeric \cr
#'   "cpDdiff" \tab absolute difference between cpD1 and cpD2 \tab numeric \cr
#'   "loglin_slope" \tab slope determined by a linear model of the data points from the minimum and maximum of the second derivative \tab numeric \cr
#'   "cpD2_range" \tab cycle difference between the maximum and the minimum of the second derivative curve \tab numeric \cr
#'   "top" \tab takeoff point. When no top can be determined, the tob value is set to the first cycle number. \tab numeric \cr
#'   "f.top" \tab fluorescence at takeoff point. When no f.tdp can be determined, the f.tdp value is set to the RFU value at the first cycle number. \tab  numeric \cr
#'   "tdp" \tab takes the maximum fluorescence subtracted by reverse values of the fluorescence and calculates then the fake takeoff point. It is so to speak the take down point (tdp). When no tdp can be determined, the tdb value is set to the last cycle number. \tab numeric \cr
#'   "f.tdp" \tab fluorescence at tdp point. When no f.tdp can be determined, the f.tdp value is set to the RFU value at the last cycle number. \tab  numeric \cr
#'   "bg.stop" \tab estimates the end (cycle) the amplification curve background based on the bg.max function and normalizes it to the total cycle number \tab numeric \cr
#'   "amp.stop" \tab estimates the end (cycle) of the amplification curve based in the bg.max function and normalizes it to the total cycle number \tab numeric \cr
#'   "b_slope" \tab Is the slope of the seven parameter model \tab  numeric \cr
#'   "b_model_param" \tab Is the b model parameter of the model optimally fitted according to the AIC \tab  numeric \cr
#'   "c_model_param" \tab Is the c model parameter of the model optimally fitted according to the AIC \tab  numeric \cr
#'   "d_model_param" \tab Is the d model parameter of the model optimally fitted according to the AIC \tab  numeric \cr
#'   "e_model_param" \tab Is the e model parameter of the model optimally fitted according to the AIC \tab  numeric \cr
#'   "f_model_param" \tab Is the f model parameter of the model optimally fitted according to the AIC \tab  numeric \cr
#'   "f_intercept" \tab Is the intercept of the seven parameter model \tab  numeric \cr
#'   "convInfo_iteratons" \tab Number of iterations needed to fit the 7 parameter model \tab numeric \cr
#'   "qPCRmodel" \tab non-linear model determined for the analysis \tab factor \cr
#'   "qPCRmodelRF" \tab non-linear model determined for the analysis of the reversed amplification curve \tab factor \cr
#'   "minRFU" \tab minimum of fluorescence amplitude \tab numeric \cr
#'   "maxRFU" \tab maximum of fluorescence amplitude \tab numeric \cr
#'   "init2" \tab initial template fluorescence from an exponential model \tab numeric \cr
#'   "fluo" \tab raw fluorescence value at the point defined by cpD2 \tab  numeric \cr
#'   "slope_bg" \tab slope of the first cycles \tab numeric \cr
#'   "k1_model_param" \tab Is the k1 model parameter of the seven parameter model \tab  numeric \cr
#'   "k2_model_param" \tab Is the k2 model parameter of the seven parameter model \tab  numeric \cr
#'   "intercept_bg" \tab intercept of the first cycles \tab numeric \cr
#'   "sigma_bg" \tab sigma of background \tab numeric \cr
#'   "sd_bg" \tab standard deviation of the background (ground phase) region (start to takeoff point) \tab numeric \cr
#'   "head2tail_ratio" \tab ratio between the signal of the background and tail region \tab numeric \cr
#'   "mblrr_intercept_bg" \tab the value of the intercept in the estimated background region of the amplification curve \tab numeric \cr
#'   "mblrr_slope_bg" \tab the value of the slope in the estimated background region of the amplification curve \tab numeric \cr
#'   "mblrr_cor_bg" \tab the value of the linear correlation coefficient in the estimated background region of the amplification curve \tab numeric \cr
#'   "mblrr_intercept_pt" \tab the value of the intercept in the estimated plateau phase of the amplification curve \tab numeric \cr
#'   "mblrr_slope_pt" \tab the value of the slope in the estimated plateau phase of the amplification curve \tab numeric \cr
#'   "mblrr_cor_pt" \tab the value of the linear correlation coefficient in the estimated plateau phase of the amplification curve \tab numeric \cr
#'   "polyarea" \tab area of a polygon given by the vertices in the vectors cycles and fluorescence \tab numeric \cr
#'   "peaks_ratio" \tab Takes the estimate approximate local minimums and maximums \tab \cr
#'   "autocorrelation" \tab is a value of autocorrelation of a gain curve from a quantitative PCR experiment \tab numeric \cr
#'   "cp_e.agglo" \tab agglomerative hierarchical estimate for multiple change points \tab numeric \cr
#'   "cp_bcp" \tab change point by Bayesian analysis methods \tab numeric \cr
#'   "amptester_shapiro" \tab tests based on the Shapiro-Wilk normality test if the amplification curve is just noise \tab binary \cr
#'   "amptester_lrt" \tab performs a cycle dependent linear regression and determines if the coefficients of determination deviates from a threshold \tab binary \cr
#'   "amptester_rgt" \tab Resids growth test (RGt) tests if fluorescence values in a linear phase are stable \tab binary \cr
#'   "amptester_tht" \tab  Threshold test (THt) takes the first 20 percent and the last 15 percent of any input data set and performs a Wilcoxon rank sum tests. \tab binary \cr
#'   "amptester_slt" \tab Signal level test compares 1. the signals by a robust "sigma" rule by median + 2 * mad and 2. by comparison of the signal/noise ratio \tab binary \cr
#'   "amptester_polygon" \tab pco test (pco) determines if the points in an amplification curve (like a polygon, in particular non-convex polygons) are in a "clockwise" order. \tab binary \cr
#'   "amptester_slope.ratio" \tab SlR uses the inder function to find the approximated first derivative maximum, second derivative minimum and the second derivative maximum. These are used for a regression analysis with the corresponding fluorescence amplitude data. \tab numeric \cr
#'   "hookreg_hook" \tab estimate of hook effect like curvature \tab binary \cr
#'   "hookreg_hook_slope" \tab estimate of slope of the hook effect like curvature \tab numeric \cr
#'   "hookreg_hook_delta" \tab Estimated value for the number of cycles from the qPCR cycle where the hook effect was determined up to the last qPCR cycle \tab numeric \cr
#'   "central_angle" \tab shows the central angle calculated from the maximum and minimum of the second derivatives, with the first derivative maximum being the center \tab numeric \cr
#'   "sd_bg" \tab shows the standard deviation of the fluorescence in the ground phase \tab numeric \cr
#'   "number_of_cycles" \tab Number of cylces \tab numeric \cr
#'   "direction" \tab test if the maximum of the first derivative is positive or negative \tab numeric \cr
#'   "range" \tab outputs the difference of fluorescence between 0.99 and 0.01 percentile. The value thus corresponds approximately to the maximum achievable signal difference of an amplification curve. \tab numeric \cr
#'   "polyarea_trapz" \tab calculates trapezoidal integration. The calculation stops when the difference from one step to the next is smaller than a tolerance value, or the iterations become too large. The value corresponds to the sum signal via the total amplification curve. \tab numeric \cr
#'   "cor" \tab is the value of the correlation coefficient from a linear correlation analysis according to Pearson between all PCR cycles and the fluorescence signals. \tab numeric \cr
#'   "res_coef_pcrfit.b" \tab is the parameter from the adjustment with a nonlinear (sigmoid) four-parametric model which describes the Hill’s slope of the curve (i.e. this is related to the steepness of the curve point e) \tab numeric \cr
#'   "res_coef_pcrfit.c" \tab is the parameter from the adjustment with a nonlinear (sigmoid) four-parametric model which describes the maximum value that can be obtained (i.e. what happens at infinite number of cycles) \tab numeric \cr
#'   "res_coef_pcrfit.d" \tab is the parameter from the adjustment with a nonlinear (sigmoid) four-parametric model which describes the minimum value that can be obtained (i.e. what happens at 0 cycles) \tab numeric \cr
#'   "res_coef_pcrfit.e" \tab is the parameter from the adjustment with a nonlinear (sigmoid) four-parametric model, which describes the point of inflection (i.e. the point on the sigmoid curve halfway between d and c) \tab numeric \cr
#'   "fitAIC" \tab is the value of the Akaike's second-order corrects Information Criterion, which was determined on a non-linear (sigmoid) four-parameter model \tab numeric \cr
#'   "fitIter" \tab Number of iterations needed to fit the 4 parameter model \tab numeric \cr
#'   "segment_x" \tab Adjusts a regression model with segmented (linear) relationships between fluorescence and PCR cycles. This segment describes the baseline in an amplification curve. \tab numeric \cr
#'   "segment_U1.x" \tab Adjusts a regression model with segmented (linear) relationships between fluorescence and PCR cycles. This segment describes the slope in an amplification curve. \tab numeric \cr
#'   "segment_U2.x" \tab Adjusts a regression model with segmented (linear) relationships between fluorescence and PCR cycles. This segment describes the plateau in an amplification curve. \tab numeric \cr
#'   "segment_psi1.x" \tab Adjusts a regression model with segmented (linear) relationships between fluorescence and PCR cycles. The value is based on the break-point(s) fixed at the values and describes the transition from the baseline phase to the exponential phase. \tab numeric \cr
#'   "segment_psi2.x" \tab Adjusts a regression model with segmented (linear) relationships between fluorescence and PCR cycles. The value is based on the break-point(s) fixed at the values and describes the transition from the exponential phase to the plateau phase. \tab numeric \cr
#'   "sumdiff" \tab describes proportion of cycles x in which the fluorescence signal of x is smaller than in x+1 \tab numeric \cr
#'   "poly_1" \tab is a value of a third-order polynomial a + b*x + c*x^2 + d*x^3 is fitted to the curve data, where the intercept correspond to the baseline and the three predictor-dependent terms deliver a approximation to the sigmoidal curve structure. This value describes the intercept "a" (cutting point of the amplification curve model with the ordinate). \tab numeric \cr
#'   "poly_2" \tab is a value of a third-order polynomial a + b*x + c*x^2 + d*x^3 is fitted to the curve data, where the intercept correspond to the baseline and the three predictor-dependent terms deliver a approximation to the sigmoidal curve structure. This value describes the linear part "b*x". \tab numeric \cr
#'   "poly_3" \tab is a value of a third-order polynomial a + b*x + c*x^2 + d*x^3 is fitted to the curve data, where the intercept correspond to the baseline and the three predictor-dependent terms deliver a approximation to the sigmoidal curve structure. This value describes the quadratic part "c*x^2". \tab numeric \cr
#'   "poly_4" \tab is a value of a third-order polynomial a + b*x + c*x^2 + d*x^3 is fitted to the curve data, where the intercept correspond to the baseline and the three predictor-dependent terms deliver a approximation to the sigmoidal curve structure. This value describes the cubic part "d*x^3". \tab numeric \cr
#'   "window_Win_1" \tab The complete curve trajectory is segmented into 10 equidistant windows by fitting an interpolating smoothing spline with smoothing factor 0.5 to the curve, interpolating exactly 50 curve points, and then cutting these into 10 windows of five values each, with a subsequent calculation of the MAD/Median ratio for each of these windows. This is the 1st window for cycles (0.961,4.9]. \tab numeric \cr
#'   "window_Win_2" \tab The complete curve trajectory is segmented into 10 equidistant windows by fitting an interpolating smoothing spline with smoothing factor 0.5 to the curve, interpolating exactly 50 curve points, and then cutting these into 10 windows of five values each, with a subsequent calculation of the MAD/Median ratio for each of these windows. This is the 2nd window for cycles (4.9,8.8]. \tab numeric \cr
#'   "window_Win_3" \tab The complete curve trajectory is segmented into 10 equidistant windows by fitting an interpolating smoothing spline with smoothing factor 0.5 to the curve, interpolating exactly 50 curve points, and then cutting these into 10 windows of five values each, with a subsequent calculation of the MAD/Median ratio for each of these windows. This is the 3rd window for cycles (8.8,12.7]. \tab numeric \cr
#'   "window_Win_4" \tab The complete curve trajectory is segmented into 10 equidistant windows by fitting an interpolating smoothing spline with smoothing factor 0.5 to the curve, interpolating exactly 50 curve points, and then cutting these into 10 windows of five values each, with a subsequent calculation of the MAD/Median ratio for each of these windows. This is the 4th window for cycles (12.7,16.6]. \tab numeric \cr
#'   "window_Win_5" \tab The complete curve trajectory is segmented into 10 equidistant windows by fitting an interpolating smoothing spline with smoothing factor 0.5 to the curve, interpolating exactly 50 curve points, and then cutting these into 10 windows of five values each, with a subsequent calculation of the MAD/Median ratio for each of these windows. This is the 5th window for cycles (16.6,20.5]. \tab numeric \cr
#'   "window_Win_6" \tab The complete curve trajectory is segmented into 10 equidistant windows by fitting an interpolating smoothing spline with smoothing factor 0.5 to the curve, interpolating exactly 50 curve points, and then cutting these into 10 windows of five values each, with a subsequent calculation of the MAD/Median ratio for each of these windows. This is the 6th window for cycles (20.5,24.4]. \tab numeric \cr
#'   "window_Win_7" \tab The complete curve trajectory is segmented into 10 equidistant windows by fitting an interpolating smoothing spline with smoothing factor 0.5 to the curve, interpolating exactly 50 curve points, and then cutting these into 10 windows of five values each, with a subsequent calculation of the MAD/Median ratio for each of these windows. This is the 7th window for cycles (24.4,28.3]. \tab numeric \cr
#'   "window_Win_8" \tab The complete curve trajectory is segmented into 10 equidistant windows by fitting an interpolating smoothing spline with smoothing factor 0.5 to the curve, interpolating exactly 50 curve points, and then cutting these into 10 windows of five values each, with a subsequent calculation of the MAD/Median ratio for each of these windows. This is the 8th window for cycles (28.3,32.2]. \tab numeric \cr
#'   "window_Win_9" \tab The complete curve trajectory is segmented into 10 equidistant windows by fitting an interpolating smoothing spline with smoothing factor 0.5 to the curve, interpolating exactly 50 curve points, and then cutting these into 10 windows of five values each, with a subsequent calculation of the MAD/Median ratio for each of these windows. This is the 9th window for cycles (32.2,36.1]. \tab numeric \cr
#'   "window_Win_10" \tab The complete curve trajectory is segmented into 10 equidistant windows by fitting an interpolating smoothing spline with smoothing factor 0.5 to the curve, interpolating exactly 50 curve points, and then cutting these into 10 windows of five values each, with a subsequent calculation of the MAD/Median ratio for each of these windows. This is the 10th window for cycles (36.1,40]. \tab numeric \cr
#'   "sd_plateau" \tab describes the standard deviation in the late phase of an amplification curve (last five cycles). With ideal PCRs, this corresponds to the plateau phase. \tab numeric \cr
#' }
#'
#' @return gives a \code{data.frame} (S3 class, type of \code{list}) as output 
#' for the curve features
#' 
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
#'  \code{\link[segmented]{segmented}}
#'
#' @rdname pcrfit_single
#' @export pcrfit_single

pcrfit_single <- function(x) {

  # Note: the input variable is x. This is changed within the following lines
  # into y, make the code more readable.

  # Remove exogenous missing values of the amplification data
  non_nas <- cumsum(!is.na(x))
  y <- x[1L:which.max(non_nas == max(non_nas))]

  # Determine the sequence of the cycles
  x <- 1L:length(y)

  # Impute endogenous missing values of the amplification data
  y <- armor(chipPCR::fixNA(x, y))

  # y => baselined
  y <- y - median(y[3:7], na.rm = TRUE)

  # Make every data set to have 50 data points in total

  res_spline <- smooth.spline(x, y, spar = 0.5)
  xseq <- seq(1, length(y), length.out = 50)
  ypred <- predict(res_spline, xseq)$y

  # Combine x and y by columns
  xy_tmp <- cbind(cyc = x, fluo = y)

  ######### Will be removed in later versions ##########
  # Normalize RFU values to the alpha percentile (0.99)
  # This normalization was exchanged by a base-lining in the following lines
  # to make calculations of the features like the efficiency more reliable
  y_99_quantile <- quantile(y, 0.99, na.rm = TRUE)
  if (y_99_quantile == 0) {
    y_quantile <- y / max(y, 0.99, na.rm = TRUE)
  } else {
    y_quantile <- y / quantile(y, 0.99, na.rm = TRUE)
  }
  ######### Will be removed in later versions ##########

  # Determine number of cycles used for the analysis

  length_cycle <- length(y)

  # Smooth data with moving average smoothing filter for further data
  # analysis steps.
  dat_smoothed <- chipPCR::smoother(x, y, method = "mova")

  # Guess direction of curve
  guess_direction <- function(y) ifelse(median(head(y, 5)) > (median(tail(y, 5)) + mad(tail(y, 5))), "max", "min")
  res_guess_direction <- armor(guess_direction(dat_smoothed))
  if (res_guess_direction == "min") res_guess_direction_bin <- 0 else res_guess_direction_bin <- 1

  # Calculate the first derivative
  res_diffQ <- armor(MBmca::diffQ(xy_tmp[-c(1:3), ], verbose = TRUE)$xy, 2)

  # Determine highest and lowest amplification curve values
  fluo_range <- stats::quantile(y, c(0.01, 0.99), na.rm = TRUE)
  fluo_range_abs <- round(as.numeric(diff(quantile(y, c(0.01, 0.99), na.rm = TRUE))), 2)

  # for qpcR
  dat_reverse <- cbind(cyc = x, fluo = (max(dat_smoothed) - rev(dat_smoothed)))
  # inefficient kludge to find l4 model
  data("sysdata", package = "qpcR", envir = parent.frame())

  # Estimate amplification start and stop via bg.max from chipPCR
  res_bg.max_tmp <- armor(chipPCR::bg.max(x, y), n = 3)
  if (length(res_bg.max_tmp) == 3) {
    res_bg.max <- c(
      bg.start = 0,
      bg.stop = 0,
      amp.stop = 0
    )
  } else {
    start_vals <- c(head(y_quantile, 5))
    stat_threshold <- mad(start_vals)
    if (stat_threshold > 0.2) {
      res_bg.max_tmp@bg.start <- 4
    } else {
      res_bg.max_tmp@bg.start <- 2
    }
    res_bg.max <- c(
      bg.start = res_bg.max_tmp@bg.start,
      bg.stop = res_bg.max_tmp@bg.stop,
      amp.stop = res_bg.max_tmp@amp.stop
    )
  }

  # Determine the head to tail ratio
  res_head2tail_ratio <- armor(PCRedux::head2tailratio(y))
  if (res_head2tail_ratio > 50 | !is.finite(res_head2tail_ratio)) res_head2tail_ratio <- 0

  # Determine the slope from the cycles 2 to 11
  res_lm_fit <- armor(PCRedux::earlyreg(x, y), 3)

  # Try to estimate the slope and the intercept of the tail region,
  # which might be indicative of a hook effect (strong negative slope)
  res_hookreg_simple <- armor(PCRedux::hookreg(x, y, normalize = FALSE), 10)

  res_hookregNL <- suppressMessages(PCRedux::hookregNL(x = x, y))
  res_hookregNL <- armor(PCRedux::hookregNL(x, y), 4)

  res_hookreg <- ifelse(res_hookreg_simple["hook"] == 1 || res_hookregNL["hook"] == 1, 1, 0)

  # Calculates the area of the amplification curve
  res_polyarea <- armor(round(pracma::polyarea(xseq, ypred), 3), 1)

  # Calculates the area of the amplification curve
  res_polyarea_trapz <- round(armor(pracma::trapz(xseq, ypred)), 3)

  # Calculate change points

  # Smooth data with moving average smoothing filter for further data
  # analysis steps.
  y_double <- c(y[-c(1:3)], y[-c(1:3)])
  cycles_double <- 1L:length(y_double)
  # Smooth data with moving average smoothing filter for further data
  # analysis steps.
  dat_smoothed_double <- cbind(cycles_double, chipPCR::smoother(cycles_double, y_double, method = "mova"))
  # Calculate the first derivative
  res_diffQ_double <- try(suppressMessages(MBmca::diffQ(dat_smoothed_double, verbose = TRUE)$xy), silent = TRUE)

  # Agglomerative hierarchical estimation algorithm for multiple change point analysis
  res_cp_e.agglo <- armor(length(ecp::e.agglo(as.matrix(res_diffQ_double[["d(F) / dT"]]))$estimates))

  # Bayesian analysis of change points
  res_bcp_tmp <- armor(bcp::bcp(res_diffQ_double[["d(F) / dT"]])$posterior.prob)
  if (length(res_bcp_tmp) != 1) {
    res_bcp <- armor(length(which(as.factor(res_bcp_tmp >= 0.95) == TRUE)))
  } else {
    res_bcp <- 0
  }

  # Median based local robust regression (mblrr)
  res_mblrr <- armor(PCRedux::mblrr(x, y), 6)

  # Calculate amptester results
  res_amptester <- armor(chipPCR::amptester(y), 7)
  if ((length(res_amptester) != 1) && as.character(class(res_amptester)) == "amptest") {
    amptester_shapiro <- res_amptester@decisions["shap.noisy"][[1]]
    amptester_lrt <- res_amptester@decisions["lrt.test"][[1]]
    amptester_rgt <- res_amptester@decisions["rgt.dec"][[1]]
    amptester_tht <- res_amptester@decisions["tht.dec"][[1]]
    amptester_slt <- res_amptester@decisions["slt.dec"][[1]]
    amptester_polygon <- res_amptester@polygon / length_cycle
    amptester_slope.ratio <- res_amptester@slope.ratio / length_cycle
    amptester <- c(
      amp.shap = as.numeric(amptester_shapiro[1]),
      amp.lrt = as.numeric(amptester_lrt[1]),
      amp.rgt = as.numeric(amptester_rgt[1]),
      amp.tht = as.numeric(amptester_tht[1]),
      amp.slt = as.numeric(amptester_slt[1]),
      amp.pol = as.numeric(amptester_polygon[1]),
      amp.ratio = as.numeric(amptester_slope.ratio)
    )
  } else {
    amptester <- c(
      amp.shap = 0,
      amp.lrt = 0,
      amp.rgt = 0,
      amp.tht = 0,
      amp.slt = 0,
      amp.pol = 0,
      amp.ratio = 0
    )
  }

  if (is.na(amptester["amp.shap"])) amptester["amp.shap"] <- 0
  if (is.na(amptester["amp.lrt"])) amptester["amp.lrt"] <- 0
  if (is.na(amptester["amp.rgt"])) amptester["amp.rgt"] <- 0
  if (is.na(amptester["amp.tht"])) amptester["amp.tht"] <- 0
  if (is.na(amptester["amp.pol"])) amptester["amp.pol"] <- 0
  if (is.na(amptester["amp.ratio"])) amptester["amp.ratio"] <- 0

  # Estimate the spread of the approximate local minima and maxima of the curve data

  res_mcaPeaks <- armor(MBmca::mcaPeaks(res_diffQ[, 1], res_diffQ[, 2]))
  if (length(res_mcaPeaks) == 2) {
    if (nrow(res_mcaPeaks[[1]]) == 0 || nrow(res_mcaPeaks[[2]]) == 0) {
      res_peaks_ratio <- 0
    } else {
      res_peaks_ratio <- base::diff(range(res_mcaPeaks$p.max[, 2])) / base::diff(range(res_mcaPeaks$p.min[, 2]))
      if (is.na(res_peaks_ratio) || res_peaks_ratio == "-Inf" || res_peaks_ratio == "Inf") {
        res_peaks_ratio <- 0
      }
    }
  } else {
    res_peaks_ratio <- res_mcaPeaks
  }

  # Estimate the slope between the minimum and the maximum of the second derivative
  res_diffQ2 <- armor(MBmca::diffQ2(cbind(x[-c(1:3)], dat_smoothed[-c(1:3)]),
    verbose = FALSE,
    fct = get(res_guess_direction, pos = "package:base"),
    inder = TRUE
  ))

  # The diffQ2 makes calculates the approximate Cq values approximately. In
  # some cases the second derivative maximum can be larger then the first
  # derivative maximum and in some cases the second derivative minimum can be
  # lower then the first derivative maximum. To account for this all values are
  # checked and set to the values at the approximate first derivative maximum
  # accordingly.

  if (res_diffQ2[["xTm1.2.D2"]][1] > res_diffQ2[["Tm"]]) {
    res_diffQ2[["xTm1.2.D2"]][1] <- res_diffQ2[["Tm"]]
    res_diffQ2[["yTm1.2.D2"]][1] <- res_diffQ2[["fluoTm"]]
  }
  if (res_diffQ2[["xTm1.2.D2"]][2] < res_diffQ2[["Tm"]]) {
    res_diffQ2[["xTm1.2.D2"]][2] <- res_diffQ2[["Tm"]]
    res_diffQ2[["yTm1.2.D2"]][2] <- res_diffQ2[["fluoTm"]]
  }

  # difference between the minimum and the maximum of the approximate second derivative.
  if (res_diffQ2[[3]][1] < res_diffQ2[1] && res_diffQ2[1] < res_diffQ2[[3]][2] && res_diffQ2[1] < length_cycle) {
    if (length(res_diffQ2) != 1) {
      cpD2_range <- as.numeric(armor(diff(res_diffQ2[[3]])))
      if (cpD2_range > 200) cpD2_range <- 0
      ROI <- armor(round(c(res_diffQ2[[1]], res_diffQ2[[3]]))[c(2, 1, 3)])
      res_loglin_slope <- armor(coefficients(lm(y[ROI] ~ ROI))[["ROI"]])
      if (is.na(res_loglin_slope)) res_loglin_slope <- 0
      if (inherits(res_loglin_slope, "try-error")) res_loglin_slope <- 0
    }
  } else {
    cpD2_range <- res_loglin_slope <- 0
  }

  # Perform an autocorrelation analysis

  res_autocorrelation <- armor(PCRedux::autocorrelation_test(y = dat_smoothed[-c(1:3)], n = 8))

  # Fitting of non-linear multiparameter models
  # Fit sigmoidal models to curve data

  b_val <- 0 # -10000
  c_val <- 0 # -10000
  d_val <- 0 # 10000
  e_val <- 0 # 10000
  f_val <- 0 # -10000
  k_val <- 0 # -10000
  k1_val <- 0 # -10000
  k2_val <- 0 # -10000

  # Get the parameters from the four-parameter model

  # The l4 model
  res_pcrfit <- armor(qpcR::pcrfit(xy_tmp, 1, 2, model = l4))
  if (length(res_pcrfit) != 1) {
    res_coef_pcrfit <- coefficients(res_pcrfit)
    res_AICc_l4 <- qpcR::AICc(res_pcrfit)
    res_iterations <- res_pcrfit[["convInfo"]][["finIter"]]
    # TakeOff Point
    # Calculates the first significant cycle of the exponential region
    # (takeoff point) using externally studentized residuals as described
    # in Tichopad et al. (2003).
    res_takeoff <- unlist(qpcR::takeoff(res_pcrfit))
    if (is.na(res_takeoff["top"])) res_takeoff["top"] <- res_takeoff["f.top"] <- 0
    # sliwin qPCR efficiency
    # Calculation of the qPCR efficiency by the 'window-of-linearity' method
    res_sliwin <- armor(qpcR::sliwin(res_pcrfit, plot = FALSE, verbose = FALSE)$eff)
    res_eff_tmp <- unlist(qpcR::efficiency(res_pcrfit, plot = FALSE))[c(1, 3, 5, 7, 8)]
    cpD2_ratio_tmp <- armor(as.numeric(res_eff_tmp[["cpD2"]] / res_diffQ2[[3]][1]))
  } else {
    res_coef_pcrfit <- c(0, 0, 0, 0)
    names(res_coef_pcrfit) <- c("b", "c", "d", "e")
    res_AICc_l4 <- res_iterations <- res_sliwin <- cpD2_ratio_tmp <- 0
    res_takeoff <- c(0, 0)
    res_eff_tmp <- c(0, 0, 0, 0, 0)
  }


  if (is.na(res_takeoff[[1]])) {
    bg_dat <- c(head(y_quantile, 8))
    res_sd_bg <- sd(bg_dat)
  } else {
    bg_dat <- y_quantile[1L:res_takeoff[[1]]]
    res_sd_bg <- armor(sd(bg_dat))
  }

  if (is.na(res_sd_bg) || res_sd_bg == "-Inf" || res_sd_bg == "-Inf") res_sd_bg <- 0

  pcrfit_model_l4 <- res_pcrfit
  res_coef_model_l4 <- armor(coefficients(pcrfit_model_l4))
  if (length(res_coef_model_l4) == 1) {
    res_coef_model_l4 <- c(b = b_val, c = c_val, d = d_val, e = e_val)
  }

  # Get the parameters from the five-parameter model
  pcrfit_model_l5 <- armor(qpcR::pcrfit(xy_tmp, 1, 2, model = l5))
  res_coef_model_l5 <- armor(coefficients(pcrfit_model_l5))
  if (length(res_coef_model_l5) == 1) {
    res_coef_model_l5 <- c(b = b_val, c = c_val, d = d_val, e = e_val, f = f_val)
  }

  # Get the parameters from the five-parameter model
  pcrfit_model_l6 <- armor(qpcR::pcrfit(xy_tmp, 1, 2, model = l6))
  res_coef_model_l6 <- armor(coefficients(pcrfit_model_l6))
  if (length(res_coef_model_l6) == 1) {
    res_coef_model_l6 <- c(b = b_val, c = c_val, d = d_val, e = e_val, f = f_val, k = k_val)
  }

  # Get the parameters from the seven-parameter model
  pcrfit_model_l7 <- armor(qpcR::pcrfit(xy_tmp, 1, 2, model = l7))
  res_coef_model_l7 <- armor(coefficients(pcrfit_model_l7))
  if (length(res_coef_model_l7) == 1) {
    res_coef_model_l7 <- c(b = b_val, c = c_val, d = d_val, e = e_val, f = f_val, k1 = k1_val, k2 = k2_val)
  }

  res_convInfo_iteratons <- armor(pcrfit_model_l7[["convInfo"]][["finIter"]])
  if (length(res_convInfo_iteratons) == 1) {
    res_convInfo_iteratons <- res_convInfo_iteratons
  }

  # Determine the optimal model based on the AICc

  res_AICc <- list(
    l4 = if (length(armor(pcrfit_model_l4)) == 1) {
      NA
    } else {
      armor(qpcR::AICc(pcrfit_model_l4))
    },
    l5 = if (armor(length(armor(pcrfit_model_l5)) == 1)) {
      NA
    } else {
      armor(qpcR::AICc(pcrfit_model_l5))
    },
    l6 = if (armor(length(armor(pcrfit_model_l6)) == 1)) {
      NA
    } else {
      armor(qpcR::AICc(pcrfit_model_l6))
    },
    l7 = if (armor(length(armor(pcrfit_model_l7)) == 1)) {
      NA
    } else {
      armor(qpcR::AICc(pcrfit_model_l7))
    }
  )

  res_fit_model <- as.factor(which(res_AICc == min(unlist(res_AICc), na.rm = TRUE)))

  if (length(res_fit_model) == 0) {
    res_fit_model <- as.factor("l0")
    names(res_fit_model) <- "l0"
  }

  if (names(res_fit_model) == "l4") res_fit <- pcrfit_model_l4
  if (names(res_fit_model) == "l5") res_fit <- pcrfit_model_l5
  if (names(res_fit_model) == "l6") res_fit <- pcrfit_model_l6
  if (names(res_fit_model) == "l7") res_fit <- pcrfit_model_l7
  if (names(res_fit_model) == "l0") res_fit <- "try-error"

  if (names(res_fit_model) == "l4") {
    res_optimal_coefficients <- c(res_coef_model_l4, f = f_val)
  }
  if (names(res_fit_model) == "l5") {
    res_optimal_coefficients <- c(res_coef_model_l5)
  }
  if (names(res_fit_model) == "l6") {
    res_optimal_coefficients <- c(res_coef_model_l6)
  }
  if (names(res_fit_model) == "l7") {
    res_optimal_coefficients <- c(res_coef_model_l7)
  }
  if (names(res_fit_model) == "l0") {
    res_optimal_coefficients <- c(b = b_val, c = c_val, d = d_val, e = e_val, f = f_val, k1 = k1_val, k2 = k2_val)
  }


  # Fit the "reverse" model
  pcrfit_startmodel_reverse <- armor(qpcR::pcrfit(dat_reverse, 1, 2, model = l4))

  res_fit_reverse <- armor(qpcR::mselect(
    pcrfit_startmodel_reverse,
    verbose = FALSE, fctList = list(l4, l5, l6, l7)
  ))

  # Determine the model for the reverse data suggested by the
  # mselect function based on the AICc

  if (length(res_fit_reverse) == 0) {
    res_fit_model_reverse <- as.factor("l0")
  } else {
    res_fit_model_reverse <- armor(names(which(res_fit_reverse[["retMat"]][, "AICc"] == min(res_fit_reverse[["retMat"]][, "AICc"]))))
  }

  if (res_fit_model_reverse != "l0") {
    # TakeOff Point from the reverse data
    # Calculates the first significant cycle of the exponential region
    #
    #   Take Down Point tdp
    #
    # (takeoff point) using externally studentized residuals as described
    # in Tichopad et al. (2003).
    # Calculate the standard deviation of the fluorescence starting from
    # cylce 2 to the takeoff point (res_sd_bg)
    res_takeoff_reverse <- armor(qpcR::takeoff(res_fit_reverse, nsig = 5), 2)
    if (is.na(res_takeoff_reverse[[1]])) {
      res_takeoff_reverse[[1]] <- res_takeoff_reverse[[2]] <- 0
      plateau_dat <- c(tail(y_quantile, 5))
      res_sd_plateau <- sd(plateau_dat)
    } else {
      res_takeoff_reverse[[1]] <- length_cycle - res_takeoff_reverse[[1]]
      res_takeoff_reverse[[2]] <- armor(predict(res_pcrfit, data.frame(Cycles = res_takeoff_reverse[[1]]))[[1]])
      plateau_dat <- y_quantile[res_takeoff_reverse[[1]]:length_cycle]
      res_sd_plateau <- sd(plateau_dat)
    }
  } else {
    res_takeoff_reverse[[1]] <- res_takeoff_reverse[[2]] <- 0
    plateau_dat <- c(tail(y_quantile, 5))
    res_sd_plateau <- sd(plateau_dat)
  }

  if (is.na(res_sd_plateau)) res_sd_plateau <- sd(tail(y_quantile, 5))

  names(res_takeoff_reverse) <- c("tdp", "f.tdp")

  # Cq of the amplification curve
  # Determine the Cq and other parameters
  res_efficiency_tmp <- armor(
    qpcR::efficiency(pcrfit_model_l4, plot = FALSE)[c(
      "eff",
      "cpD1",
      "cpD2",
      "fluo",
      "init2"
    )]
  )

  if (length(res_efficiency_tmp) == 5) {
    res_cpDdiff <- armor(res_efficiency_tmp[["cpD1"]] - res_efficiency_tmp[["cpD2"]])
  } else {
    res_efficiency_tmp <- list(
      eff = 0,
      cpD1 = 0,
      cpD2 = 0,
      fluo = 0,
      init2 = 0
    )
    res_cpDdiff <- 0
  }

  if (is.na(res_efficiency_tmp["init2"])) res_efficiency_tmp["init2"] <- 0

  # Calculate the ratio between the the approximate second derivative maximum
  # cpD2_approx and the second derivative maximum cpD2

  cpD2_ratio_tmp <- res_efficiency_tmp[["cpD2"]] / res_diffQ2[[3]][1]
  if (cpD2_ratio_tmp == Inf) cpD2_ratio <- 0
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

  if (is.na(res_angle)) res_angle <- 0

  ## Method 17
  res_cor <- armor(cor(x, y))

  ## Method 18
  ##############################################################################
#   # code taken literally from Package segmented version 1.2-0
#   # to avoid error message:
#   # Undefined global functions or variables:      seg.control
#   seg.control <- function(n.boot = 10, display = FALSE, tol = 1e-05, it.max = 30,
#                           fix.npsi = TRUE, K = 10, quant = TRUE, maxit.glm = 25, h = 1,
#                           size.boot = NULL, jt = FALSE, nonParam = TRUE, random = TRUE,
#                           seed = NULL, fn.obj = NULL, digits = NULL, conv.psi = FALSE,
#                           alpha = 0.02, min.step = 1e-04, powers = c(1, 1), last = TRUE,
#                           stop.if.error = NULL, gap = FALSE, fc = 0.95) {
#     list(
#       toll = tol, it.max = it.max, visual = display, stop.if.error = stop.if.error,
#       K = K, last = last, maxit.glm = maxit.glm, h = h, n.boot = n.boot,
#       size.boot = size.boot, gap = gap, jt = jt, nonParam = nonParam,
#       random = random, pow = powers, seed = seed, quant = quant,
#       fn.obj = fn.obj, digits = digits, conv.psi = conv.psi,
#       alpha = alpha, fix.npsi = fix.npsi, min.step = min.step,
#       fc = fc
#     )
#   }
  ##############################################################################

  control_fct <- segmented::seg.control(display = FALSE, n.boot = 20, size.boot = 20, it.max = 50)

  segmenter <- function(x, y, control_fct) {
    res_lm_segmenter <- lm(y ~ x)
    segments <- segmented::segmented(res_lm_segmenter, seg.Z = ~x, npsi = 2, control = control_fct)
    return(c(segments$coefficients[2:4], segments$psi[, 2]))
  }

  res_segment <- armor(segmenter(x, y), 5)

  if (length(which(is.na(res_segment) == TRUE)) >= 1) {
    res_segment <- rep(0, 5)
    names(res_segment) <- c("x", "U1.x", "U2.x", "psi1.x", "psi2.x")
  }

  ## Method 19
  sumdiffer <- function(y) sum(diff(y) > 0, na.rm = TRUE) / length(y)
  res_sumdiff <- armor(sumdiffer(y))

  ## Method 20
  polyreger <- function(x, y) {
    res_lm_polyreger <- lm(y ~ poly(x, 3))
    res_coef_polyreger <- coefficients(res_lm_polyreger)
    names(res_coef_polyreger) <- paste("poly_", 1:4, sep = "")
    return(res_coef_polyreger)
  }
  res_poly <- armor(polyreger(x, y), 4)

  # Takes any amplification curve and fits a spline to predict exactly 50 data
  # points for further processing.
  # The predicted data are cut in 10 segments containing 5 values. These values
  # are used to calculate the median and the MAD of these regions of interest.
  windower <- function(x, y) {
    res_spline <- smooth.spline(x, y, spar = 0.5)
    xseq <- seq(1, length(y), length.out = 50)
    ypred <- predict(res_spline, xseq)$y
    cut_y_values <- cut(xseq, 10)
    #     res_location_segments <- tapply(ypred, cut_y_values, function(x) median(x, na.rm = TRUE))
    res_variation_segments <- tapply(ypred, cut_y_values, function(x) mad(x, na.rm = TRUE))
    res_windower <- res_variation_segments # / res_location_segments
    names(res_windower) <- paste("Win_", 1:10, sep = "")
    return(res_windower)
  }

  res_window <- armor(windower(x, y), 10)

  all_results <- data.frame(
    # Quantification points, derivatives, efficiencies,
    # curve fitting
    cpD1 = res_efficiency_tmp[["cpD1"]],
    cpD2 = res_efficiency_tmp[["cpD2"]],
    cpD2_approx = res_diffQ2[[3]][1],
    cpD2_ratio = cpD2_ratio,
    eff = res_eff_tmp[1], # res_efficiency_tmp[["eff"]],
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
    b_slope = res_coef_model_l7[["b"]],
    b_model_param = res_optimal_coefficients[["b"]],
    c_model_param = res_optimal_coefficients[["c"]],
    d_model_param = res_optimal_coefficients[["d"]],
    e_model_param = res_optimal_coefficients[["e"]],
    f_model_param = res_optimal_coefficients[["f"]],
    f_intercept = res_coef_model_l7[["f"]],
    k1_model_param = res_coef_model_l7[["k1"]],
    k2_model_param = res_coef_model_l7[["k2"]],
    convInfo_iteratons = res_convInfo_iteratons,
    qPCRmodel = factor(names(res_fit_model)), # res_fit_model[[1]],
    qPCRmodelRF = factor(res_fit_model_reverse[[1]]),
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
    amptester_shapiro = amptester["amp.shap"], # res_amptester@decisions["shap.noisy"][[1]],
    amptester_lrt = amptester["amp.lrt"], #  res_amptester@decisions["lrt.test"][[1]],
    amptester_rgt = amptester["amp.rgt"], #  res_amptester@decisions["rgt.dec"][[1]],
    amptester_tht = amptester["amp.tht"], #  res_amptester@decisions["tht.dec"][[1]],
    amptester_slt = amptester["amp.slt"], #  res_amptester@decisions["slt.dec"][[1]],
    amptester_polygon = amptester["amp.pol"], #  res_amptester@"polygon" / length_cycle,
    amptester_slope.ratio = amptester["amp.ratio"], #  ifelse(is.na(res_amptester@"slope.ratio"), 0, res_amptester@"slope.ratio"),
    # Curvature
    hookreg_hook = res_hookreg,
    hookreg_hook_slope = res_hookreg_simple[["slope"]],
    hookreg_hook_delta = res_hookreg_simple[["hook.delta"]],
    # Angle
    central_angle = res_angle,
    # Number of Cycles
    number_of_cycles = length_cycle,
    # Identifier
    # October 2020 update
    direction = res_guess_direction_bin,
    range = fluo_range_abs,
    polyarea_trapz = res_polyarea_trapz,
    cor = res_cor,
    res_coef_pcrfit.b = res_coef_pcrfit["b"],
    res_coef_pcrfit.c = res_coef_pcrfit["c"],
    res_coef_pcrfit.d = res_coef_pcrfit["d"],
    res_coef_pcrfit.e = res_coef_pcrfit["e"],
    fitAIC = res_AICc_l4,
    fitIter = res_iterations,
    segment_x = res_segment[[1]],
    segment_U1.x = res_segment[[2]],
    segment_U2.x = res_segment[[3]],
    segment_psi1.x = res_segment[[4]],
    segment_psi2.x = res_segment[[5]],
    sumdiff = res_sumdiff,
    poly_1 = res_poly[["poly_1"]],
    poly_2 = res_poly[["poly_2"]],
    poly_3 = res_poly[["poly_3"]],
    poly_4 = res_poly[["poly_4"]],
    window_Win_1 = res_window[["Win_1"]],
    window_Win_2 = res_window[["Win_2"]],
    window_Win_3 = res_window[["Win_3"]],
    window_Win_4 = res_window[["Win_4"]],
    window_Win_5 = res_window[["Win_5"]],
    window_Win_6 = res_window[["Win_6"]],
    window_Win_7 = res_window[["Win_7"]],
    window_Win_8 = res_window[["Win_8"]],
    window_Win_9 = res_window[["Win_9"]],
    window_Win_10 = res_window[["Win_10"]],
    sd_plateau = res_sd_plateau,
    row.names = "results"
  )
  return(all_results)
}
