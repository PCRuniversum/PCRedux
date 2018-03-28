#' @title Visualizing the content of data from an analysis with the pcrfit_single function
#' @description gives a `ggplot` or interactive `ggplotly` object (default) of a 'pcrfit_single' data frame. The function is based on the \link[visdat]{vis_dat} function by Tierney (2017). The colored indicators show the class (e.g., missing values, factors).
#' @param data contains the result from an analysis with the \link[PCRedux]{pcrfit_single} function
#' @param type specifies of all or only a subset of results should be shown,
#' Default: 'all'. Alternatives 'qpcR' or 'amptester'
#' @param interactive is a logical parameter, which indicates if the plot
#' should be interactive, Default: TRUE
#' @author Stefan Roediger, Michal Burdukiewcz
#' @references N. Tierney, visdat: Visualising Whole Data Frames, \emph{The
#' Journal of Open Source Software}. 2 (2017). doi:10.21105/joss.00355.
#' @details
#' type 'all' shows all results from the analysis by the
#' \link[PCRedux]{pcrfit_single} function.
#'
#' @examples
#'
#' # Calculate curve features of an amplification curve data. Note that not all
#' # available CPU cores are used. If need set "all" to use all available cores.
#' library(qpcR)
#' # take the samples F1.1 (positive) and F1.3 (negative) for this example.
#'
#' test_data <- testdat[, c(1,2,4)]
#'
#' # Plot the amplification curves
#'
#' matplot(test_data[, 1], test_data[, -1], xlab="Cycle", ylab="RFU",
#'         main="testdat data set", type="l", lty=1, lwd=2, col=1:2)
#' legend("topleft", paste(colnames(test_data)[-1], c("pos", "neg")),
#'        pch=19, col=1:2)
#'
#' # Analyze the amplification curves with the pcrfit_single function
#' res_1 <- cbind(runs="F1.1", pcrfit_single(test_data[, 2]))
#' res_2 <- cbind(runs="F1.3", pcrfit_single(test_data[, 3]))
#' res <- rbind(F1.1=res_1, F1.3=res_2)
#'
#' # Show all results in an interactive plot
#' visdat_pcrfit(res)
#' @seealso
#'  \code{\link[plotly]{ggplotly}}

#'  \code{\link[visdat]{vis_dat}}
#' @rdname visdat_pcrfit
#' @export
#' @importFrom plotly ggplotly
#' @importFrom visdat vis_dat

visdat_pcrfit <- function(data, type="all", interactive=TRUE) {
  fetch_all <- c(
    "runs", "eff", "cpD1", "cpD2", "fluo", "init2", "top",
#   "runs", "eff", "cpD1", "cpD2", "fluo", "init1", "init2", "top",
    "f.top", "sliwin", "cpDdiff", "slope_bg",
    "intercept_bg", "polyarea", "changepoint_e.agglo", "changepoint_bcp",
    "qPCRmodel", "amptester_shapiro", "amptester_lrt", "amptester_rgt",
    "amptester_tht", "amptester_slt", "amptester_polygon",
    "amptester_slope.ratio", "minRFU", "maxRFU",
    "bg.stop", "amp.stop", "head2tail_ratio",
    "autocorellation", "mblrr_intercept_bg", "mblrr_slope_bg",
    "mblrr_cor_bg", "mblrr_intercept_pt", "mblrr_slope_pt",
    "mblrr_cor_pt", "hookreg_hook", "peaks_ratio",
    "loglin_slope", "cpD2_range"
  )


  fetch_qpcR <- c(
    "runs", "eff", "cpD1", "cpD2", "fluo", "init2", "top",
#   "runs", "eff", "cpD1", "cpD2", "fluo", "init1", "init2", "top",
    "f.top", "sliwin", "cpDdiff", "qPCRmodel"
  )

  fetch_amptester <- c(
    "runs", "slope_bg",
    "intercept_bg", "amptester_shapiro", "amptester_lrt", "amptester_rgt",
    "amptester_tht", "amptester_slt", "amptester_polygon",
    "amptester_slope.ratio", "minRFU", "maxRFU",
    "bg.stop", "amp.stop", "head2tail_ratio",
    "mblrr_intercept_bg", "mblrr_slope_bg",
    "mblrr_cor_bg", "mblrr_intercept_pt", "mblrr_slope_pt",
    "mblrr_cor_pt", "hookreg_hook", "peaks_ratio",
    "loglin_slope", "cpD2_range"
  )



  switch(type,
    all = {
      data_tmp <- data[, fetch_all]
    },
    qpcR = {
      data_tmp <- data[, fetch_qpcR]
    },
    amptester = {
      data_tmp <- data[, fetch_amptester]
    },
    stop("Enter all, qpcR or amptester.")
  )

  plot_vis_dat <- visdat::vis_dat(data_tmp)

  if (interactive) {
    suppressMessages(plotly::ggplotly(plot_vis_dat))
  } else {
    plot_vis_dat
  }
}
