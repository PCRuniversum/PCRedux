#' @title Visualises the content of data from an analysis with the pcrfit_parallel function
#' @description gives a `ggplot` or interactive `ggplotly` object (default) of a 'pcrfit_parallel' data frame. The function is based on the \link[visdat]{vis_dat} function by Tierney (2017). The colored indicators show the class (e.g., missing values, factors).
#' @param data contains the result from an analysis with the \link[PCRedux]{pcrfit_parallel} function
#' @param type specifies of all or only a subset of results should be shown, 
#' Default: 'all'. Alternatives 'qpcR' or 'amptester'
#' @param interactive is a logical parameter, which indicates if the plot 
#' should be interactive, Default: TRUE
#' @author Stefan Roediger, Michal Burdukiewcz
#' @references N. Tierney, visdat: Visualising Whole Data Frames, \emph{The 
#' Journal of Open Source Software}. 2 (2017). doi:10.21105/joss.00355.
#' @details 
#' type 'all' shows all results from the analysis by the 
#' \link[PCRedux]{pcrfit_parallel} function. 
#'
#' @examples
#'
#' # Calculate curve features of an amplification curve data. Note that not all 
#' # available CPU cores are used. If need set "all" to use all available cores.
#' library(qpcR)
#' res_pcrfit_parallel <- pcrfit_parallel(boggy[, 1:7])
#' res_pcrfit_parallel
#'
#' # Show all results in an interactive plot
#' visdat_pcrfit_parallel(res_pcrfit_parallel, type="all")
#' @seealso 
#'  \code{\link[plotly]{ggplotly}}

#'  \code{\link[visdat]{vis_dat}}
#' @rdname visdat_pcrfit_parallel
#' @export 
#' @importFrom plotly ggplotly
#' @importFrom visdat vis_dat

visdat_pcrfit_parallel <- function(data, type="all", interactive=TRUE){
    
    fetch_all <- c("runs", "eff", "cpD1", "cpD2", "fluo", "init1", "init2", "top", 
                    "f.top", "resLRE", "ressliwin", "cpDdiff", "slope_background", 
                    "intercept_background", "polyarea", "changepoint.e.agglo", "changepoint.bcp", 
                    "qPCRmodel", "amptester_shap.noisy", "amptester_lrt.test", "amptester_rgt.dec", 
                    "amptester_tht.dec", "amptester_slt.dec", "amptester_polygon", 
                    "amptester_slope.ratio", "minRFU", "maxRFU", "bg.start_normalized", 
                    "bg.stop_normalized", "amp.stop_normalized", "head_to_tail_ratio", 
                    "autocorellation", "mblrr_intercept_less", "mblrr_slope_less", 
                    "mblrr_cor_less", "mblrr_intercept_more", "mblrr_slope_more", 
                    "mblrr_cor_more", "hookreg_hook", "mcaPeaks_minima_maxima_ratio", 
                    "diffQ2_slope", "diffQ2_Cq_range", "detection_chemistry", "device")
                    
 
    fetch_qpcR <- c("runs", "eff", "cpD1", "cpD2", "fluo", "init1", "init2", "top", 
                    "f.top", "resLRE", "ressliwin", "cpDdiff", "qPCRmodel")
    
    fetch_amptester <- c("runs", "slope_background", 
                    "intercept_background", "amptester_shap.noisy", "amptester_lrt.test", "amptester_rgt.dec", 
                    "amptester_tht.dec", "amptester_slt.dec", "amptester_polygon", 
                    "amptester_slope.ratio", "minRFU", "maxRFU", "bg.start_normalized", 
                    "bg.stop_normalized", "amp.stop_normalized", "head_to_tail_ratio", 
                    "mblrr_intercept_less", "mblrr_slope_less", 
                    "mblrr_cor_less", "mblrr_intercept_more", "mblrr_slope_more", 
                    "mblrr_cor_more", "hookreg_hook", "mcaPeaks_minima_maxima_ratio", 
                    "diffQ2_slope", "diffQ2_Cq_range")
    
    
    
    switch(type,
        all={
             data_tmp <- data[, fetch_all]
        },
        qpcR={
             data_tmp <- data[, fetch_qpcR]
        },
        amptester={
             data_tmp <- data[, fetch_amptester]
        },
        stop("Enter all, qpcR or amptester.")
    )
         
    plot_vis_dat <- visdat::vis_dat(data_tmp)
    
    if(interactive) {
        suppressMessages(plotly::ggplotly(plot_vis_dat))
        } else {
          plot_vis_dat
        }
}
