my_DT <- function(x)
  datatable(x, escape = FALSE, extensions = 'Buttons',
            filter = "top", rownames = FALSE, 
            options = list(scrollX = TRUE,
                           dom = "Brtip",
                           buttons = c("copy", "csv", "excel", "print"),
                           pageLength = 50,
                           digits = 4)) %>% 
  formatStyle(1L:ncol(x), color = "black")

# shiny_encu -------------------------------------------
# a shiny-compatible version of encu https://github.com/devSJR/PCRedux/blob/master/R/encu.R

shiny_encu <- function(data, detection_chemistry = NA, device = NA) {
  # Determine the number of available cores and register them
  
  # Prepare the data for further processing
  # Normalize RFU values to the alpha quantiles (0.999)
  cycles <- data.frame(cycles = data[, 1])
  data_RFU <- data.frame(data[, -1])
  ncol_data_RFU <- ncol(data_RFU)
  data_RFU_colnames <- colnames(data_RFU)
  data_RFU <- sapply(1L:ncol_data_RFU, function(i) {
    data_RFU[, i] / quantile(data_RFU[, i], 0.999, na.rm = TRUE)
  })
  colnames(data_RFU) <- data_RFU_colnames

  run_res <- do.call(rbind, lapply(1L:ncol_data_RFU, function(ith_run) {
    incProgress(1/ncol_data_RFU)
    pcrfit_single(data_RFU[, ith_run])
    }))

  res <- cbind(runs = colnames(data_RFU), run_res)
  rownames(res) <- NULL
  res
}
