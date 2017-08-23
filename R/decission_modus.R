#' A function to get a final decission (modus) from a vector of classes
#' 
#' \code{decission_modus} is a function to to get a final decission (modus) 
#' from a vector of classes.
#' 
#' @param data is a table contining the classes.
#' @param variables is the class to look for.
#' @param max_freq delivers either the most ocurring class or a summary.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords decission modus
#' @examples
#' 
#' # #'
#' library(data.table)
#' library(magrittr)
#' filename <- system.file("decission_res_testdat.csv", package = "PCRedux")
#  
#' my_data <- as.data.frame(fread(filename))
#' head(my_data)
#'
#' dec <- lapply(1L:nrow(my_data), function(i) {
#'        decission_modus(my_data[i, 2:4])
#' }) %>% unlist
#' 
#' names(dec) <- my_data[, 1]
#' dec
#' @export decission_modus

decission_modus <- function(data, variables=c("a", "n", "y"), max_freq=TRUE){
    unlisted_data <- unlist(data)
    unique_variables <- unlisted_data %>% unique
    sum_unique_variables <- sapply(1L:length(unique_variables), function(i) {
        unlisted_data %in% unique_variables[i] %>% sum
    })
    
    if(max_freq) {
        tmp <- data.frame(variable=as.character(unique_variables), freq=sum_unique_variables)
        tmp[tmp[, "freq"] == max(tmp[, "freq"]), "variable"]
    } else{
        data.frame(variable=as.character(unique_variables), freq=sum_unique_variables)
    }
}
