#' A function to get a decision (modus) from a vector of classes
#' 
#' \code{decision_modus} is a function that can be used to find the most frequent
#' (modus) decision. The classes can be defined by the user (e.g., a", "n", "y" 
#' -> "ambiguos", "negaive", "positive"). This function is usefuel if large 
#' collections of varying decision (e.g., "a", "a", "a", "n", "n") need to be 
#' condensed to a single decision (3 x "a", 2 x "n" -> "a").
#' 
#' @param data is a table contining the classes.
#' @param variables is the class to look for.
#' @param max_freq delivers either the most ocurring class or a summary.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords decision modus
#' @examples
#' 
#' #
#' library(data.table)
#' library(magrittr)
#' filename <- system.file("decision_res_testdat.csv", package = "PCRedux")
#  
#' my_data <- as.data.frame(fread(filename))
#' head(my_data)
#'
#' dec <- lapply(1L:nrow(my_data), function(i) {
#'        decision_modus(my_data[i, 2:4])
#' }) %>% unlist
#' 
#' names(dec) <- my_data[, 1]
#' dec
#' 
#' @export decision_modus

decision_modus <- function(data, variables=c("a", "n", "y"), max_freq=TRUE){
    # read in data and unlist them for the processing
    unlisted_data <- unlist(data)
    #  find the unique elements
    unique_variables <- unlisted_data %>% unique
    # traverse over the vector with the decision elements and apply the sum 
    # function to the to get the total number for each decision element.
    sum_unique_variables <- sapply(1L:length(unique_variables), function(i) {
        unlisted_data %in% unique_variables[i] %>% sum
    })
    # Perform a logical operation on the summarized decision elements.
    # Either report the most common element or total statistics about the decision
    if(max_freq) {
        # Make a data frame wiht the decision and their frequencies.
        tmp <- data.frame(variable=as.character(unique_variables), 
                          freq=sum_unique_variables)
        # Report the most frequent decision only
        tmp[tmp[, "freq"] == max(tmp[, "freq"]), "variable"]
    } else{
        data.frame(variable=as.character(unique_variables), freq=sum_unique_variables)
    }
}
