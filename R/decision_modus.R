#' A function to get a decision (modus) from a vector of classes
#'
#' \code{decision_modus} is a function that can be used to find the most frequent
#' (modus) decision. The classes can be defined by the user (e.g., a", "n", "y"
#' -> "ambiguous", "negative", "positive"). This function is useful if large
#' collections of varying decision (e.g., "a", "a", "a", "n", "n") need to be
#' condensed to a single decision (3 x "a", 2 x "n" -> "a").
#' 
#' @return gives a \code{factor} (S3 class, type of \code{integer}) as output 
#' for a decision
#'
#' @param data is a table containing the classes.
#' @param variables is the class to look for.
#' @param max_freq is a logical parameter (default == TRUE) delivers either the
#' most occurring class or a summary.
#' @author Stefan Roediger, Michal Burdukiewcz
#' @keywords decision modus
#' @examples
#'
#' # First example
#' # Enter a string of arbritary of "a","a","y","n"
#' # Result:
#' # [1] a
#' # Levels: a b n y
#'
#' decision_modus(c("a","a","y","n","b"))
#'
#' # Second example
#' # Analyze data from the decision_res_testdat.csv data file
#' filename <- system.file("decision_res_testdat.csv", package = "PCRedux")
#
#' my_data <- read.csv(filename)
#' head(my_data)
#'
#' dec <- unlist(lapply(1L:nrow(my_data), function(i) {
#'        decision_modus(my_data[i, 2:4])
#' }))
#'
#' names(dec) <- my_data[, 1]
#' dec
#'
#' @export decision_modus

decision_modus <- function(data, variables=c("a", "n", "y"), max_freq=TRUE) {
  # read in data and unlist them for the processing
  unlisted_data <- unlist(data)
  #  find the unique elements
  unique_variables <- unique(unlisted_data)
  # traverse over the vector with the decision elements and apply the sum
  # function to the to get the total number for each decision element.
  sum_unique_variables <- sapply(1L:length(unique_variables), function(i) {
    sum(unlisted_data %in% unique_variables[i])
  })
  
  # Make a data frame wiht the decision and their frequencies.
  freq_df <- data.frame(
    variable = as.character(unique_variables),
    freq = sum_unique_variables, 
    stringsAsFactors = TRUE
  )
  
  # Perform a logical operation on the summarized decision elements.
  # Either report the most common element or total statistics about the decision
  if (max_freq) {
    
    # Report the most frequent decision only
    most_frequent_decision <- freq_df[freq_df[, "freq"] == max(freq_df[, "freq"]), "variable"]
    
    if(length(most_frequent_decision) > 1) {
      return(NA)
    } else {
      return(most_frequent_decision)
    }
  } else {
    freq_df
  }
}
