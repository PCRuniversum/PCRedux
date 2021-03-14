# Add this code example code to 'shiny' user .Rprofile on Shiny server
# to activate source of this file and use remote 'pcrfit_single'.
# Also set ssh access to remote server by key without passphrase.
# library(future)
# initRemote <- function() {
#   login <- tweak(remote, workers = "user@server")
#   plan(list(
#     login,
#     multiprocess
#   ))
# }

# Run remote server connection 
initRemote()

# Overload of calcRunRes() function to use inside shiny_encu()
library(listenv)
calcRunRes <- function(ncol_data_RFU, data_RFU) {
  value(
    future(
    {
      y <- listenv()
      for (ith_run in 1L:ncol_data_RFU) {
        y[[ith_run]] %<-% {
          pcrfit_single(data_RFU[, ith_run])
        }
      }
      Reduce(rbind, y)
    })
  )
}
