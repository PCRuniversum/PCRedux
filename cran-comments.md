# Note

- The previous version of the package was taken from CRAN because an e-mail address in the chiPPCR package had reached EOL. This caused the arching of the MBmca, the chipPCR and the PCRedux package. Both the MBmca and the chipPCR package are back on CRAN.
- The package is feature complete, parts of it were published (peer-reviewed) and used by peers use it. Therefore we change to version 1.0+
- Fixes of this submission include 
    - Grammar and spelling corrections
    - The PCRedux package now has a cran-comments.md file
    - The PCRedux package has been tested with rhub
    - In DESCRIPTION, all authors now have an ORCID
    - The documentation of the pcrfit_single() function was improved
    - RoxygenNote is now on version 7.1.1 (previously 5.0.1)
    - Citations in the vignettes were revised
    - We make sure that we do not change the user's options, par or working directory
    - NEWS updated

## Test environments
* local R installation, R 4.0.3
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

# check_rhub(env_vars=c(R_COMPILE_AND_INSTALL_PACKAGES = "always"))
✔  checking for file ‘/home/tux/Work/paper/PCRedux_all/PCRedux/DESCRIPTION’
─  preparing ‘PCRedux’: (678ms)
✔  checking DESCRIPTION meta-information
─  installing the package to build vignettes
✔  creating vignettes (9s)
─  checking for LF line-endings in source and make files and shell scripts
─  checking for empty or unneeded directories
   Removed empty directory ‘PCRedux/inst/PCRedux-app’
─  building ‘PCRedux_1.1.tar.gz’
   
─  Uploading package
─  Preparing build, see status at
   https://builder.r-hub.io/status/PCRedux_1.1.tar.gz-c435c6b161e74ed99b89a493aa3a9fdd
   https://builder.r-hub.io/status/PCRedux_1.1.tar.gz-4783440ec72c4c21a6e33783a6bfc542
   https://builder.r-hub.io/status/PCRedux_1.1.tar.gz-4eb0289c199e4ccc8f45a8986d888e5a
─  Build started
─  Creating new user
─  Downloading and unpacking package file
─  Querying package dependencies
─  Installing package dependencies
─  Running R CMD check
   setting _R_CHECK_FORCE_SUGGESTS_ to false
   setting R_COMPILE_AND_INSTALL_PACKAGES to never
   setting _R_CHECK_THINGS_IN_CHECK_DIR_ to false
   setting R_REMOTES_STANDALONE to true
   setting R_REMOTES_NO_ERRORS_FROM_WARNINGS to true
   setting R_COMPILE_AND_INSTALL_PACKAGES to always
─  using log directory 'C:/Users/USEREbGBlhkWMf/PCRedux.Rcheck'
─  using R Under development (unstable) (2021-02-15 r80013)
─  using platform: x86_64-w64-mingw32 (64-bit)
─  using session charset: ISO8859-1
─  using option '--as-cran'
✔  checking for file 'PCRedux/DESCRIPTION'
─  checking extension type ... Package
─  this is package 'PCRedux' version '1.1'
─  package encoding: UTF-8
N  checking CRAN incoming feasibility
   Maintainer: 'Stefan Roediger <stefan.roediger@b-tu.de>'
   
   New submission
   
   Package was archived on CRAN
   
   CRAN repository db overrides:
     X-CRAN-Comment: Archived on 2020-12-16 as requires archived package
       'chipPCR'.
✔  checking package namespace information
✔  checking package dependencies
✔  checking if this is a source package
✔  checking if there is a namespace
✔  checking for executable files
✔  checking for hidden files and directories
✔  checking for portable file names
✔  checking whether package 'PCRedux' can be installed
✔  checking installed package size
✔  checking package directory
✔  checking for future file timestamps
✔  checking 'build' directory
✔  checking DESCRIPTION meta-information
✔  checking top-level files
✔  checking for left-over files
✔  checking index information
✔  checking package subdirectories
✔  checking R files for non-ASCII characters
✔  checking R files for syntax errors
✔  checking whether the package can be loaded
✔  checking whether the package can be loaded with stated dependencies
✔  checking whether the package can be unloaded cleanly
✔  checking whether the namespace can be loaded with stated dependencies
✔  checking whether the namespace can be unloaded cleanly
✔  checking loading without being on the library search path
✔  checking use of S3 registration
✔  checking dependencies in R code
✔  checking S3 generic/method consistency
✔  checking replacement functions
✔  checking foreign function calls
✔  checking R code for possible problems
✔  checking Rd files
✔  checking Rd metadata
✔  checking Rd line widths
✔  checking Rd cross-references
✔  checking for missing documentation entries
✔  checking for code/documentation mismatches
✔  checking Rd \usage sections
✔  checking Rd contents
✔  checking for unstated dependencies in examples
✔  checking contents of 'data' directory
✔  checking data for non-ASCII characters
✔  checking data for ASCII and uncompressed saves
N  checking sizes of PDF files under 'inst/doc'
   Unable to find GhostScript executable to run checks on size reduction
✔  checking installed files from 'inst/doc'
✔  checking files in 'vignettes'
✔  checking examples
✔  checking for unstated dependencies in 'tests'
─  checking tests
✔  Running 'spelling.R'
✔  Running 'testthat.R'
✔  checking for unstated dependencies in vignettes
✔  checking package vignettes in 'inst/doc'
✔  checking re-building of vignette outputs
✔  checking PDF version of manual
✔  checking for non-standard things in the check directory
✔  checking for detritus in the temp directory
   
─  Done with R CMD check
─  Cleaning up files and user
    

── PCRedux 1.1: NOTE

  Build ID:   PCRedux_1.1.tar.gz-c435c6b161e74ed99b89a493aa3a9fdd
  Platform:   Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  Submitted:  6m 37.6s ago
  Build time: 6m 22.5s

❯ checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Stefan Roediger <stefan.roediger@b-tu.de>'
  
  New submission
  
  Package was archived on CRAN
  
  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2020-12-16 as requires archived package
      'chipPCR'.

❯ checking sizes of PDF files under 'inst/doc' ... NOTE
  Unable to find GhostScript executable to run checks on size reduction

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

── PCRedux 1.1: CREATED

  Build ID:   PCRedux_1.1.tar.gz-4783440ec72c4c21a6e33783a6bfc542
  Platform:   Ubuntu Linux 20.04.1 LTS, R-release, GCC
  Submitted:  6m 37.6s ago


── PCRedux 1.1: CREATED

  Build ID:   PCRedux_1.1.tar.gz-4eb0289c199e4ccc8f45a8986d888e5a
  Platform:   Fedora Linux, R-devel, clang, gfortran
  Submitted:  6m 37.6s ago
