[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/PCRedux)](https://cran.r-project.org/package=PCRedux)
[![Downloads](http://cranlogs.r-pkg.org/badges/PCRedux)](https://cran.r-project.org/package=PCRedux)
[![R-CMD-check](https://github.com/PCRuniversum/PCRedux/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PCRuniversum/PCRedux/actions/workflows/R-CMD-check.yaml)

![PCRedux](https://raw.githubusercontent.com/PCRuniversum/PCRedux/master/vignettes/Logo.png)

# PCRedux

Quantitative PCR Machine Learning Toolkit

The software extracts features from amplification curve data of quantitative Polymerase Chain Reactions (qPCR) ([Pabinger S. et al., 2014]( https://doi.org/10.1016%2Fj.bdq.2014.08.002)) for machine learning purposes. Helper functions prepare the amplification curve data for processing as functional data (e.g., Hausdorff distance) or enable the plotting of amplification curve classes (negative, ambiguous, positive). The hookreg() and hookregNL() functions ([Burdukiewicz M. et al., 2018](https://doi.org/10.1016%2Fj.bdq.2018.08.001)) can be used to predict amplification curves with an hook effect-like curvature. The pcrfit_single() function can be used to extract features from an amplification curve.

## Installation

*PCRedux* is available [on CRAN](https://cran.r-project.org/package=PCRedux). However, you 
can install the latest development version of the code using the following code:

```R
library("devtools")
install_github("PCRuniversum/PCRedux")
```

## Manual

The manual is available [online](https://PCRuniversum.github.io/PCRedux/).
