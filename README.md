# Dirac Delta Regression (DDR)

This is an R package implementing DDR, an algorithm that transforms the response variable of regression into a set of asymptotically Dirac delta functions using kernel density functions. This allows the user to convert a non-linear regressor to a conditional density estimator. This implementation uses kernel ridge regression with a reproducing Gaussian kernel as the underlying regressor.

The academic article describing DDR in detail can be found here. Please cite the article if you use any of the code in this repository.

# Installation

Please install the `pracma`, `doParallel`, and `foreach` packages on CRAN. Then:

> install_github("ericstrobl/DDR")

> library(DDR)

> DDR(rnorm(1000),rnorm(1000),rnorm(1000))

> RCoT(rnorm(1000),rnorm(1000),rnorm(1000))

 The DDR algorithm accepts matrices with rows and columns corresponding to samples and features, respectiely.
