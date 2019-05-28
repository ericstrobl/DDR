# Dirac Delta Regression (DDR)

This is an R package implementing DDR, an algorithm that transforms the response variable of regression into a set of asymptotically Dirac delta functions using kernel density functions. This allows the user to convert a non-linear regressor to a conditional density estimator. We use kernel ridge regression with a Gaussian reproducing kernel as the underlying regressor in this implementation.

The academic article describing DDR in detail can be found [here](https://arxiv.org/abs/1905.10330). Please cite the article if you use any of the code in this repository.

# Installation

Please install the `pracma`, `doParallel`, and `foreach` packages on CRAN. Then:

> library(devtools)

> install_github("ericstrobl/DDR")

> library(DDR)

> numCores <- detectCores()-1; registerDoParallel(numCores)

> cd_est = DDR(matrix(rnorm(200),100,2),rnorm(100),matrix(rnorm(20),10,2))
