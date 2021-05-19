# Dirac Delta Regression (DDR)

This is an R package implementing DDR, an algorithm that transforms the response variable of regression into a set of asymptotically Dirac delta functions using kernel density functions. This allows the user to convert a non-linear regressor to a conditional density estimator. We use kernel ridge regression as the underlying regressor in this implementation.

# Installation

Please install the `FNN`, `pracma`, `doParallel` and `foreach` packages on CRAN. Then:

> library(devtools)

> install_github("ericstrobl/DDR")

> library(DDR)

> numCores <- detectCores()-1; registerDoParallel(numCores) # set up parallel computing

> cd_est = DDR(matrix(rnorm(400),200,2),rnorm(200),matrix(rnorm(20),10,2)) # run DDR

> plot(cd_est$y,cd_est$dens[1,],type="l") # plot the conditional density estimate of the first test sample
