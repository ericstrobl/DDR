# Dirac Delta Regression (DDR)

This is an R package implementing DDR, an algorithm that transforms the response variable of regression into a set of asymptotically Dirac delta functions using kernel density functions. This allows the user to convert a non-linear regressor to a conditional density estimator. We use kernel ridge regression as the underlying regressor in this implementation.

# Installation

Please install the `FNN`, `pracma`, `doParallel`, `Rfast` and `foreach` packages on CRAN. Then:

> library(devtools)

> install_github("ericstrobl/DDR")

> library(DDR)

# Unbounded Response

> numCores <- detectCores()-1; registerDoParallel(numCores) # set up parallel computing

> cd_est = DDR(matrix(rnorm(400),200,2),rnorm(200),matrix(rnorm(20),10,2)) # run DDR on toy Gaussian data

> plot(cd_est$y,cd_est$dens[1,],type="l") # plot the conditional density estimate of the first test sample

> lines(cd_est$y,dnorm(cd_est$y),col="red") # plot ground truth in red

# Bounded Response

Recommended if you know that the response variable Y is bounded on an interval [lb,ub]

> numCores <- detectCores()-1; registerDoParallel(numCores) # set up parallel computing

> cd_est = DDR(matrix(rnorm(400),200,2),rbeta(200,0.5,0.5),matrix(rnorm(20),10,2),lb=0,ub=1) # run DDR on data sampled from a beta distribution (alpha=0.5, beta=0.5)

> plot(cd_est$y,cd_est$dens[1,],type="l") # plot the conditional density estimate of the first test sample

> lines(cd_est$y,dbeta(cd_est$y,0.5,0.5),col="red") # plot ground truth in red
