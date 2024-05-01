
<!-- README.md is generated from README.Rmd. Please edit that file -->

# arnetworks

The goal of arnetworks is to simulate, estimate and predict from
autoregressive (AR) network models with edge dependence, following the
framework specified in Chang et al. (2024+).

## Installation

You can install the development version of arnetworks from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
# devtools::install_github("peterwmacd/arnetworks")
```

## Example 1: Transitivity Model

The package provides a detailed implementation for fitting a particular
AR network model with transitivity effects (see **Chang et al. (2024+),
Section 3.3**). This is a basic example which shows how to simulate,
estimate and predict with this model.

``` r
library(arnetworks)
# Transitivity model

p = 100; n = 50
xi = rep(0.7, p); eta = rep(0.8, p)
a = 30; b = 15

# Simulate data using simulateTransitivity function
data1 = simulateTransitivity(p, n, xi, eta, a, b)
X = data1$X
U = data1$U
V = data1$V

# in addition to the data, simulateTransitivity returns the sufficient statistics, the
# scaled number of common and uncommon neighbours for each edge as U and V respectively

fit1 = estTransitivity(X, U, V, rSeqGlob=c(50,10), rSeqLoc=c(0.5,0.1))

# the model optimization requires tuning parameters rSeqGlob and rSeqLoc 
# (see Chang et al. (2024+), Section 4.3). These parameters control the radius of search for the
# estimator refinement. Default values are provided if not specified (see the documenation for
# estTransitivity)

# global parameters for edge formation and dissoluton
print(fit1$gVal)
#> [1] 28.87053 15.27376
# local parameters for edge formation (quartiles)
print(quantile(fit1$xi))
#>        0%       25%       50%       75%      100% 
#> 0.6239415 0.7053004 0.7267167 0.7647692 0.8599142
# local parameters for edge formation (quartiles)
print(quantile(fit1$eta))
#>        0%       25%       50%       75%      100% 
#> 0.6621665 0.7676051 0.7907868 0.8242973 0.8891144

# note that the same global parameters are shared by the edge formation and dissolution models

# Predict the next two network snapshots

# Specify snapshot to predict from
Xnew = data1$X[,,n]
# Prediction with predictNet
pred1 = predictTransitivity(fit1,Xnew,nStep=2)

# predictTransitivity predicts from the estimated model parameters recursively; that is, after
# predicting the n+1 snapshot, it plugs those edge probabilities into the model again to predict the
# n+2 snapshot.
# The resulting predictions are returned as a p x p x nStep array
```

## Example 2: Persistence Model

The package also allows users to **specify their own AR network models**
with both local and global parameters. This is a basic example which
shows how to simulate, estimate and predict with the **persistence
model** (see Chang et al. (2024+), Section 3.2). For an additional
example on how to estimate and predict with the **transitivity model**,
please see the documenation for `estNet`.

``` r
library(arnetworks)
# Persistence model

# Set model parameters
p = 100; n = 50
xi = runif(p, 0.5, 0.9)
eta = runif(p, 0.5, 0.9)
a = 0.5
b = 0.5
# Simulation with simulatePersistence
data2 = simulatePersistence(p, n, xi, eta, a, b)

# Specify model components for fitting with EstNet

# Specify data and sufficient statistics
X = data2$X
statsAlpha =  1-X[,,4:n-2]+ ( 1-X[,,4:n-2])*( 1-X[,,4:n-3])
statsBeta =  X[,,4:n-2]+ ( X[,,4:n-2])*( X[,,4:n-3])
X = X[,, 3: n]
# NOTE: estimation for the persistence model implicitly conditions on the first two network snapshots,
# so the data used for fitting has n-2 network snapshots

# In general statsAlpha and statsBeta are p x p x d x n' arrays, where d is the dimension of the 
# sufficient statistic for each edge, and n' is the number of snapshots used for fitting. 
# In this case d = 1 and n' = n-2 so the array have only 3 dimensions

# Define edge formation and dissolution functions
fij = function(global, stats) {
  return (exp(-1 -global*stats))
}
gij <- function(global, stats) {
  return (exp(-1 -global*stats))
}

# These functions take the global parameters and sufficient statistics as input and return
# the global part of the edge growth (fij) and dissolution (gij) probabilities.
# In this case there is one global parameter governing growth, and one governing dissolution, so
# both functions take a scalar and an array valued parameter, and return p x p x n' arrays, where
# n' is the number of snapshots used for fitting.

# Model fitting with estNet
fit2 <- estNet(X, fij, gij, 
                 statsAlpha, statsBeta, 
                 globInitAlpha=1, globInitBeta=1,
                 shrGPrm = 0, maxIter = 100)

# global parameter for edge formation
print(fit2$gAlphaVal)
#> [1] 0.5321285
# global parameter for edge dissolution
print(fit2$gBetaVal)
#> [1] 0.5177779
# local parameters for edge formation (quartiles)
print(quantile(fit2$xi))
#>        0%       25%       50%       75%      100% 
#> 0.4166115 0.5810328 0.6607679 0.7851298 1.0219436
# local parameters for edge formation (quartiles)
print(quantile(fit2$eta))
#>        0%       25%       50%       75%      100% 
#> 0.4453189 0.5823875 0.6773083 0.7655460 0.9533456

# Predict the next network snapshot

# Specify snapshot to predict from
Xnew = data2$X[,,n]
# Specify the sufficient statistics for prediction (as arrays)
statsAlphaNew =  1-data2$X[,,n-1,drop=FALSE]+ ( 1-data2$X[,,n-1,drop=FALSE])*( 1-data2$X[,,n-2,drop=FALSE])
statsBetaNew =  data2$X[,,n-1,drop=FALSE]+ ( data2$X[,,n-1,drop=FALSE])*( data2$X[,,n-2,drop=FALSE])
# Prediction with predictNet
pred2 = predictNet(fit2,Xnew,statsAlphaNew,statsBetaNew,fij,gij)

# predictNet uses the same objects as estNet to predict the next network snapshot.
# Note that Xnew (the final snapshot of the observed network data) is a p x p matrix, while 
# statsAlphaNew and statsBetaNew are p x p x 1 arrays
```
