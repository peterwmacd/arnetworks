
<!-- README.md is generated from README.Rmd. Please edit that file -->

# arnetworks

The goal of arnetworks is to simulate, estimate and predict from
autoregressive (AR) network models with edge dependence, following the
framework specified in Fang et al. (2024+).

## Installation

You can install the development version of arnetworks from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
# devtools::install_github("peterwmacd/arnetworks")
```

## Example 1: Transitivity Model

The package provides a detailed implementation for fitting a particular
AR network model with transitivity effects (see Fang et al. (2024+),
Section X.X). This is a basic example which shows how to simulate,
estimate and predict with this model.

``` r
library(arnetworks)
# transitivity model
```

## Example 2: Persistence Model

The package also allows users to specify their own AR network models
with both local and global parameters. This is a basic example which
shows how to simulate, estimate and predict with the persistence model
(see Fang et al. (2024+), Section Y.Y).

``` r
library(arnetworks)
# Persistence Model

# Set model parameters
p = 30; n = 20
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
# NOTE: estimation implicitly conditions on the first two network snapshots

# Define edge formation and dissolution functions
fij = function(global, stats) {
  return (exp(-1 -global*stats))
}
gij <- function(global, stats) {
  return (exp(-1 -global*stats))
}
# Model fitting with estNet
fit2 <- estNet(X, fij, gij, 
                 statsAlpha, statsBeta, 
                 globInitAlpha=1, globInitBeta=1,
                 shrGPrm = 0, maxIter = 100)

# global parameter for edge formation
print(fit2$gAlphaVal)
#> [1] 0.7964095
# global parameter for edge dissolution
print(fit2$gBetaVal)
#> [1] 0.7806978
# local parameters for edge formation (quartiles)
print(quantile(fit2$xi))
#>        0%       25%       50%       75%      100% 
#> 0.4282513 0.7076076 0.8253161 1.0036155 1.3879133
# local parameters for edge formation (quartiles)
print(quantile(fit2$eta))
#>        0%       25%       50%       75%      100% 
#> 0.4111616 0.7799851 0.8964201 0.9683512 1.4531107

# Predict the next network snapshot

# Specify snapshot to predict from
Xnew = data2$X[,,n]
# Specify the sufficient statistics for prediction (as arrays)
statsAlphaNew =  1-data2$X[,,n-1,drop=FALSE]+ ( 1-data2$X[,,n-1,drop=FALSE])*( 1-data2$X[,,n-2,drop=FALSE])
statsBetaNew =  data2$X[,,n-1,drop=FALSE]+ ( data2$X[,,n-1,drop=FALSE])*( data2$X[,,n-2,drop=FALSE])
# Prediction with predictNet
pred2 = predictNet(fit2,Xnew,statsAlphaNew,statsBetaNew,fij,gij)
```
