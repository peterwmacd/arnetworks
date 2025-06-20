# refitting with btergm

# load data
setwd('~/packages/arnetworks/real_data/')
# load packages (incl. arnetworks)
library(devtools)
library(pROC)
library(latex2exp)
library(btergm)

# source code from arnetworks package
load_all()
source('realdataFunctions.R')
#source('Functions.R')
#source('Estim.R')

# black and white indicator for plots
bw <- FALSE

# manufacturing network
X_man <- readRDS('data/X_man.rds')
p <- dim(X_man)[1]
n <- dim(X_man)[3]

# calculate U/V stats
# splitting into two regimes for fitting
# first regime snapshots 1 to 13
X_man1 <- X_man[,,1:13]
UV_man1 <- arnetworks::statsTransitivity(X_man1) # scaled by p-2
n1 <- dim(X_man1)[3]

# second regime snapshots 14 to 39
X_man2 <- X_man[,,14:39]
UV_man2 <- arnetworks::statsTransitivity(X_man2) # scaled by p-2
n2 <- dim(X_man2)[3]

# tergm for first regime

# data
net1 <- list()
U1 <- list()
V1 <- list()
for (ii in 1:n1) { 
  # store networks
  if(ii > 1){
    nw <- network::network(X_man1[,,ii-1],directed=FALSE)  # create network object
    net1[[(ii-1)]] <- nw          # add network to the list
  }
  if(ii < n1){
    U1[[ii]] <- UV_man1$U[,,ii]
    V1[[ii]] <- UV_man1$V[,,ii]
  }
}

# fit with just U (common neighbours) effect
fit1 <- mtergm(net1 ~ edges + edgecov(U1)) 
summary(fit1)

# tergm for second regime

# data
net2 <- list()
U2 <- list()
V2 <- list()
for (ii in 1:n2) { 
  # store networks
  if(ii > 1){
    nw <- network::network(X_man2[,,ii-1],directed=FALSE)  # create network object
    net2[[(ii-1)]] <- nw          # add network to the list
  }
  if(ii < n2){
    U2[[ii]] <- UV_man2$U[,,ii]
    V2[[ii]] <- UV_man2$V[,,ii]
  }
}

# fit with just U (common neighbours) effect
fit2 <- mtergm(net2 ~ edges + edgecov(U2))
summary(fit2)

# similarly find strong transitivity effect of previous common neighbours on edge presence, although
# this TERGM does not isolate the effects on formation and dissolution conditional on the previous edge 
# status, nor can it incorporate node-specific popularity effects (xis and etas in our model)


