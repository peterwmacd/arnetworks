# load data
setwd('~/packages/arnetworks/real_data/')
# load package
library(devtools)
load_all()

# manufacturing network
X_man <- readRDS('data/X_man.rds')
p <- dim(X_man)[1]
n <- dim(X_man)[3]
# transitivity statistics
UV_man <- transitivity_stats(X_man)
# mean(U) ~ 1.42, 42% of entries 0
# mean(V) ~ 8, 1% of entries 0

fit_man <- estim_transitivity(X_man,verbose=TRUE)
# save fitted model
saveRDS(fit_man,file='data/fit_man.rds')

#### COMMENTS on errors ####
# initializers ab_init = c(1,1) works, other attempts threw fn=Inf error in global estimation...
# returns (a,b) = (0.193,0.177)

# with ab_init=c(0,1), no UB, fails at initialization stage, despite nearly identical input
# parameters (although a,b, get *slightly* larger, seems that this breaks globalMLE_ab)

# with ab_init=c(0,1), ab_max=4 works, returns (a,b)=(0.192,0.178)

# with ab_init=c(0,0), ab_max=4 works, returns (a,b)=(0.193,0.178)

# with ab_init=c(2,0), rough ab_max, returns (a,b)=c(0.193,0.177)

# initializers ab_init=c(1,1) ***works w/o upper bound, fails with (too large) upper bound,
# with upper box constraints, need to evaluate the boundary values?

# note that large values for a,b lead to NaN values for globalMLE_ab, why?
# eg
# globalMLE_ab(c(7,7),A1,B1,A2,B2,U,V,thetaE,etaE)
# returns NaN,
# occurs when either argument increases

# globalMLE_ab is a (negative) sum of O(p^2) ~ 1e4 terms
# each of those terms is a sum of 2 means, each with O(n-1) ~ 1e1 terms
# ie likely that Infs are occuring inside

# NaN w/o division implies there must be a sum Inf + (-Inf) somewhere...
# easy to see Inf arising: exp(800) = Inf, could occur with large b,V_ij

# ERROR:
# there is an ij such that:
# - exp(b*Vij) = Inf
# - (1- thetai*thetaj) < 0
# - exp(a*Uij) = Inf

# can fix this by either:
# - bounding a,b so that exp(bV) and exp(aU) < Inf ### currently implemented
# - bounding theta so that max(theta) < 1 (but might still get infinite fn values with large a,b)

#### Model interpretation ####

# load model
fit_man <- readRDS('data/fit_man.rds')

# a is larger (slightly, but entries of V are generally larger); not as
# extreme as the other dataset

# histogram of fitted thetas and etas
hist(fit_man$theta,20,
     main='Histogram of theta-hat',xlab='',col='blue')
abline(v=mean(fit_man$theta),lty=2)
# quite dispersed, mean around 0.5

hist(fit_man$eta,20,
     main='Histogram of eta-hat',xlab='',col='red')
abline(v=mean(fit_man$eta),lty=2)
# less dispersed, mean around 0.9

# xy plot of fitted thetas and etas
plot(fit_man$theta,fit_man$eta,
     main='Plot of theta-hat vs eta-hat',xlab='theta-hat',ylab='eta-hat',
     xlim=c(0,1.5),ylim=c(0,1.5))
abline(lm(fit_man$eta~fit_man$theta),lty=2)
# negative relationship: smallest values of eta-hat have large values of theta-hat
# implies people less likely to dissolve connections, also more likely to form new
# connections
# think: overall degree heterogeneity, rather than a split between more/less
# dynamic

# overall model summary with 4 numbers
summary_man <- rep(NA,4)
# mean of theta_i*theta_j
summary_man[1] <- mean(tcrossprod(fit_man$theta))
# mean of exp(aU)/(1 + exp(aU) + exp(bV))
summary_man[2] <- mean(exp(fit_man$a*UV_man$U)/(1 + exp(fit_man$a*UV_man$U) + exp(fit_man$b*UV_man$V)))
# mean of eta_i*eta_j
summary_man[3] <- mean(tcrossprod(fit_man$eta))
# mean of exp(bU)/(1 + exp(aU) + exp(bV))
summary_man[4] <- mean(exp(fit_man$b*UV_man$V)/(1 + exp(fit_man$a*UV_man$U) + exp(fit_man$b*UV_man$V)))

# in words:
# average about 5% chance of forming a new edge, sometimes much higher ~ 20% for largest thetas
# average about 50% chance of dissolving an existing edge, sometimes much lower ~ 15% for smallest etas

# how close are edge probabilities to form/dissolve at random?
sd(tcrossprod(fit_man$theta) * summary_man[2])/prod(summary_man[1:2]) # 0.839
sd(tcrossprod(fit_man$eta) * summary_man[4])/prod(summary_man[3:4])  # 0.258
# CoV is ~3x higher for formation probabilities vs dissolution probabilities, more
# 'signal' to predict, but still good signal for dissolution

#### analysis against metadata ####

#### GOF ####

# - formal goodness of fit?
# - validation with link prediction? better than a naive model?

# fitted flip probabilities
probs_man <- arnetworks:::model_probs(fit_man,UV_man,X_man)

# fitted residuals
resid_man <- arnetworks:::model_residuals(probs_man,X_man)

# backburner for now...

# validation for density
undir_scale <- (p^2) / (2*choose(p,2))
obs_dens <- undir_scale*apply(X,3,mean)[-1]
pred_dens <- undir_scale*apply(probs_man$gamma,3,mean)
#pred_dens <- undir_scale*apply(X,3,mean)[-n]
plot(pred_dens,obs_dens,xlim=c(.05,.15),ylim=c(.05,.15))
abline(a=0,b=1,lty=2)
abline(h=mean(obs_dens),lty=3)
abline(v=mean(pred_dens),lty=3)
# or better to plot a scatter over time
plot(1:(n-1),obs_dens,xlab='Time',ylab='Edge density')
lines(1:(n-1),pred_dens,col='red',type='p')

# could do this with other statistics? triangles? etc.

# validation for degree sequences
par(mfrow=c(3,2))
for(t in 1:(n-1)){
  # plot observed degree sequence
  plot(density(rowSums(X[,,t+1]),from=0,to=100),
       main=paste0('Obs/Pred degree distribution, time ',t+1),
       xlim=c(0,100),ylim=c(0,.1))
  # plot estimated degree sequence
  lines(density(rowSums(probs_man$gamma[,,t]),from=0,to=100),col='red')
  # null comparison against bootstrap resampling (can use t+1, 'oracle' same distn,
  # or t 'naive' prediction of degree distn)
  #lines(density(sample(rowSums(X[,,t]),replace=TRUE),from=0,to=100),col='red')
}

#### link prediction ####

# vary the size of training set
# predict entire matrix then subset to just the upper triangle

# model-based prediction: recursive step forecasting with gammas
# naive stationary prediction: use the overall mean for that node pair
# naive one-step prediction: use the previous edge value (only get a sens/spec point
# in the ROC space)

# recall pROC package functionality (see old HOSEA code)
# plotting an ROC curve, plotting additional points in the space



