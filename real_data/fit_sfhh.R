# load data
setwd('~/packages/arnetworks/real_data/')
# load package
library(devtools)
load_all()

# conference network
X_sfhh <- readRDS('data/X_sfhh.rds')
# transitivity statistics
UV_sfhh <- transitivity_stats(X_sfhh)
# mean(U) ~ 0.08, 95% of entries are zero
# mean(V) ~ 2.2, 26% of entries are zero (still quite small, perhaps shows more
# modularity?)

fit_sfhh <- estim_transitivity(X_sfhh,verbose=TRUE)
# save fitted model
saveRDS(fit_sfhh,file='data/fit_sfhh.rds')

#### Model interpretation ####

# load model
fit_sfhh <- readRDS('data/fit_sfhh.rds')

# b is larger, as are the entries of V

# histogram of fitted thetas and etas
hist(fit_sfhh$theta,20,
     main='Histogram of theta-hat',xlab='',col='blue')
abline(v=mean(fit_sfhh$theta),lty=2)
# generally small (mean 0.18) with a long right tail

hist(fit_sfhh$eta,20,
     main='Histogram of eta-hat',xlab='',col='red')
abline(v=mean(fit_sfhh$eta),lty=2)
# larger, less dispersed, mean around 1.075

# xy plot of fitted thetas and etas
plot(fit_sfhh$theta,fit_sfhh$eta,
     main='Plot of theta-hat vs eta-hat',xlab='theta-hat',ylab='eta-hat',
     xlim=c(0,1),ylim=c(.5,1.5))
abline(lm(fit_sfhh$eta~fit_sfhh$theta),lty=2)
# weak negative relationship, implies forming more edges has
# little effect on the dissolution, but if anything forming more edge makes one
# less likely to dissolve edges (again, consistent with degree heterogeneity)

# overall model summary with 4 numbers
summary_sfhh <- rep(NA,4)
# mean of theta_i*theta_j
summary_sfhh[1] <- mean(tcrossprod(fit_sfhh$theta))
# mean of exp(aU)/(1 + exp(aU) + exp(bV))
summary_sfhh[2] <- mean(exp(fit_sfhh$a*UV_sfhh$U)/(1 + exp(fit_sfhh$a*UV_sfhh$U) + exp(fit_sfhh$b*UV_sfhh$V)))
# mean of eta_i*eta_j
summary_sfhh[3] <- mean(tcrossprod(fit_sfhh$eta))
# mean of exp(bU)/(1 + exp(aU) + exp(bV))
summary_sfhh[4] <- mean(exp(fit_sfhh$b*UV_sfhh$V)/(1 + exp(fit_sfhh$a*UV_sfhh$U) + exp(fit_sfhh$b*UV_sfhh$V)))

# in words:
# average about 1% chance of forming a new edge, sometimes much higher ~ 10% for largest thetas
# average about 50% chance of dissolving an existing edge, sometimes much lower ~ 15% for smallest etas

# how close are edge probabilities to form/dissolve at random?
sd(tcrossprod(fit_sfhh$theta) * summary_sfhh[2])/prod(summary_sfhh[1:2]) # 0.908
sd(tcrossprod(fit_sfhh$eta) * summary_sfhh[4])/prod(summary_sfhh[3:4])  # 0.066
# CoV is much higher for formation probabilities vs dissolution probabilities, more
# 'signal' to predict, very little for dissolution process
# Q: will this be seen in link prediction results?

#### GOF ####
