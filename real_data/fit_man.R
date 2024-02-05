# load data
setwd('~/packages/arnetworks/real_data/')
# load packages (incl. arnetworks)
library(devtools)
library(pROC)
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

#### checking for empirical AR structure ####

rho_man <- apply(X_man,c(1,2),
                 function(v){cor(v[-1],v[-n])})

acor_obs <- mean(abs(rho_man),na.rm=TRUE)
acor_boot <- rep(0,100)
for(bb in 1:100){
  acor_boot[bb] <- mean(abs(apply(X_man[,,sample(1:n,n)],c(1,2),
                         function(v){cor(v[-1],v[-n])})),na.rm=TRUE)
}

hist(acor_boot,xlim=c(0,1))
abline(v=acor_obs,lty=2)

# appears to be a significant amount of autocorrelation in the data

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
pdf(file='fit_plots_man/theta_hist_man.pdf')
hist(fit_man$theta,20,
     main='Histogram of theta-hat',xlab='',col='blue')
abline(v=mean(fit_man$theta),lty=2)
dev.off()
# quite dispersed, mean around 0.5

pdf(file='fit_plots_man/eta_hist_man.pdf')
hist(fit_man$eta,20,
     main='Histogram of eta-hat',xlab='',col='red')
abline(v=mean(fit_man$eta),lty=2)
dev.off()
# less dispersed, mean around 0.9

# xy plot of fitted thetas and etas
#pdf(file='fit_plots_man/theta_scatter.pdf')
plot(fit_man$theta,fit_man$eta,
     main='Plot of theta-hat vs eta-hat',xlab='theta-hat',ylab='eta-hat',
     xlim=c(0,1.5),ylim=c(0,1.5))
abline(lm(fit_man$eta~fit_man$theta),lty=2)
#dev.off()
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

levels <- readRDS('data/levels_man.rds')
print(table(levels))

# same scatter plot of theta and eta parameters for different hierarchical levels
pdf(file='fit_plots_man/theta_scatter_man.pdf')
plot(fit_man$theta,fit_man$eta,
     main='Plot of theta-hat vs eta-hat, by org level',xlab='theta-hat',ylab='eta-hat',
     xlim=c(0,1.5),ylim=c(0,1.5),col=levels,cex=.8*levels)
abline(lm(fit_man$eta~fit_man$theta),lty=2)
dev.off()

# no particularly strange behavior for CEO or others high in the hierarchy

# managers tend to have higher value of theta (more likely to form new edges)
mean(fit_man$theta[levels==1])
mean(fit_man$theta[levels>1])
# also (slightly) smaller values of eta (less likely to dissolve an existing edge)
mean(fit_man$eta[levels==1])
mean(fit_man$eta[levels>1])

#### GOF ####

# - formal goodness of fit?
# - validation with link prediction? better than a naive model?

# fitted flip probabilities
probs_man <- arnetworks:::model_probs(fit_man,UV_man,X_man)
# fitted residuals
resid_man <- arnetworks:::model_residuals(probs_man,X_man)

# validation for density
undir_scale <- (p^2) / (2*choose(p,2))
obs_dens <- undir_scale*apply(X_man,3,mean)[-1]
pred_dens <- undir_scale*apply(probs_man$gamma,3,mean)
#pred_dens <- undir_scale*apply(X,3,mean)[-n]
par(mfrow=c(1,1))
plot(pred_dens,obs_dens,xlim=c(.05,.15),ylim=c(.05,.15))
abline(a=0,b=1,lty=2)
abline(h=mean(obs_dens),lty=3)
abline(v=mean(pred_dens),lty=3)
#abline(lm(obs_dens~pred_dens),lty=4)
# or better to plot a scatter over time
pdf('fit_plots_man/valid_dens_man.pdf')
plot(1:(n-1),obs_dens,xlab='Time',ylab='Edge density')
lines(1:(n-1),pred_dens,col='orange',type='p')
dev.off()

# could do this with other statistics? triangles? etc.

# validation for degree sequences
par(mfrow=c(3,2))
for(t in 1:(n-1)){
  # plot observed degree sequence
  plot(density(rowSums(X_man[,,t+1]),from=0,to=100),
       main=paste0('Obs/Pred degree distribution, time ',t+1),
       xlim=c(0,100),ylim=c(0,.1))
  # plot estimated degree sequence
  lines(density(rowSums(probs_man$gamma[,,t]),from=0,to=100),col='orange')
  # null comparison against bootstrap resampling (can use t+1, 'oracle' same distn,
  # or t 'naive' prediction of degree distn)
  #lines(density(sample(rowSums(X[,,t]),replace=TRUE),from=0,to=100),col='red')
}

#### link prediction ####
library(pROC)

# helper to take above the diagonal of a square matrix
ut <- function(M){
  c(M[upper.tri(M,diag=FALSE)])
}

# vary the size of training set, prediction horizon
# set training sets n_train=10,20,30
# set prediction horizon n_out=1,2,3

# AR model: recursive forecasts from the fitted AR network model
# Degree parameters: rank-one approximation of the time-averaged adjacency matrix
# Edge means: time-averaged adjacency matrix (more parameters than AR model)
# Previous edge: use the most recent observed edge (only get one spec/sens pair)

# set training sets n_train=10,20,30
# set prediction horizon n_out=1,2,3
n_train <- 10*1:3

# # reduced model fits
# fit_man_reduce <- list()
# length(fit_man_reduce) <- length(n_train)
# for(i in 1:3){
#   fit_man_reduce[[i]] <- estim_transitivity(X_man[,,1:n_train[i]],verbose=TRUE)
# }
# # save reduced model fits
# saveRDS(fit_man_reduce,file='data/fit_man_reduce.rds')

# load reduced model fits
fit_man_reduce <- readRDS(file='data/fit_man_reduce.rds')

# predict and plot ROCs
pdf('fit_plots_man/roc_curves_man.pdf')
par(mfrow=c(1,1))
for(i in 1:length(n_train)){
  # training set edge means or 1-dimensional approx (similar number of parameters)
  Xmean <- apply(X_man[,,1:n_train[i]],c(1,2),mean)
  eigXmean <- eigen(Xmean)
  pred_stationary <- eigXmean$values[1]*tcrossprod(eigXmean$vectors[,1])
  # model probabilities
  pred_model <- model_predict(3,fit_man_reduce[[i]],X_man[,,n_train[i]])
  for(n_out in 1:3){
    # ROC for naive stationary prediction
    temp <- roc(response=ut(X_man[,,n_train[i]+n_out]),predictor=ut(pred_stationary))
    plot(temp,col='blue',
         main=paste0('ROC, ',n_train[i],' training snapshots, horizon ',n_out))
    # ROC for model-based prediction
    temp2 <- roc(response=ut(X_man[,,n_train[i]+n_out]),predictor=ut(pred_model[,,n_out]))
    plot(temp2,col='orange',add=TRUE)
    # ROC for edge means
    temp3 <- roc(response=ut(X_man[,,n_train[i]+n_out]),predictor=ut(Xmean))
    plot(temp3,col='red',add=TRUE)
    # calculate point for naive one-step
    naive_class <- table(ut(X_man[,,n_train[i]]),ut(X_man[,,n_train[i]+n_out]))
    fpr <- naive_class[2,1]/(naive_class[2,1] + naive_class[1,1])
    tpr <- naive_class[2,2]/(naive_class[2,2] + naive_class[1,2])
    points(1-fpr,tpr,col='green',pch=15,cex=1.2)
    # legend
    legend(x=.75,y=.2,ncol=2,cex=.7,
           legend=c('AR model','Degree parameters','Edge means','Previous edge'),
           lty=c(1,1,1,NA),pch=c(NA,NA,NA,15),col=c('orange','blue','red','green'))
  }
}
dev.off()

#### link prediction V02 ####

# vary the size of training set, prediction horizon
# set training sets n_train=10,...,22, starting from time 14
# set prediction horizon n_out=1,2,3

# AR model: recursive forecasts from the fitted AR network model
# Degree parameters: rank-one approximation of the time-averaged adjacency matrix
# Edge means: time-averaged adjacency matrix (more parameters than AR model)
# Previous edge: use the most recent observed edge (only get one spec/sens pair)
# simple AR model: recursive forecasts from simple AR model

X_man_sub <- X_man[,,-(1:13)]
n_train <- 10:23

# # reduced model fits
#fit_man_reduce <- fit_man_simple <- fit_man_edge <- list()
#length(fit_man_reduce) <- length(fit_man_simple) <- length(fit_man_edge) <- length(n_train)
for(i in 1:length(n_train)){
  #fit_man_reduce[[i]] <- estim_transitivity(X_man_sub[,,1:n_train[i]],verbose=TRUE)
  fit_man_simple[[i]] <- simple_ar_fit(X_man_sub[,,1:n_train[i]])
  fit_man_edge[[i]] <- edge_ar_fit(X_man_sub[,,1:n_train[i]])
}
# save reduced model fits
#saveRDS(fit_man_reduce,file='data/fit_man_reduce.rds')
saveRDS(fit_man_simple,file='data/fit_man_simple.rds')
saveRDS(fit_man_edge,file='data/fit_man_edge.rds')

# load reduced model fits
#fit_man_reduce <- readRDS(file='data/fit_man_reduce.rds')
fit_man_simple <- readRDS(file='data/fit_man_simple.rds')
fit_man_edge <- readRDS(file='data/fit_man_edge.rds')

# predict and plot ROCs
response_combined <- list(NULL,NULL,NULL)
pred_degree_combined <- list(NULL,NULL,NULL)
pred_model_combined <- list(NULL,NULL,NULL)
pred_mean_combined <- list(NULL,NULL,NULL)
pred_simple_combined <- list(NULL,NULL,NULL)
pred_edgear_combined <- list(NULL,NULL,NULL)
pred_naive_combined <- list(matrix(0,2,2),matrix(0,2,2),matrix(0,2,2))

for(i in 1:length(n_train)){
  # training set edge means or 1-dimensional approx (similar number of parameters)
  Xmean <- apply(X_man_sub[,,1:n_train[i]],c(1,2),mean)
  eigXmean <- eigen(Xmean)
  pred_stationary <- eigXmean$values[1]*tcrossprod(eigXmean$vectors[,1])
  # model probabilities
  #pred_model <- model_predict(3,fit_man_reduce[[i]],X_man_sub[,,n_train[i]])
  # simple AR probabilities
  pred_simple <- simple_ar_predict(3,fit_man_simple[[i]],X_man_sub[,,n_train[i]])
  # edgewise AR probabilities
  pred_edgear <- edge_ar_predict(3,fit_man_edge[[i]],X_man_sub[,,n_train[i]])
  for(n_out in 1:3){
    # store response
    response_combined[[n_out]] <- c(response_combined[[n_out]],ut(X_man_sub[,,n_train[i]+n_out]))
    # store pred for naive stationary prediction
    pred_degree_combined[[n_out]] <- c(pred_degree_combined[[n_out]],ut(pred_stationary))
    # store pred for model-based prediction
    #pred_model_combined[[n_out]] <- c(pred_model_combined[[n_out]],ut(pred_model[,,n_out]))
    # store pred for mean prediction
    pred_mean_combined[[n_out]] <- c(pred_mean_combined[[n_out]],ut(Xmean))
    # store pred for simple AR prediction
    pred_simple_combined[[n_out]] <- c(pred_simple_combined[[n_out]],ut(pred_simple[,,n_out]))
    # store pred for edgewise AR prediction
    pred_edgear_combined[[n_out]] <- c(pred_edgear_combined[[n_out]],ut(pred_edgear[,,n_out]))
    # store pred for naive prediction
    pred_naive_combined[[n_out]] <- pred_naive_combined[[n_out]] + as.matrix(table(ut(X_man_sub[,,n_train[i]]),ut(X_man_sub[,,n_train[i]+n_out])))
  }
}

pdf('fit_plots_man/roc_curves_man_combined.pdf')
par(mfrow=c(1,1))
for(n_out in 1:3){
    temp <- roc(response=response_combined[[n_out]],predictor=pred_degree_combined[[n_out]])
    plot(temp,col='blue',
         main=paste0('ROC, prediction horizon ',n_out))
    # ROC for model-based prediction
    #temp2 <- roc(response=response_combined[[n_out]],predictor=pred_model_combined[[n_out]])
    #plot(temp2,col='orange',add=TRUE)
    # ROC for edge means
    temp3 <- roc(response=response_combined[[n_out]],predictor=pred_mean_combined[[n_out]])
    plot(temp3,col='red',add=TRUE)
    # ROC for simple AR model
    temp4 <- roc(response=response_combined[[n_out]],predictor=pred_simple_combined[[n_out]])
    plot(temp4,col='purple',add=TRUE)
    # ROC for edgewise AR model
    temp5 <- roc(response=response_combined[[n_out]],predictor=pred_edgear_combined[[n_out]])
    plot(temp5,col='cyan',add=TRUE)
    # calculate point for naive one-step
    naive_class <- pred_naive_combined[[n_out]]
    fpr <- naive_class[2,1]/(naive_class[2,1] + naive_class[1,1])
    tpr <- naive_class[2,2]/(naive_class[2,2] + naive_class[1,2])
    points(1-fpr,tpr,col='green',pch=15,cex=1.2)
    # legend
    legend(x=.75,y=.2,ncol=2,cex=.7,
           legend=c('AR model','Simple AR model','Edgewise AR model', 'Degree parameters','Edge means','Previous edge'),
           lty=c(1,1,1,1,1,NA),pch=c(NA,NA,NA,NA,NA,15),col=c('orange','purple','cyan','blue','red','green'))
}
dev.off()

######

# Likelihood/AIC/BIC-based model comparison

# want to do this before and after the split

# minimize:
# AIC -2ll + 2{# param}
# BIC -2ll + log(k_bic){# param}
n_bic <- (n-1)*choose(p,2)

ics <- matrix(NA,6,2)
rownames(ics) <- c('Transitivity AR model',
                   'Simple AR model',
                   'Edge-wise AR model',
                   'Edge means model',
                   'Degree parameter model',
                   'Global mean model')
colnames(ics) <- c('AIC','BIC')

# 1. Transitivity AR model
# 2p+2 parameters
gamma_tar <- model_probs(fit_man,UV_man,X_man)$gamma
ll_tar <- ar_loglike(X_man[,,-1],gamma_tar)
ics[1,1] <- 2*(2*p + 2) - 2*ll_tar
ics[1,2] <- log(n_bic)*(2*p + 2) - 2*ll_tar

# 2. Simple AR model
# 2 parameters
fit_sar <- simple_ar_fit(X_man)
gamma_sar <- fit_sar[1] + X_man[,,-n]*(1 - sum(fit_sar))
ll_sar <- ar_loglike(X_man[,,-1],gamma_sar)
ics[2,1] <- 2*2 - 2*ll_sar
ics[2,2] <- log(n_bic)*2 - 2*ll_sar

# 3. Edgewise AR model
# 2(p \choose 2) parameters
fit_ear <- edge_ar_fit(X_man)
gamma_ear <- array(rep(fit_ear$A,n-1),c(p,p,n-1)) + X_man[,,-n]*(1 - array(rep(fit_ear$A+fit_ear$B,n-1),c(p,p,n-1)))
ll_ear <- ar_loglike(X_man[,,-1],gamma_ear)
ics[3,1] <- 2*(2*choose(p,2)) - 2*ll_ear
ics[3,2] <- log(n_bic)*(2*choose(p,2)) - 2*ll_ear

# 4. Edge means model
# (p \choose 2) parameters
X_man_mean <- apply(X_man[,,-1],c(1,2),mean)
gamma_mean <- array(rep(X_man_mean,n-1),c(p,p,n-1))
ll_mean <- ar_loglike(X_man[,,-1],gamma_mean)
ics[4,1] <- 2*choose(p,2) - 2*ll_mean
ics[4,2] <- log(n_bic)*choose(p,2) - 2*ll_mean

# 5. Degree parameters/1-dim RDPG model
# p parameters
eigXmean <- eigen(X_man_mean)
X_man_deg <- pmin(pmax(eigXmean$values[1]*tcrossprod(eigXmean$vectors[,1]),0),1)
gamma_deg <- array(rep(X_man_deg,n-1),c(p,p,n-1))
ll_deg <- ar_loglike(X_man[,,-1],gamma_deg)
ics[5,1] <- 2*p - 2*ll_deg
ics[5,2] <- log(n_bic)*p - 2*ll_deg

# 6. ER/global mean
p_mean <- mean(X_man[,,-1])
gamma_global <- array(p_mean,c(p,p,n-1))
ll_global <- ar_loglike(X_man[,,-1],gamma_global)
ics[6,1] <- 2 - 2*ll_global
ics[6,2] <- log(n_bic) - 2*ll_global

saveRDS(ics,file='data/ics_man.rds')

