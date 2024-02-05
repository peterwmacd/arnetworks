# load data
setwd('~/packages/arnetworks/real_data/')
# load package
library(devtools)
library(pROC)
load_all()

# conference network
X_sfhh <- readRDS('data/X_sfhh.rds')
p <- dim(X_sfhh)[1]
n <- dim(X_sfhh)[3]
# transitivity statistics
UV_sfhh <- transitivity_stats(X_sfhh)
# mean(U) ~ 0.08, 95% of entries are zero
# mean(V) ~ 2.2, 26% of entries are zero (still quite small, perhaps shows more
# modularity?)

# fit_sfhh <- estim_transitivity(X_sfhh,verbose=TRUE)
# # save fitted model
# saveRDS(fit_sfhh,file='data/fit_sfhh.rds')

#### AR structure checking ###

rho_sfhh <- apply(X_sfhh,c(1,2),
                 function(v){cor(v[-1],v[-n])})

acor_obs <- mean(abs(rho_sfhh),na.rm=TRUE)
acor_boot <- rep(0,100)
for(bb in 1:100){
  acor_boot[bb] <- mean(abs(apply(X_sfhh[,,sample(1:n,n)],c(1,2),
                                  function(v){cor(v[-1],v[-n])})),na.rm=TRUE)
}

hist(acor_boot,xlim=c(0,1))
abline(v=acor_obs,lty=2)


#### Model interpretation ####

# load model
fit_sfhh <- readRDS('data/fit_sfhh.rds')

# b is larger, as are the entries of V

# histogram of fitted thetas and etas
pdf(file='fit_plots_sfhh/theta_hist_sfhh.pdf')
hist(fit_sfhh$theta,20,
     main='Histogram of theta-hat',xlab='',col='blue')
abline(v=mean(fit_sfhh$theta),lty=2)
dev.off()
# generally small (mean 0.18) with a long right tail

pdf(file='fit_plots_sfhh/eta_hist_sfhh.pdf')
hist(fit_sfhh$eta,20,
     main='Histogram of eta-hat',xlab='',col='red')
abline(v=mean(fit_sfhh$eta),lty=2)
dev.off()
# larger, less dispersed, mean around 1.075

# xy plot of fitted thetas and etas
pdf(file='fit_plots_sfhh/theta_scatter_sfhh.pdf')
plot(fit_sfhh$theta,fit_sfhh$eta,
     main='Plot of theta-hat vs eta-hat',xlab='theta-hat',ylab='eta-hat',
     xlim=c(0,1),ylim=c(.5,1.5))
abline(lm(fit_sfhh$eta~fit_sfhh$theta),lty=2)
dev.off()
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

# fitted flip probabilities
probs_sfhh <- arnetworks:::model_probs(fit_sfhh,UV_sfhh,X_sfhh)
# fitted residuals
resid_sfhh <- arnetworks:::model_residuals(probs_sfhh,X_sfhh)

# validation for density
undir_scale <- (p^2) / (2*choose(p,2))
obs_dens <- undir_scale*apply(X_sfhh,3,mean)[-1]
pred_dens <- undir_scale*apply(probs_sfhh$gamma,3,mean)
#pred_dens <- undir_scale*apply(X,3,mean)[-n]
par(mfrow=c(1,1))
plot(pred_dens,obs_dens,xlim=c(0,.05),ylim=c(0,.05))
abline(a=0,b=1,lty=2)
abline(h=mean(obs_dens),lty=3)
abline(v=mean(pred_dens),lty=3)
# or better to plot a scatter over time
pdf(file='fit_plots_sfhh/valid_dens_sfhh.pdf')
plot(1:(n-1),obs_dens,xlab='Time',ylab='Edge density')
lines(1:(n-1),pred_dens,col='orange',type='p')
dev.off()

#### link prediction ####

# vary the size of training set, prediction horizon
# set training sets n_train=10,20,30
# set prediction horizon n_out=1,2,3

# AR model: recursive forecasts from the fitted AR network model
# Degree parameters: rank-one approximation of the time-averaged adjacency matrix
# Edge means: time-averaged adjacency matrix (more parameters than AR model)
# Previous edge: use the most recent observed edge (only get one spec/sens pair)

# set training sets n_train=6,12,18
# set prediction horizon n_out=1,2,3
n_train <- 6*1:3

# reduced model fits
# fit_sfhh_reduce <- list()
# length(fit_sfhh_reduce) <- length(n_train)
# for(i in 1:3){
#   fit_sfhh_reduce[[i]] <- estim_transitivity(X_sfhh[,,1:n_train[i]],verbose=TRUE)
# }
# # save reduced model fits
# saveRDS(fit_sfhh_reduce,file='data/fit_sfhh_reduce.rds')

# load reduced model fits
fit_sfhh_reduce <- readRDS(file='data/fit_sfhh_reduce.rds')

# predict and plot ROCs
pdf('fit_plots_sfhh/roc_curves_sfhh.pdf')
par(mfrow=c(1,1))
for(i in 1:length(n_train)){
  # training set edge means or 1-dimensional approx (similar number of parameters)
  Xmean <- apply(X_sfhh[,,1:n_train[i]],c(1,2),mean)
  eigXmean <- eigen(Xmean)
  pred_stationary <- eigXmean$values[1]*tcrossprod(eigXmean$vectors[,1])
  # model probabilities
  pred_model <- model_predict(3,fit_sfhh_reduce[[i]],X_sfhh[,,n_train[i]])
  for(n_out in 1:3){
    # ROC for naive stationary prediction
    temp <- roc(response=ut(X_sfhh[,,n_train[i]+n_out]),predictor=ut(pred_stationary))
    plot(temp,col='blue',
         main=paste0('ROC, ',n_train[i],' training snapshots, horizon ',n_out))
    # ROC for model-based prediction
    temp2 <- roc(response=ut(X_sfhh[,,n_train[i]+n_out]),predictor=ut(pred_model[,,n_out]))
    plot(temp2,col='orange',add=TRUE)
    # ROC for edge means
    temp3 <- roc(response=ut(X_sfhh[,,n_train[i]+n_out]),predictor=ut(Xmean))
    plot(temp3,col='red',add=TRUE)
    # calculate point for naive one-step
    naive_class <- table(ut(X_sfhh[,,n_train[i]]),ut(X_sfhh[,,n_train[i]+n_out]))
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
