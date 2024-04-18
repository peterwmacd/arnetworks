# load data
setwd('~/packages/arnetworks/real_data/')
# load packages
library(devtools)
library(pROC)
library(latex2exp)
# source code
load_all()
source('realdataFunctions.R')
# source('Functions.R')
# source('Estim.R')

# colorblind palete
cbp <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# conference network
X_sfhh <- readRDS('data/X_sfhh.rds')
p <- dim(X_sfhh)[1]
n <- dim(X_sfhh)[3]
# transitivity statistics
UV_sfhh <- arnetworks::statsTransitivity(X_sfhh)
# mean(U) ~ 0.08, 95% of entries are zero
# mean(V) ~ 2.2, 26% of entries are zero (still quite small, perhaps shows more
# modularity?)

fit_sfhh <- arnetworks::estTransitivity(X_sfhh,verbose=TRUE)

# save fitted model
saveRDS(fit_sfhh,file='data/fit_sfhh.rds')

#### Model interpretation ####

# load model
fit_sfhh <- readRDS('data/fit_sfhh.rds')

# estimates of a,b
print(fit_sfhh$gVal)

# means of xi, eta
print(mean(fit_sfhh$xi))
print(mean(fit_sfhh$eta))

# ranges of xi, eta
print(quantile(fit_sfhh$xi))
print(quantile(fit_sfhh$eta))

# histogram of fitted thetas and etas
pdf(file='fit_plots_sfhh/theta_hist_sfhh.pdf')
hist(fit_sfhh$xi,20,
     main=TeX('Histogram of $\\{ \\hat{\\xi}_i \\}_{i=1}^{200}$'),
     xlab='',col=cbp[3],xlim=c(0,1),cex.main=1.8,cex.lab=1.3)
abline(v=mean(fit_sfhh$xi),lty=2)
dev.off()
# generally small (mean 0.18) with a long right tail

pdf(file='fit_plots_sfhh/eta_hist_sfhh.pdf')
hist(fit_sfhh$eta,20,
     main=TeX('Histogram of $\\{ \\hat{\\eta}_i \\}_{i=1}^{200}$'),
     xlab='',col=cbp[7],xlim=c(.8,1.4),cex.main=1.8,cex.lab=1.3)
abline(v=mean(fit_sfhh$eta),lty=2)
dev.off()
# larger, less dispersed, mean around 1.08

# xy plot of fitted thetas and etas
pdf(file='fit_plots_sfhh/theta_scatter_sfhh.pdf')
par(mar=c(4,5.5,2,2))
plot(fit_sfhh$xi,fit_sfhh$eta,
     main='',xlab=TeX('$\\hat{\\xi}_i$'),ylab=TeX('$\\hat{\\eta}_i$'),
     xlim=c(0,1),ylim=c(.8,1.4),cex.lab=2)
abline(lm(fit_sfhh$eta~fit_sfhh$xi),lty=2)
dev.off()
# weak negative relationship, implies forming more edges has
# little effect on the dissolution, but if anything forming more edge makes one
# less likely to dissolve edges (again, consistent with degree heterogeneity)

#### AIC/BIC comparison ####

# minimize:
# AIC -2ll + 2{# param}
# BIC -2ll + log(k_bic){# param}

# first regime
n_bic <- (n-1)*choose(p,2)

ics <- matrix(NA,6,2)
rownames(ics) <- c('Transitivity AR model',
                    'Global AR model',
                    'Edge-wise AR model',
                    'Edge-wise mean model',
                    'Degree parameter mean model',
                    'Global mean model')
colnames(ics) <- c('AIC','BIC')

# 1. Transitivity AR model
# 2p+2 parameters
gamma_tar <- model_probs(fit_sfhh,UV_sfhh,X_sfhh)$gamma
ll_tar <- ar_loglike(X_sfhh[,,-1],gamma_tar)
ics[1,1] <- 2*(2*p + 2) - 2*ll_tar
ics[1,2] <- log(n_bic)*(2*p + 2) - 2*ll_tar

# 2. Simple AR model
# 2 parameters
fit_sar <- simple_ar_fit(X_sfhh)
gamma_sar <- fit_sar[1] + X_sfhh[,,-n]*(1 - sum(fit_sar))
ll_sar <- ar_loglike(X_sfhh[,,-1],gamma_sar)
ics[2,1] <- 2*2 - 2*ll_sar
ics[2,2] <- log(n_bic)*2 - 2*ll_sar

# 3. Edgewise AR model
# 2(p \choose 2) parameters
fit_ear <- edge_ar_fit(X_sfhh)
gamma_ear <- array(rep(fit_ear$A,n-1),c(p,p,n-1)) + X_sfhh[,,-n]*(1 - array(rep(fit_ear$A+fit_ear$B,n-1),c(p,p,n-1)))
ll_ear <- ar_loglike(X_sfhh[,,-1],gamma_ear)
ics[3,1] <- 2*(2*choose(p,2)) - 2*ll_ear
ics[3,2] <- log(n_bic)*(2*choose(p,2)) - 2*ll_ear

# 4. Edge means model
# (p \choose 2) parameters
X_sfhh_mean <- apply(X_sfhh[,,-1],c(1,2),mean)
gamma_mean <- array(rep(X_sfhh_mean,n-1),c(p,p,n-1))
ll_mean <- ar_loglike(X_sfhh[,,-1],gamma_mean)
ics[4,1] <- 2*choose(p,2) - 2*ll_mean
ics[4,2] <- log(n_bic)*choose(p,2) - 2*ll_mean

# 5. Degree parameters/1-dim RDPG model
# p parameters
eigXmean <- eigen(X_sfhh_mean)
X_sfhh_deg <- pmin(pmax(eigXmean$values[1]*tcrossprod(eigXmean$vectors[,1]),0),1)
gamma_deg <- array(rep(X_sfhh_deg,n-1),c(p,p,n-1))
ll_deg <- ar_loglike(X_sfhh[,,-1],gamma_deg)
ics[5,1] <- 2*p - 2*ll_deg
ics[5,2] <- log(n_bic)*p - 2*ll_deg

# 6. ER/global mean
p_mean <- mean(X_sfhh[,,-1])
gamma_global <- array(p_mean,c(p,p,n-1))
ll_global <- ar_loglike(X_sfhh[,,-1],gamma_global)
ics[6,1] <- 2 - 2*ll_global
ics[6,2] <- log(n_bic) - 2*ll_global

# reorder for presentation
ics <- ics[c(1,2,3,5,6,4),]

saveRDS(ics,file='data/ics_sfhh.rds')

#### link prediction ####

# # vary the size of training set, prediction horizon
# # set training sets n_train=10,20,30
# # set prediction horizon n_out=1,2,3
#
# # AR model: recursive forecasts from the fitted AR network model
# # Degree parameters: rank-one approximation of the time-averaged adjacency matrix
# # Edge means: time-averaged adjacency matrix (more parameters than AR model)
# # Previous edge: use the most recent observed edge (only get one spec/sens pair)
#
# # set training sets n_train=6,12,18
# # set prediction horizon n_out=1,2,3
# n_train <- 6*1:3
#
# # reduced model fits
# # fit_sfhh_reduce <- list()
# # length(fit_sfhh_reduce) <- length(n_train)
# # for(i in 1:3){
# #   fit_sfhh_reduce[[i]] <- estim_transitivity(X_sfhh[,,1:n_train[i]],verbose=TRUE)
# # }
# # # save reduced model fits
# # saveRDS(fit_sfhh_reduce,file='data/fit_sfhh_reduce.rds')
#
# # load reduced model fits
# fit_sfhh_reduce <- readRDS(file='data/fit_sfhh_reduce.rds')
#
# # predict and plot ROCs
# pdf('fit_plots_sfhh/roc_curves_sfhh.pdf')
# par(mfrow=c(1,1))
# for(i in 1:length(n_train)){
#   # training set edge means or 1-dimensional approx (similar number of parameters)
#   Xmean <- apply(X_sfhh[,,1:n_train[i]],c(1,2),mean)
#   eigXmean <- eigen(Xmean)
#   pred_stationary <- eigXmean$values[1]*tcrossprod(eigXmean$vectors[,1])
#   # model probabilities
#   pred_model <- model_predict(3,fit_sfhh_reduce[[i]],X_sfhh[,,n_train[i]])
#   for(n_out in 1:3){
#     # ROC for naive stationary prediction
#     temp <- roc(response=ut(X_sfhh[,,n_train[i]+n_out]),predictor=ut(pred_stationary))
#     plot(temp,col='blue',
#          main=paste0('ROC, ',n_train[i],' training snapshots, horizon ',n_out))
#     # ROC for model-based prediction
#     temp2 <- roc(response=ut(X_sfhh[,,n_train[i]+n_out]),predictor=ut(pred_model[,,n_out]))
#     plot(temp2,col='orange',add=TRUE)
#     # ROC for edge means
#     temp3 <- roc(response=ut(X_sfhh[,,n_train[i]+n_out]),predictor=ut(Xmean))
#     plot(temp3,col='red',add=TRUE)
#     # calculate point for naive one-step
#     naive_class <- table(ut(X_sfhh[,,n_train[i]]),ut(X_sfhh[,,n_train[i]+n_out]))
#     fpr <- naive_class[2,1]/(naive_class[2,1] + naive_class[1,1])
#     tpr <- naive_class[2,2]/(naive_class[2,2] + naive_class[1,2])
#     points(1-fpr,tpr,col='green',pch=15,cex=1.2)
#     # legend
#     legend(x=.75,y=.2,ncol=2,cex=.7,
#            legend=c('AR model','Degree parameters','Edge means','Previous edge'),
#            lty=c(1,1,1,NA),pch=c(NA,NA,NA,15),col=c('orange','blue','red','green'))
#   }
# }
# dev.off()
