# AR transitivity fitting and analysis for manufacturing email data

# libraries
library(arnetworks)
library(latex2exp)
library(pROC)

# helpers
source('real_data/realdataFunctions.R')

# colorblind palete
cbp <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#### (1) data preprocessing ####

# manufacturing network
X_man <- readRDS('real_data/data/X_man.rds')
p <- dim(X_man)[1]
n <- dim(X_man)[3]

# first regime snapshots 1 to 13
X_man1 <- X_man[,,1:13]
UV_man1 <- arnetworks::statsTransitivity(X_man1)
n1 <- dim(X_man1)[3]

# second regime snapshots 14 to 39
X_man2 <- X_man[,,14:39]
UV_man2 <- arnetworks::statsTransitivity(X_man2)
n2 <- dim(X_man2)[3]

#### (2) main model fitting ####

# fit AR transitivity to first period
fit_man1 <- arnetworks::estTransitivity(X_man1,
                                        tauSeq_a = 0.3, tauSeq_b = 0.3,
                                        tauSeq_xi = 0.05, tauSeq_eta = 0.05,
                                        verbose=TRUE,doInference=FALSE)
# NOTE: runs in ~5 minutes
# save fitted model
saveRDS(fit_man1,file='real_data/data/fit_man1.rds')

# fit AR transitivity to second period
fit_man2 <- arnetworks::estTransitivity(X_man2,
                                        tauSeq_a = 0.3, tauSeq_b = 0.3,
                                        tauSeq_xi = 0.05, tauSeq_eta = 0.05,
                                        verbose=TRUE,doInference=FALSE)
# NOTE: runs in ~6 minutes
# save fitted model
saveRDS(fit_man2,file='real_data/data/fit_man2.rds')

#### (3) model interpretation ####

# load models
fit_man1 <- readRDS('real_data/data/fit_man1.rds')
fit_man2 <- readRDS('real_data/data/fit_man2.rds')

# first regime

# global parameters
print(c(fit_man1$gVal))
# a is larger (slightly, but entries of V are generally larger); not as
# extreme as the other dataset

# local parameter summaries
print(mean(fit_man1$xi))
print(mean(fit_man1$eta))

# second regime

# global parameters
print(c(fit_man2$gVal))
# a is larger (slightly, but entries of V are generally larger); not as
# extreme as the other dataset

# local parameter summaries
print(mean(fit_man2$xi))
print(mean(fit_man2$eta))

# load metadata (organizational levels)
levels <- readRDS('real_data/data/levels_man.rds')
print(table(levels))

# first regime

# produce degree_scatter_man1.pdf (Figure 1, left panel)
# scatter plot of theta and eta parameters for different hierarchical levels
pdf(file='real_data/fit_plots_man/degree_scatter_man1.pdf')
par(mar=c(4,5.5,4,4))
plot(fit_man1$xi,fit_man1$eta,
     main='Period 1',xlab=TeX('$\\hat{\\xi}_i'),ylab=TeX('$\\hat{\\eta}_i'),
     xlim=c(0,2.6),ylim=c(0,1.5),col=levels,cex=.8*levels,cex.lab=1.8,cex.main=1.5)
abline(lm(fit_man1$eta~fit_man1$xi),lty=2)
dev.off()

# means of xi-hat for managers and non-managers
print(mean(fit_man1$xi[levels==1]))
print(mean(fit_man1$xi[levels>1]))

# second regime

# produce degree_scatter_man1.pdf (Figure 1, right panel)
# scatter plot of theta and eta parameters for different hierarchical levels
pdf(file='real_data/fit_plots_man/degree_scatter_man2.pdf')
par(mar=c(4,5.5,4,4))
plot(fit_man2$xi,fit_man2$eta,
     main='Period 2',xlab=TeX('$\\hat{\\xi}_i'),ylab=TeX('$\\hat{\\eta}_i'),
     xlim=c(0,2.6),ylim=c(0,1.5),col=levels,cex=.8*levels,cex.lab=1.8,cex.main=1.5)
abline(lm(fit_man1$eta~fit_man1$xi),lty=2)
dev.off()

# means of xi-hat for managers and non-managers
print(mean(fit_man2$xi[levels==1]))
print(mean(fit_man2$xi[levels>1]))

#### (4) link prediction: fitting ####

# only for second regime
# set training sets n_train=10,...,22, starting from time 14
# set prediction horizon n_out=1,2,3

# AR model: recursive forecasts from the fitted AR network model
# Degree parameters: rank-one approximation of the time-averaged adjacency matrix
# Edge means: time-averaged adjacency matrix (more parameters than AR model)
# Previous edge: use the most recent observed edge (only get one spec/sens pair)
# simple AR model: recursive forecasts from simple AR model

# subset to second regime
n_train <- 10:23

# reduced model fits
fit_man_reduce <- fit_man_simple <- fit_man_edge <- fit_man_edgemod <- list()
length(fit_man_reduce) <- length(fit_man_simple) <- length(fit_man_edge) <- length(n_train)
for(i in 1:length(n_train)){
  fit_man_reduce[[i]] <- arnetworks::estTransitivity(X_man2[,,1:n_train[i]],
                                                     tauSeq_a = 0.3, tauSeq_b = 0.3,
                                                     tauSeq_xi = 0.05, tauSeq_eta = 0.05,
                                                     verbose=TRUE,doInference=FALSE)
  fit_man_simple[[i]] <- simple_ar_fit(X_man2[,,1:n_train[i]])
  fit_man_edge[[i]] <- edge_ar_fit(X_man2[,,1:n_train[i]])
  print(paste0('models ',i,' of ',length(n_train)))
}
# NOTE: runs in ~60 minutes

# temporary code to update just one list
# fit_man_edgemod <- list()
# for(i in 1:length(n_train)){
#   fit_man_edgemod[[i]] <- edge_ar_fit_modified(X_man2[,,1:n_train[i]])
# }

# save reduced model fits
saveRDS(fit_man_reduce,file='real_data/data/fit_man2_reduce.rds')
saveRDS(fit_man_simple,file='real_data/data/fit_man2_simple.rds')
saveRDS(fit_man_edge,file='real_data/data/fit_man2_edgewise.rds')

#### (5) link prediction: results ####

# load reduced model fits
fit_man_reduce <- readRDS(file='real_data/data/fit_man2_reduce.rds')
fit_man_simple <- readRDS(file='real_data/data/fit_man2_simple.rds')
fit_man_edge <- readRDS(file='real_data/data/fit_man2_edgewise.rds')

# predict and plot ROCs
response_combined <- list(NULL,NULL,NULL)
pred_degree_combined <- list(NULL,NULL,NULL)
pred_model_combined <- list(NULL,NULL,NULL)
pred_mean_combined <- list(NULL,NULL,NULL)
pred_simple_combined <- list(NULL,NULL,NULL)
pred_edgear_combined <- list(NULL,NULL,NULL)
pred_naive_combined <- list(matrix(0,2,2),matrix(0,2,2),matrix(0,2,2))

# subset to regime 2
n_train <- 10:23
# loop over training set sizes, compute predictions
for(i in 1:length(n_train)){
  # training set edge means or 1-dimensional approx (similar number of parameters)
  Xmean <- apply(X_man2[,,1:n_train[i]],c(1,2),mean)
  eigXmean <- eigen(Xmean)
  pred_stationary <- eigXmean$values[1]*tcrossprod(eigXmean$vectors[,1])
  # model probabilities
  pred_model <- arnetworks::predictTransitivity(fit_man_reduce[[i]],X_man2[,,n_train[i]],nStep=3)
  #pred_model <- model_predict(3,fit_man_reduce[[i]],X_man2[,,n_train[i]])
  # simple AR probabilities
  pred_simple <- simple_ar_predict(3,fit_man_simple[[i]],X_man2[,,n_train[i]])
  # edgewise AR probabilities
  #pred_edgear <- edge_ar_predict(3,fit_man_edge[[i]],X_man2[,,n_train[i]])
  pred_edgear <- edge_ar_predict(3,fit_man_edge[[i]],X_man2[,,n_train[i]])
  for(n_out in 1:3){
    # store response
    response_combined[[n_out]] <- c(response_combined[[n_out]],ut(X_man2[,,n_train[i]+n_out]))
    # store pred for naive stationary prediction
    pred_degree_combined[[n_out]] <- c(pred_degree_combined[[n_out]],ut(pred_stationary))
    # store pred for model-based prediction
    pred_model_combined[[n_out]] <- c(pred_model_combined[[n_out]],ut(pred_model[,,n_out]))
    # store pred for mean prediction
    pred_mean_combined[[n_out]] <- c(pred_mean_combined[[n_out]],ut(Xmean))
    # store pred for simple AR prediction
    pred_simple_combined[[n_out]] <- c(pred_simple_combined[[n_out]],ut(pred_simple[,,n_out]))
    # store pred for edgewise AR prediction
    pred_edgear_combined[[n_out]] <- c(pred_edgear_combined[[n_out]],ut(pred_edgear[,,n_out]))
    # store pred for naive prediction
    pred_naive_combined[[n_out]] <- pred_naive_combined[[n_out]] + as.matrix(table(ut(X_man2[,,n_train[i]]),ut(X_man2[,,n_train[i]+n_out])))
  }
}

# produce roc_curves_man2_hor1-3.pdf, (Figure 2, all panels)
for(n_out in 1:3){
  pdf(paste0('real_data/fit_plots_man/roc_curves_man2_hor',n_out,'.pdf'),width=6,height=6)
  par(mfrow=c(1,1))
  temp <- roc(response=response_combined[[n_out]],predictor=pred_degree_combined[[n_out]])
  plot(temp,col=cbp[3],
       main=TeX(paste0('ROC, $n_{{step}}=',n_out,'$')),
       cex.lab=1.5,cex.main=1.5)
  # ROC for model-based prediction
  temp2 <- roc(response=response_combined[[n_out]],predictor=pred_model_combined[[n_out]])
  plot(temp2,col=cbp[2],add=TRUE)
  # ROC for edge means
  temp3 <- roc(response=response_combined[[n_out]],predictor=pred_mean_combined[[n_out]])
  plot(temp3,col=cbp[7],add=TRUE)
  # ROC for simple AR model
  temp4 <- roc(response=response_combined[[n_out]],predictor=pred_simple_combined[[n_out]])
  plot(temp4,col=cbp[4],add=TRUE)
  # ROC for edgewise AR model
  temp5 <- roc(response=response_combined[[n_out]],predictor=pred_edgear_combined[[n_out]])
  plot(temp5,col=cbp[6],add=TRUE)
  # calculate point for naive one-step
  naive_class <- pred_naive_combined[[n_out]]
  fpr <- naive_class[2,1]/(naive_class[2,1] + naive_class[1,1])
  tpr <- naive_class[2,2]/(naive_class[2,2] + naive_class[1,2])
  points(1-fpr,tpr,col=cbp[1],pch=15,cex=1.2)
  # legend
  if(n_out==1){
    legend(x=.585,y=.42,ncol=1,cex=1.2,lwd=2,
           legend=c('Transitivity model','Global AR model','Edgewise AR model', 'Degree mean model','Edgewise mean model','Previous edge'),
           lty=c(1,1,1,1,1,NA),pch=c(NA,NA,NA,NA,NA,15),col=cbp[c(2,4,6,3,7,1)])
  }
  dev.off()
}

#### (6) AIC/BIC comparison ####

# minimize:
# AIC -2ll + 2{# param}
# BIC -2ll + log(k_bic){# param}

# first regime
n1_bic <- (n1-1)*choose(p,2)

# initialize and populate ics1 (Table 2, columns 2-3)
ics1 <- matrix(NA,6,2)
rownames(ics1) <- c('Transitivity AR model',
                    'Global AR model',
                    'Edge-wise AR model',
                    'Edge-wise mean model',
                    'Degree parameter mean model',
                    'Global mean model')
colnames(ics1) <- c('AIC','BIC')

# 1. Transitivity AR model
# 2p+2 parameters
gamma_tar <- model_probs(fit_man1,UV_man1,X_man1)$gamma
ll_tar <- ar_loglike(X_man1[,,-1],gamma_tar)
ics1[1,1] <- 2*(2*p + 2) - 2*ll_tar
ics1[1,2] <- log(n1_bic)*(2*p + 2) - 2*ll_tar

# 2. Simple AR model
# 2 parameters
fit_sar <- simple_ar_fit(X_man1)
gamma_sar <- fit_sar[1] + X_man1[,,-n1]*(1 - sum(fit_sar))
ll_sar <- ar_loglike(X_man1[,,-1],gamma_sar)
ics1[2,1] <- 2*2 - 2*ll_sar
ics1[2,2] <- log(n1_bic)*2 - 2*ll_sar

# 3. Edgewise AR model
# 2(p \choose 2) parameters
fit_ear <- edge_ar_fit(X_man1)
gamma_ear <- array(rep(fit_ear$A,n1-1),c(p,p,n1-1)) + X_man1[,,-n1]*(1 - array(rep(fit_ear$A+fit_ear$B,n1-1),c(p,p,n1-1)))
ll_ear <- ar_loglike(X_man1[,,-1],gamma_ear)
ics1[3,1] <- 2*(2*choose(p,2)) - 2*ll_ear
ics1[3,2] <- log(n1_bic)*(2*choose(p,2)) - 2*ll_ear

# 4. Edge means model
# (p \choose 2) parameters
X_man_mean <- apply(X_man1[,,-1],c(1,2),mean)
gamma_mean <- array(rep(X_man_mean,n1-1),c(p,p,n1-1))
ll_mean <- ar_loglike(X_man1[,,-1],gamma_mean)
ics1[4,1] <- 2*choose(p,2) - 2*ll_mean
ics1[4,2] <- log(n1_bic)*choose(p,2) - 2*ll_mean

# 5. Degree parameters/1-dim RDPG model
# p parameters
eigXmean <- eigen(X_man_mean)
X_man_deg <- pmin(pmax(eigXmean$values[1]*tcrossprod(eigXmean$vectors[,1]),0),1)
gamma_deg <- array(rep(X_man_deg,n1-1),c(p,p,n1-1))
ll_deg <- ar_loglike(X_man1[,,-1],gamma_deg)
ics1[5,1] <- 2*p - 2*ll_deg
ics1[5,2] <- log(n1_bic)*p - 2*ll_deg

# 6. ER/global mean
p_mean <- mean(X_man1[,,-1])
gamma_global <- array(p_mean,c(p,p,n1-1))
ll_global <- ar_loglike(X_man1[,,-1],gamma_global)
ics1[6,1] <- 2 - 2*ll_global
ics1[6,2] <- log(n1_bic) - 2*ll_global

# reorder for presentation
ics1 <- ics1[c(1,2,3,5,6,4),]

# save ics1 (Table 2, columns 2-3)
saveRDS(ics1,file='real_data/data/ics_man1.rds')

# second regime
n2_bic <- (n2-1)*choose(p,2)

# initialize and populate ics2 (Table 2, columns 4-5)
ics2 <- matrix(NA,6,2)
rownames(ics2) <- c('Transitivity AR model',
                    'Global AR model',
                    'Edge-wise AR model',
                    'Edge-wise mean model',
                    'Degree parameter mean model',
                    'Global mean model')
colnames(ics2) <- c('AIC','BIC')

# 1. Transitivity AR model
# 2p+2 parameters
gamma_tar <- model_probs(fit_man2,UV_man2,X_man2)$gamma
ll_tar <- ar_loglike(X_man2[,,-1],gamma_tar)
ics2[1,1] <- 2*(2*p + 2) - 2*ll_tar
ics2[1,2] <- log(n2_bic)*(2*p + 2) - 2*ll_tar

# 2. Simple AR model
# 2 parameters
fit_sar <- simple_ar_fit(X_man2)
gamma_sar <- fit_sar[1] + X_man2[,,-n2]*(1 - sum(fit_sar))
ll_sar <- ar_loglike(X_man2[,,-1],gamma_sar)
ics2[2,1] <- 2*2 - 2*ll_sar
ics2[2,2] <- log(n2_bic)*2 - 2*ll_sar

# 3. Edgewise AR model
# 2(p \choose 2) parameters
fit_ear <- edge_ar_fit(X_man2)
gamma_ear <- array(rep(fit_ear$A,n2-1),c(p,p,n2-1)) + X_man2[,,-n2]*(1 - array(rep(fit_ear$A+fit_ear$B,n2-1),c(p,p,n2-1)))
ll_ear <- ar_loglike(X_man2[,,-1],gamma_ear)
ics2[3,1] <- 2*(2*choose(p,2)) - 2*ll_ear
ics2[3,2] <- log(n2_bic)*(2*choose(p,2)) - 2*ll_ear

# 4. Edge means model
# (p \choose 2) parameters
X_man_mean <- apply(X_man2[,,-1],c(1,2),mean)
gamma_mean <- array(rep(X_man_mean,n2-1),c(p,p,n2-1))
ll_mean <- ar_loglike(X_man2[,,-1],gamma_mean)
ics2[4,1] <- 2*choose(p,2) - 2*ll_mean
ics2[4,2] <- log(n2_bic)*choose(p,2) - 2*ll_mean

# 5. Degree parameters/1-dim RDPG model
# p parameters
eigXmean <- eigen(X_man_mean)
X_man_deg <- pmin(pmax(eigXmean$values[1]*tcrossprod(eigXmean$vectors[,1]),0),1)
gamma_deg <- array(rep(X_man_deg,n2-1),c(p,p,n2-1))
ll_deg <- ar_loglike(X_man2[,,-1],gamma_deg)
ics2[5,1] <- 2*p - 2*ll_deg
ics2[5,2] <- log(n2_bic)*p - 2*ll_deg

# 6. ER/global mean
p_mean <- mean(X_man2[,,-1])
gamma_global <- array(p_mean,c(p,p,n2-1))
ll_global <- ar_loglike(X_man2[,,-1],gamma_global)
ics2[6,1] <- 2 - 2*ll_global
ics2[6,2] <- log(n2_bic) - 2*ll_global

# reorder for presentation
ics2 <- ics2[c(1,2,3,5,6,4),]

# save ics2 (Table 2, columns 4-5)
saveRDS(ics2,file='real_data/data/ics_man2.rds')
