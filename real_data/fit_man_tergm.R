# refitting with btergm

# load data
setwd('~/packages/arnetworks/real_data/')
# load packages (incl. arnetworks)
library(devtools)
library(pROC)
library(latex2exp)
library(btergm)
library(statnet) # includes package tergm

# source code from arnetworks package
load_all()
source('realdataFunctions.R')
#source('Functions.R')
#source('Estim.R')

# manufacturing network
X_man <- readRDS('data/X_man.rds')
p <- dim(X_man)[1]
n <- dim(X_man)[3]
UV_man <- arnetworks::statsTransitivity(X_man)

# calculate U/V stats
# splitting into two regimes for fitting
# first regime snapshots 1 to 131
X_man1 <- X_man[,,1:13]
UV_man1 <- arnetworks::statsTransitivity(X_man1) # scaled by p-2
n1 <- dim(X_man1)[3]

# second regime snapshots 14 to 39
X_man2 <- X_man[,,14:39]
UV_man2 <- arnetworks::statsTransitivity(X_man2) # scaled by p-2
n2 <- dim(X_man2)[3]

# store as complete network list
nets <- U <- V <- list()
for(ii in 1:n){
  nets[[ii]] <- network::network(X_man[,,ii],directed=FALSE)
  U[[ii]] <- UV_man$U[,,ii]
  V[[ii]] <- UV_man$V[,,ii]
}

#### btergm ####

# tergm for first regime

# data
net1 <- list()
U1 <- list()
V1 <- list()
for (ii in 2:n1) {
  # store networks
  nw <- network::network(X_man1[,,ii],directed=FALSE)  # create network object
  net1[[(ii-1)]] <- nw          # add network to the list
  U1[[(ii-1)]] <- UV_man1$U[,,ii-1]
  V1[[(ii-1)]] <- UV_man1$V[,,ii-1]
}

# DON'T RUN: btergm version
# # fit with just U (common neighbours) effect
# fit1 <- mtergm(net1 ~ edges + edgecov(U1))
# summary(fit1)

# tergm for second regime

# data
net2 <- list()
U2 <- list()
V2 <- list()
for (ii in 2:n2) {
  # store networks
  nw <- network::network(X_man2[,,ii],directed=FALSE)  # create network object
  net2[[(ii-1)]] <- nw          # add network to the list
  U2[[(ii-1)]] <- UV_man2$U[,,ii-1]
  V2[[(ii-1)]] <- UV_man2$V[,,ii-1]
}

# DON'T RUN: btergm version
# # fit with just U (common neighbours) effect
# fit2 <- mtergm(net2 ~ edges + edgecov(U2))
# summary(fit2)

# similarly find strong/significant transitivity effect of previous common neighbours on edge presence, although
# this TERGM does not isolate the effects on formation and dissolution conditional on the previous edge
# status, nor can it incorporate node-specific popularity effects (xis and etas in our model)

#### tergm ####

# fitting a sequence of STERGMs with edge + edgecov

# allocate space for coefficients
coef_ecov <- pval_ecov <- se_ecov <- matrix(NA,n-1,6)
# fit sequence of models
for(ii in 1:(n-1)){
  temp <- tergm(list(nets[[ii]],nets[[(ii+1)]]) ~ Form(~edges + edgecov(U[[ii]]) + edgecov(V[[ii]])) + Persist(~edges + edgecov(U[[ii]])+ edgecov(V[[ii]])),estimate='CMLE')
  coef_ecov[ii,] <- coef(temp)
  se_ecov[ii,] <- summary(temp)$coefficients[,2]
  pval_ecov[ii,] <- summary(temp)$coefficients[,5]
}
# NOTE: seems like edgecov breaks when the model is specified as Form/Diss instead of Form/Persist
# NOTE: networks are too small and sparse to fit node-specific sociality effects to individual transitions, MPLE
# does not exist

# remove 13 and 14 to account for the period change

# plot over time
pdf('fit_plots_man/tergm_sequence.pdf',width=8,height=8)
matplot(c(1:12,15:(n-1)),coef_ecov[c(1:12,15:(n-1)),],ylim=c(-30,60),xlab='Week',ylab='Coef. estimate',
        lty=1,pch=16,type='p',main='TERGM transition model parameter estimates',
        cex.lab=1.4,cex.main=1.4)
#matplot(coef_ecov - 2*se_ecov,type='p',lty=1,pch=12,add=TRUE)
abline(v=13.5,lty=2) # cover up transition between regimes
# add legend and label regimes
legend(x=20,y=-7,ncol=2,col=1:6,lty=0,pch=16,cex=1.4,
       legend=c(TeX('$\\hat{\\beta}_0^+'),
                TeX('$\\hat{\\beta}_U^+'),
                TeX('$\\hat{\\beta}_V^+'),
                TeX('$\\hat{\\beta}_0^-'),
                TeX('$\\hat{\\beta}_U^-'),
                TeX('$\\hat{\\beta}_V^-')))
text(x=c(3,22.5),y=60,labels=c('Period 1','Period 2'),pos=4,cex=1.3)
dev.off()

# fitting two STERGMs with edge + sociality

netsn1 <- NetSeries(net1); netsn2 <- NetSeries(net2)

fit1soc <- tergm(net1 ~ Form(~edges + sociality) + Persist(~edges + sociality),estimate='CMLE')
fit2soc <- tergm(net2 ~ Form(~edges + sociality) + Persist(~edges + sociality),estimate='CMLE')

# NOTE: fitting with sociality and triangles, no longer fast with pseudolikelihood, requires slow
# MCMLE approach ... running started ~ 11:40AM for small (n=12) regime 1, still on iteration 1 after 5 minutes

# load and compare to arnetworks heterogeneity parameters

# period 1
pdf('fit_plots_man/tergm_degree.pdf',width=10,height=8)
par(mfrow=c(2,2),mar=c(4,5,4,1))
fit1_arnet <- readRDS('data/fit_man1_imom.rds')
# plot xi-hat against formation socialitys
plot(fit1_arnet$xi,c(0,coef(fit1soc)[2:106]),
     xlab=TeX('AR network $\\hat{\\xi}_i$'),
     ylab=TeX('STERGM $\\hat{\\alpha}_i^+$'),
     main='Formation parameter estimates, period 1',
     cex.lab=1.3,cex.main=1.3)
# plot eta-hat against persistence sociality
plot(fit1_arnet$eta,c(0,coef(fit1soc)[108:212]),
     xlab=TeX('AR network $\\hat{\\eta}_i$'),
     ylab=TeX('STERGM $\\hat{\\alpha}_i^-$'),
     main='Dissolution parameter estimates, period 1',
     cex.lab=1.3,cex.main=1.3)

# NOTE: bad node is 80, it is extremely outlying in this plot. It is farther from the bulk of
# degree vs sociality than degree vs xi, and it is indeed below the bulk of degree vs triangle participation,
# ie overall it has fewer common neighbors than other nodes of similar degree

# period 2
fit2_arnet <- readRDS('data/fit_man2_imom.rds')
# plot xi-hat against formation sociality
plot(fit2_arnet$xi[-86],c(0,coef(fit2soc)[2:106])[-86],
     xlab=TeX('AR network $\\hat{\\xi}_i$'),
     ylab=TeX('STERGM $\\hat{\\alpha}_i^+$'),
     main='Formation parameter estimates, period 2',
     cex.lab=1.3,cex.main=1.3)
# plot eta-hat against persistence sociality
plot(fit2_arnet$eta[-86],c(0,coef(fit2soc)[108:212])[-86],
     xlab=TeX('AR network $\\hat{\\eta}_i$'),
     ylab=TeX('STERGM $\\hat{\\alpha}_i^-$'),
     main='Dissolution parameter estimates, period 2',
     cex.lab=1.3,cex.main=1.3)
dev.off()

# NOTE: MPLE does not exist for period 2, when refit with problem node 86 removed,
# remaining parameters are essentially the same, leave as is, removing node 86 from
# the plots as its sociality parameter is unreliable

# DO NOT RUN: does not coverge and/or optimization takes too long to run
# fitting two STERGMs with edge + sociality + triangle (+ threepath)

# fit1tri <- tergm(netsn1 ~ Form(~edges + triangle) + Persist(~edges + triangle),estimate='CMLE') # error MCMLE stuck
# fit2tri <- tergm(netsn2 ~ Form(~edges + triangle) + Persist(~edges + triangle),estimate='CMLE')
# fit1trithree <- tergm(netsn1 ~ Form(~edges + triangle + threetrail) + Persist(~edges + triangle + threetrail),estimate='CMLE')
# fit2trithree <- tergm(netsn2 ~ Form(~edges + triangle + threetrail) + Persist(~edges + triangle + threetrail),estimate='CMLE')
# save(fit1tri,fit2tri,fit1trithree,fit2trithree,file='data/fit_tergms.RData')

# fit1tri terminates w/error
# fit2tri terminates w/error
# fit1trithree does not converge (17/60 iterations in ~4 days)
# fit2trithree not run
