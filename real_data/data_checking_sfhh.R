# Summary plots for SFHH conference proximity data

#### helper functions ####

hollowize <- function(M){
  M - diag(diag(M))
}

# colorblind palete
cbp <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# black and white indicator
bw <- FALSE

#### importing data ####
setwd('~/packages/arnetworks/real_data/')

# importing SFHH conference data
raw_data <- read.csv('data/tij_SFHH.dat_.gz',
                     sep=' ',header=FALSE,col.names = c('t','i','j'))
# convert times to walltimes (decimal hours)
raw_data$walltime <- paste0(1 + (raw_data$t %/% (60*60*24)),'-',
                            (raw_data$t %% (60*60*24)) %/% (60*60),':',
                            (raw_data$t %% (60*60)) %/% 60,':',
                            (raw_data$t %% 60))
# network construction parameters
dT <- 1200 # 20 minute windowing
tmin <- min(raw_data$t) + 2*60*60
tmax <- min(raw_data$t) + 9*60*60
dynamic_thresh <- 50

##### network construction #####

# construct network snapshots for a given time window

# construct snapshots by reading edgelist
# get nodes
nodes <- unique(c(raw_data$i,raw_data$j))
p_raw <- length(nodes)
# get snapshot windows
dn <- seq(tmin,tmax+dT,dT)
n <- length(dn)-1
X_raw <- array(0,c(p_raw,p_raw,n))

for(rr in 1:nrow(raw_data)){
  temp <- which(raw_data$t[rr] >= dn)
  if(length(temp) > 0){
    tt <- max(which(raw_data$t[rr] >= dn))
    if(tt < (n+1)){
      ii <- which(nodes==raw_data$i[rr])
      jj <- which(nodes==raw_data$j[rr])
      if(ii != jj){
        X_raw[ii,jj,tt] <- 1
        # symmetrize
        X_raw[jj,ii,tt] <- 1
      }
    }
  }
}
#print(dim(X_raw))

#### Dynamic metrics ####

# metric to decide whether to keep -> total dissolves and grows for a given node
grow <- X_raw[,,-1]*(1 - X_raw[,,-n])
diss <- X_raw[,,-n]*(1 - X_raw[,,-1])
dynamic_count <- apply(grow+diss,1,sum)

# pdf('data_plots_sfhh/sfhh_nodehetraw.pdf',width=7,height=5)
# hist(dynamic_count,40,
#      main='Total dynamic activity by node (raw data)')
# dev.off()

node_keep <- which(dynamic_count >= dynamic_thresh)
X <- X_raw[node_keep,node_keep,]
p <- dim(X)[1]

# save clean data
saveRDS(X,file='data/X_sfhh.rds')

# subset grow/dissolve arrays
grow_sub <- grow[node_keep,node_keep,]
diss_sub <- diss[node_keep,node_keep,]

# now edge density on the subnetwork
dens <- apply(X,3,mean)

pdf('data_plots_sfhh/sfhh_dens.pdf',width=7,height=5)
plot(dens,type='b',main='Edge density',
     xlab='20 minute interval',ylab='Edge density',
     cex.main=1.5,cex.lab=1.4)
dev.off()

# scaled growing/dissolving
grt <- apply(grow_sub,3,sum)/p
dst <- apply(diss_sub,3,sum)/p

# pdf('data_plots_sfhh/sfhh_grds.pdf',width=7,height=5)
# plot(grt,type='b',col=cbp[3],main='Dynamic activity',
#      ylim=c(0,max(grt)+1),ylab='Grown or dissolved edges / nNodes',xlab='')
# lines(dst,type='b',col=cbp[7])
# abline(h=mean(grt),col=cbp[3],lty=3)
# abline(h=mean(dst),col=cbp[7],lty=3)
# abline(h=1,lty=2)
# dev.off()

# O(p^2) scaled growth and dissolution
grt_p2 <- apply(grow_sub,3,sum)/(p*(p-1))
dst_p2 <- apply(diss_sub,3,sum)/(p*(p-1))

if(bw){
  pdf('data_plots_sfhh/sfhh_grds_p2_bw.pdf',width=7,height=5)
  plot(grt_p2,type='l',lty=1,main='Dynamic activity',
       ylim=c(0,0.05),ylab='Dynamic activity - form/dissolve',
       xlab='20 minute interval',cex.main=1.5,cex.lab=1.4)
  lines(grt_p2,type='p',pch=17)
  lines(dst_p2,type='l',lty=2)
  lines(dst_p2,type='p',pch=25,bg='white')
  #abline(h=mean(grt_p2),col=cbp[3],lty=3)
  #abline(h=mean(dst_p2),col=cbp[7],lty=3)
  #abline(h=1,lty=2)
  dev.off()
} else{
  pdf('data_plots_sfhh/sfhh_grds_p2.pdf',width=7,height=5)
  plot(grt_p2,type='b',col=cbp[3],main='Dynamic activity',
       ylim=c(0,0.05),ylab='Dynamic activity - form/dissolve',
       xlab='20 minute interval',cex.main=1.5,cex.lab=1.4)
  lines(dst_p2,type='b',col=cbp[7])
  abline(h=mean(grt_p2),col=cbp[3],lty=3)
  abline(h=mean(dst_p2),col=cbp[7],lty=3)
  #abline(h=1,lty=2)
  dev.off()
}

# normalized growth and dissolution
grow_norm <- grow_sub
grow_norm[(X[,,-n]==1)] <- NA
grt_norm <- apply(grow_norm,3,mean,na.rm=TRUE)

diss_norm <- diss_sub
diss_norm[(X[,,-n]==0)] <- NA
dst_norm <- apply(diss_norm,3,mean,na.rm=TRUE)

# pdf('data_plots_sfhh/sfhh_grds_norm.pdf',width=7,height=5)
# plot(grt_norm,type='b',col=cbp[3],main='Dynamic activity',
#      ylim=c(0,1),ylab='Relative freq.')
# lines(dst_norm,type='b',col=cbp[7])
# abline(h=mean(grt_norm),col=cbp[3],lty=3)
# abline(h=mean(dst_norm),col=cbp[7],lty=3)
# abline(h=1,lty=2)
# dev.off()

# # growth and dissolution prob. after conditioning on density
# e_sparse <- (apply(X,c(1,2),mean) < 0.2)
# e_dense <- (apply(X,c(1,2),mean) >= 0.2)
#
# diss_norm_sparse <- diss_norm_dense <- diss_norm
# grow_norm_sparse <- grow_norm_dense <- grow_norm
# for(kk in 1:dim(diss_norm)[3]){
#   diss_norm_sparse[,,kk][e_dense] <- NA
#   grow_norm_sparse[,,kk][e_dense] <- NA
#   diss_norm_dense[,,kk][e_sparse] <- NA
#   grow_norm_dense[,,kk][e_sparse] <- NA
# }
#
# grt_norm_sparse <- apply(grow_norm_sparse,3,mean,na.rm=TRUE)
# dst_norm_sparse <- apply(diss_norm_sparse,3,mean,na.rm=TRUE)
# grt_norm_dense <- apply(grow_norm_dense,3,mean,na.rm=TRUE)
# dst_norm_dense <- apply(diss_norm_dense,3,mean,na.rm=TRUE)
#
# plot(grt_norm_sparse,type='b',col=cbp[3],main='Normalized dynamic activity, sparse nodes',
#      ylim=c(0,1),ylab='Growth or dissolution prob.')
# lines(dst_norm_sparse,type='b',col=cbp[7])
# abline(h=mean(grt_norm_sparse),col=cbp[3],lty=3)
# abline(h=mean(dst_norm_sparse),col=cbp[7],lty=3)
# abline(h=1,lty=2)
#
# plot(grt_norm_dense,type='b',col=cbp[3],main='Normalized dynamic activity, dense nodes',
#      ylim=c(0,1),ylab='Growth or dissolution prob.')
# lines(dst_norm_dense,type='b',col=cbp[7])
# abline(h=mean(grt_norm_dense),col=cbp[3],lty=3)
# abline(h=mean(dst_norm_dense),col=cbp[7],lty=3)
# abline(h=1,lty=2)

# looking for whether p(grow) = 1 - prob(dissolve) which would imply that there
# is no observed AR structure, this is not the case, both groups show evidence of
# persistence

#### Transitivity metrics ####

# look at grown edge probability, plot against number of common neighbors
# look at dissolved edge probability, plot against previous number of
# common non-neighbors

gr_cn <- NULL
for(tt in 2:n){
  prev_cn <- X[,,tt-1]%*%X[,,tt-1]
  # growing edges
  grow_ind <- as.logical((1-X[,,tt-1])*X[,,tt])
  if(sum(grow_ind) > 0){
    gr_cn <- rbind(gr_cn,cbind(1,prev_cn[grow_ind]))
  }
  # non-growing edges
  nongrow_ind <- as.logical(hollowize((1-X[,,tt-1])*(1-X[,,tt])))
  if(sum(nongrow_ind)>0){
    gr_cn <- rbind(gr_cn,cbind(0,prev_cn[nongrow_ind]))
  }
}

gr_cn_summary <- table(gr_cn[,1],gr_cn[,2])
# rows by indicator of the growth of an edge
# columns by number of previous common neighbors

# summarize this table
if(bw){
  pdf('data_plots_sfhh/sfhh_grtrans_bw.pdf',width=7,height=5)
  plot(as.numeric(colnames(gr_cn_summary)),
       gr_cn_summary[2,]/colSums(gr_cn_summary),
       type='p',pch=17,cex=(3+log(colSums(gr_cn_summary)))/6,
       main='Transitivity effects on grown edges',cex.main=1.5,
       xlab='Number of common neighbours',ylab='Relative freq. formed edge',cex.lab=1.4)
  legend(x=0,y=0.18,pch=17,title='Sample size',
         pt.cex=(3+log(c(1e4,1e2,1)))/6,legend=c(TeX('$10^4$'),TeX('$10^2$'),TeX('$1$')))
  #abline(h=0,lty=2)
  #abline(h=1,lty=2)
  dev.off()
} else{
  pdf('data_plots_sfhh/sfhh_grtrans.pdf',width=7,height=5)
  plot(as.numeric(colnames(gr_cn_summary)),
       gr_cn_summary[2,]/colSums(gr_cn_summary),
       type='p',pch=16,cex=(3+log(colSums(gr_cn_summary)))/6,col=cbp[3],
       main='Transitivity effects on grown edges',cex.main=1.5,
       xlab='Number of common neighbours',ylab='Relative freq. formed edge',cex.lab=1.4)
  legend(x=0,y=0.18,pch=16,col=cbp[3],title='Sample size',
         pt.cex=(3+log(c(1e4,1e2,1)))/6,legend=c(TeX('$10^4$'),TeX('$10^2$'),TeX('$1$')))
  #abline(h=0,lty=2)
  #abline(h=1,lty=2)
  dev.off()
}

# same for dissolving edges
ds_ncn <- NULL
for(tt in 2:n){
  prev_ncn <- X[,,tt-1] %*% hollowize(1-X[,,tt-1])
  # symmetrize
  prev_ncn <- .5*(prev_ncn + t(prev_ncn))
  # dissolving edges
  dissolve_ind <- as.logical(X[,,tt-1]*(1-X[,,tt]))
  if(sum(dissolve_ind) > 0){
    ds_ncn <- rbind(ds_ncn,cbind(1,prev_ncn[dissolve_ind]))
  }
  # non-growing edges
  nondissolve_ind <- as.logical(X[,,tt-1]*X[,,tt])
  if(sum(nondissolve_ind)>0){
    ds_ncn <- rbind(ds_ncn,cbind(0,prev_ncn[nondissolve_ind]))
  }
}

ds_ncn_summary <- table(ds_ncn[,1],ds_ncn[,2])
# rows by indicator of the growth of an edge
# columns by number of previous common neighbors

# only look at columns with enough sample size (>50)
#ds_ncn_summary <- ds_ncn_summary[,colSums(ds_ncn_summary) > 50]

# summarize this table
if(bw){
  pdf('data_plots_sfhh/sfhh_dstrans_bw.pdf',width=7,height=5)
  plot(as.numeric(colnames(ds_ncn_summary)),
       ds_ncn_summary[2,]/colSums(ds_ncn_summary),
       type='p',pch=6,bg='white',cex=(3+log(colSums(ds_ncn_summary)))/6,
       main='Transitivity effects on dissolved edges',cex.main=1.5,
       xlab='Number of disjoint neighbours',ylab='Relative freq. dissolved edge',cex.lab=1.4)
  legend(x=0,y=0.4,pch=6,bg='white',title='Sample size',
         pt.cex=(3+log(c(1e4,1e2,1)))/6,legend=c(TeX('$10^4$'),TeX('$10^2$'),TeX('$1$')))
  # abline(h=0,lty=2)
  # abline(h=1,lty=2)
  dev.off()
} else{
  pdf('data_plots_sfhh/sfhh_dstrans.pdf',width=7,height=5)
  plot(as.numeric(colnames(ds_ncn_summary)),
       ds_ncn_summary[2,]/colSums(ds_ncn_summary),
       type='p',pch=16,cex=(3+log(colSums(ds_ncn_summary)))/6,col=cbp[7],
       main='Transitivity effects on dissolved edges',cex.main=1.5,
       xlab='Number of disjoint neighbours',ylab='Relative freq. dissolved edge',cex.lab=1.4)
  legend(x=0,y=0.4,pch=16,col=cbp[7],title='Sample size',
         pt.cex=(3+log(c(1e4,1e2,1)))/6,legend=c(TeX('$10^4$'),TeX('$10^2$'),TeX('$1$')))
  # abline(h=0,lty=2)
  # abline(h=1,lty=2)
  dev.off()
}


#### Degree heterogeneity ####

# construct the matrices like in the paper

# rough degree parameters for growing edges based on summary counts
gr_nhet <- apply(grow_sub,1,sum)

# pdf('data_plots_sfhh/sfhh_grnhet.pdf',width=7,height=5)
# hist(gr_nhet,20,
#      col=cbp[3],
#      main="Growing edge node heterogeneity")
# dev.off()

# rough degree parameters for dissolving edges based on summary counts
ds_nhet <- apply(diss_sub,1,sum)

# pdf('data_plots_sfhh/sfhh_dsnhet.pdf',width=7,height=5)
# hist(ds_nhet,20,
#      col=cbp[7],
#      main="Dissolving edge node heterogeneity")
# dev.off()

#### comparison of dynamic and density ####

# recall gr_nhet and ds_nhet for the total number of grows/dissolves
# pdf('data_plots_sfhh/sfhh_nodehet.pdf',width=7,height=5)
# hist(gr_nhet + ds_nhet,20,
#      main='Total dynamic activity by node',
#      xlab='dynamic_count')
# dev.off()

# total edges
deg_nhet <- apply(X,1,sum)
# hist(deg_nhet,20,
#      main='Total degree by node')

# scatter
# pdf('data_plots_sfhh/sfhh_degvsdens.pdf',width=7,height=5)
# plot(deg_nhet,gr_nhet + ds_nhet,
#      type='p',xlab='Total edges',ylab='Total dynamic activity',
#      main='Total degree vs dynamic activity by node')
# abline(a=0,b=2,lty=2) # upper bound on dynamic activity
# lines(deg_nhet[order(deg_nhet)],
#       loess(gr_nhet + ds_nhet ~ deg_nhet,degree=1)$fitted[order(deg_nhet)],
#       col=cbp[3],lty=1,lwd=1.5)
# dev.off()

#### looking at local node densities ####

# node densities over time
node_dens <- apply(X,c(1,3),sum)
# local linear smoothing
node_dens_smooth <- apply(node_dens,1,function(x){loess(x ~ c(1:n),degree=1)$fitted})
# criterion to decide which might be nonstationary
nonstat <- apply(node_dens_smooth,2,sd)

# plot
library(scales)

ncol <- 11
colvec <- rep(1:ncol,ceiling(p/ncol))[1:p]

# pdf('data_plots_sfhh/sfhh_nodestat.pdf',width=7,height=5)
# matplot(1:n,node_dens_smooth,type='l',pch=1,lty=1,lwd=.5+nonstat/(.5+max(nonstat)),
#         col=alpha(colvec,alpha=nonstat/max(nonstat)),
#         main="Smoothed degree sequence by node",
#         ylab="Smoothed degree",
#         xlab='Time')
# dev.off()

# show the data for the five 'most nonstationary'
# pdf('data_plots_sfhh/sfhh_nodestatsub.pdf',width=7,height=5)
# matplot(1:n,t(node_dens[order(nonstat,decreasing=TRUE)[1:5],]),
#         type='b',pch=1:5,#col=colvec[order(nonstat,decreasing=TRUE)[1:5]],
#         col=1,lty=1,
#         main='Degree sequences of potential nonstationary nodes (top 5)',
#         ylab='Degree',xlab='Time')
# abline(h=0,lty=2)
# dev.off()

#### persistence summary metrics ####

# similar to transitivity metrics

# persistence model effect interpretation:

# 1. less likely to flip on depending on 3 categories:
# edge at t-2, no edge at t-2 but at t-3, no edge at t-2 AND t-3
# 2. less likely to flip off depending on 3 categories:
# no edge at t-2, edge at t-2 but not at t-3, edges at t-2 AND t-3

# needs t-3 \geq 1 -> starting from t=3, checking whether there was a flip to t=4

# show as a bar plot for the three categories, color by grow/dissolve and
# bar width ~ log sample size

Ap <- X[,,-c(1,2,n)]
Ap1 <- X[,,-c(1,2,3)]
Apm1 <- X[,,-c(1,n-1,n)]
Apm2 <- X[,,-c(n-2,n-1,n)]

# grown edges summary counts
gr_11 <- sum(Ap1*(1-Ap)*Apm1)/2
gr_1t <- sum((1-Ap)*Apm1)/2
gr_21 <- sum(Ap1*(1-Ap)*(1-Apm1)*Apm2)/2
gr_2t <- sum((1-Ap)*(1-Apm1)*Apm2)/2
gr_31 <- sum(Ap1*(1-Ap)*(1-Apm1)*(1-Apm2))/2
gr_3t <- (sum((1-Ap)*(1-Apm1)*(1-Apm2)) - p*n)/2

# dissolved edges summary counts
ds_11 <- sum((1-Ap1)*Ap*(1-Apm1))/2
ds_1t <- sum(Ap*(1-Apm1))/2
ds_21 <- sum((1-Ap1)*Ap*Apm1*(1-Apm2))/2
ds_2t <- sum(Ap*Apm1*(1-Apm2))/2
ds_31 <- sum((1-Ap1)*Ap*Apm1*Apm2)/2
ds_3t <- sum(Ap*Apm1*Apm2)/2

# plotting
library(latex2exp)

# pdf('data_plots_sfhh/sfhh_grpers.pdf',width=7,height=5)
# barplot(c(gr_11/gr_1t,gr_21/gr_2t,gr_31/gr_3t),
#         width=log(c(gr_1t,gr_2t,gr_3t)),
#         col=cbp[3],
#         main='Persistence effects on grown edges',
#         ylab='P(grow edge) @ t',
#         names.arg=c(TeX('$X_{{ij}}^{{t-2}}=1$'),
#                     TeX('$X_{{ij}}^{{t-2}}=0, X_{{ij}}^{{t-3}}=1$'),
#                     TeX('$X_{{ij}}^{{t-2}}=0, X_{{ij}}^{{t-3}}=0$')))
# dev.off()
#
# pdf('data_plots_sfhh/sfhh_dspers.pdf',width=7,height=5)
# barplot(c(ds_11/ds_1t,ds_21/ds_2t,ds_31/ds_3t),
#         width=log(c(ds_1t,ds_2t,ds_3t)),
#         col=cbp[7],
#         main='Persistence effects on dissolved edges',
#         ylab='P(dissolve edge) @ t',
#         names.arg=c(TeX('$X_{{ij}}^{{t-2}}=0$'),
#                     TeX('$X_{{ij}}^{{t-2}}=1, X_{{ij}}^{{t-3}}=0$'),
#                     TeX('$X_{{ij}}^{{t-2}}=1, X_{{ij}}^{{t-3}}=1$')))
# dev.off()

#### density-dependent model summary metrics ####

# density model effect interpretation (in undirected case):

# more likely to flip on if the total number of edges at time t-1 is larger
# more likely to flip off if the total number of edges at time t-1 is smaller

# checking stationarity?
# expit <- function(x){exp(x)/(1 + exp(x))}
#
# rho0 <- runif(1)
# a0 <- 1
# b0 <- 1
# rho <- rep(NA,1000)
# rho[1] <- rho0
# for(ii in 1:999){
#   rho[ii+1] <- expit(a0*rho[ii])*(1-rho[ii]) + (1-expit(b0*(1-rho[ii])))*rho[ii]
# }
# plot(rho,type='l',ylim=c(0,1))

# natural balancing effect?
# so check these two types of density effects (global and local) separately

# GLOBAL
# this first type should just compare previous density against rates of flipping
# as a scatter plot
# ie normalized rates of churn on the y, density or 1-density on the x

gr_rate <- apply(grow_sub,3,sum) / apply(1-X[,,-n],3,function(x){sum(x) - p})
ds_rate <- apply(diss_sub,3,sum) / apply(X[,,-n],3,sum)

# comparing grown edge rate to previous density
# pdf('data_plots_sfhh/sfhh_grgdens.pdf',width=7,height=5)
# plot(dens[-n],gr_rate,
#      main='Global density effects on grown edges',
#      col=cbp[3],pch=16,
#      xlab='Edge density @ t-1',
#      ylab='P(grow edge) @ t')
# dev.off()
#
# # comparing dissolved edge rate to previous density
# pdf('data_plots_sfhh/sfhh_dsgdens.pdf',width=7,height=5)
# plot(dens[-n],ds_rate,
#      main='Global density effects on dissolved edges',
#      col=cbp[7],pch=16,
#      xlab='Edge density @ t-1',
#      ylab='P(dissolve edge) @ t')
# dev.off()

# LOCAL
# 1. more likely to flip on if the total degree of i and j at time t-1 is larger
# 2. more likely to flip off if the total degree of i and j at time t-1 is smaller
# note these should ignore the edge i,j
pairsums <- function(v){sapply(v,function(x){sapply(v,function(y){x+y})})}

gr_localdens <- NULL
ds_localdens <- NULL
for(tt in 2:n){
  prev_localdens <- pairsums(rowSums(X[,,tt-1]))
  # growing edges
  grow_ind <- as.logical((1-X[,,tt-1])*X[,,tt])
  if(sum(grow_ind) > 0){
    gr_localdens <- rbind(gr_localdens,cbind(1,prev_localdens[grow_ind]))
  }
  # non-growing edges
  nongrow_ind <- as.logical(hollowize((1-X[,,tt-1])*(1-X[,,tt])))
  if(sum(nongrow_ind)>0){
    gr_localdens <- rbind(gr_localdens,cbind(0,prev_localdens[nongrow_ind]))
  }
  # dissolving edges
  dissolve_ind <- as.logical(X[,,tt-1]*(1-X[,,tt]))
  if(sum(dissolve_ind) > 0){
    ds_localdens <- rbind(ds_localdens,cbind(1,prev_localdens[dissolve_ind]))
  }
  # non-dissolving edges
  nondissolve_ind <- as.logical(X[,,tt-1]*X[,,tt])
  if(sum(nondissolve_ind)>0){
    ds_localdens <- rbind(ds_localdens,cbind(0,prev_localdens[nondissolve_ind]))
  }
}

gr_localdens_summary <- table(gr_localdens[,1],gr_localdens[,2])
ds_localdens_summary <- table(ds_localdens[,1],ds_localdens[,2])
# rows by indicator of the growth of an edge
# columns by number of previous common neighbors

# summarize this table
# pdf('data_plots_sfhh/sfhh_grldens.pdf',width=7,height=5)
# plot(as.numeric(colnames(gr_localdens_summary)),
#      gr_localdens_summary[2,]/colSums(gr_localdens_summary),
#      type='p',pch=16,cex=log(colSums(gr_localdens_summary))/8,col=cbp[3],
#      ylim=c(0,.2),
#      main='Local density effects on grown edges',
#      xlab='Degree i + Degree j @ t-1',ylab='P(grow edge) @ t')
# abline(h=0,lty=2)
# abline(h=1,lty=2)
# dev.off()
#
# pdf('data_plots_sfhh/sfhh_dsldens.pdf',width=7,height=5)
# plot(as.numeric(colnames(ds_localdens_summary)),
#      ds_localdens_summary[2,]/colSums(ds_localdens_summary),
#      type='p',pch=16,cex=log(colSums(ds_localdens_summary))/6,col=cbp[7],
#      main='Local density effects on dissolved edges',
#      xlab='Degree i + Degree j @ t-1',ylab='P(dissolve edge) @ t')
# abline(h=0,lty=2)
# abline(h=1,lty=2)
# dev.off()

# correlation between transitivity and local density effects?
gr_transdenscor <- cor(gr_cn[,2],gr_localdens[,2])
# conference: 0.473

ds_transdenscor <- cor(ds_ncn[,2],ds_localdens[,2])
# conference: 0.920
