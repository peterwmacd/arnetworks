# Preliminary analysis and plots for Polish manufacturing email data

# libraries
library(latex2exp)
library(networkDynamicData)

# helpers
source('real_data/realdataFunctions.R')

# colorblind palete
cbp <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#### (1) importing data ####

# importing for Polish email data
data("manufacturingEmails")
package_data <- as.data.frame(manufacturingEmails)
raw_data <- data.frame(t=ceiling((package_data$onset - 1262482809)/(60*60*24)),
                       i=package_data$tail,
                       j=package_data$head)
# network construction parameters
dT <- 7 # weekly windowing
tmin <- min(raw_data$t)+1
tmax <- max(raw_data$t)
dynamic_thresh <- 100

# 167 nodes, 82613 email events
# note converted times in seconds to times in days (ie kth day of the study),
# ranges from 1 to 272

##### (2) network construction #####

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

##### (3) network construction by day #####

# construct network snapshots for a given time window

# construct snapshots by reading edgelist
# get nodes
nodes <- unique(c(raw_data$i,raw_data$j))
p_raw <- length(nodes)
# get snapshot windows
dn_day <- seq(tmin,tmax+1,1)
n_day <- length(dn_day)-1
X_raw_day <- array(0,c(p_raw,p_raw,n_day))

for(rr in 1:nrow(raw_data)){
  temp <- which(raw_data$t[rr] >= dn_day)
  if(length(temp) > 0){
    tt <- max(which(raw_data$t[rr] >= dn_day))
    if(tt < (n_day+1)){
      ii <- which(nodes==raw_data$i[rr])
      jj <- which(nodes==raw_data$j[rr])
      if(ii != jj){
        X_raw_day[ii,jj,tt] <- 1
        # symmetrize
        X_raw_day[jj,ii,tt] <- 1
      }
    }
  }
}

#### (4) metadata ####

# calculate level in the hierarchy
org <- as.data.frame(manufacturingReportsTo)
employee <- data.frame(list(id=1:167,level=1))

for(jj in 1:10){
  for(ii in 1:nrow(org)){
    child <- org[ii,1]
    parent <- org[ii,2]
    if(child != parent){
      employee[parent,2] <- max(employee[parent,2],1+employee[child,2])
    }
  }
}
# update singleton nodes
# technical emails (code as 0)
employee[c(4,10,21,23,24,26,46),2] <- 0
# former employees (code as -1)
employee[c(51,75,87,93,111,139),2] <- -1

node_level <- rep(NA,length(nodes))
for(ii in 1:length(nodes)){
  node_level[ii] <- employee[which(employee$id==nodes[ii]),2]
}

#### (5) dynamic metrics ####

# metric to decide whether to keep -> total dissolves and grows for a given node
grow <- X_raw[,,-1]*(1 - X_raw[,,-n])
diss <- X_raw[,,-n]*(1 - X_raw[,,-1])
dynamic_count <- apply(grow+diss,1,sum)

node_keep <- which(dynamic_count >= dynamic_thresh)
X <- X_raw[node_keep,node_keep,]
X_day <- X_raw_day[node_keep,node_keep,]
p <- dim(X)[1]
# hierarchical level of active nodes
levels <- node_level[node_keep]
# note didn't keep any technical accounts or former employees

# save clean data and metadata
saveRDS(X,file='real_data/data/X_man.rds')
saveRDS(X_day,file='real_data/data/X_man_day.rds')
saveRDS(levels,file='real_data/data/levels_man.rds')

# subset grow/dissolve arrays
grow_sub <- grow[node_keep,node_keep,]
diss_sub <- diss[node_keep,node_keep,]

# now edge density on the subnetwork
dens <- apply(X,3,mean)

# produce man_dens.pdf (Figure S2, left panel)
pdf('real_data/data_plots_man/man_dens.pdf',width=7,height=5)
plot(dens,type='b',main='Edge density',
     ylab='Edge density',xlab='Week',cex.main=1.5,cex.lab=1.4)
dev.off()

# edge density on the subnetwork by day
dens <- apply(X_day,3,mean)

# produce man_dens_day.pdf (Figure S4, top panel)
pdf('real_data/data_plots_man/man_day_dens.pdf',width=10,height=5)
plot(dens,type='l',main='Edge density (daily)',
     ylab='Edge density',xlab='Day',cex.main=1.5,cex.lab=1.4,cex=0.8)
lines(dens,cex=0.75,type='p')
dev.off()

# scaled growing/dissolving
grt <- apply(grow_sub,3,sum)/p
dst <- apply(diss_sub,3,sum)/p

# O(p^2) scaled growth and dissolution
grt_p2 <- apply(grow_sub,3,sum)/(p*(p-1))
dst_p2 <- apply(diss_sub,3,sum)/(p*(p-1))

# produce man_grds_p2.pdf (Figure S2, right panel)
pdf('real_data/data_plots_man/man_grds_p2.pdf',width=7,height=5)
plot(grt_p2,type='b',col=cbp[3],main='Dynamic activity',
     ylim=c(0,0.1),ylab='Dynamic activity - form/dissolve',xlab='Week',cex.main=1.5,cex.lab=1.4)
lines(dst_p2,type='b',col=cbp[7])
abline(h=mean(grt_p2),col=cbp[3],lty=3)
abline(h=mean(dst_p2),col=cbp[7],lty=3)
#abline(h=1,lty=2)
dev.off()

#### (6) dynamic metrics by day ####

grow_day <- X_day[,,-1]*(1 - X_day[,,-n_day])
diss_day <- X_day[,,-n_day]*(1 - X_day[,,-1])

grt_p2_day <- apply(grow_day,3,sum)/(p*(p-1))
dst_p2_day <- apply(diss_day,3,sum)/(p*(p-1))

# produce man_gr_day.pdf (Figure S4, middle panel)
pdf('real_data/data_plots_man/man_gr_day.pdf',width=10,height=5)
plot(grt_p2_day,type='l',col=cbp[3],main='Edge formation (daily)',
     ylim=c(0,0.1),ylab='Formation rate',xlab='Day',cex.main=1.5,cex.lab=1.4,cex=0.8)
lines(grt_p2_day,col=cbp[3],type='p',cex=0.75)
#lines(dst_p2_day,type='p',col=cbp[7],cex=0.8)
#abline(h=mean(grt_p2_day),col=cbp[3],lty=3)
#abline(h=mean(dst_p2_day),col=cbp[7],lty=3)
#abline(h=1,lty=2)
dev.off()

# produce man_ds_day.pdf (Figure S4, bottom panel)
pdf('real_data/data_plots_man/man_ds_day.pdf',width=10,height=5)
plot(dst_p2_day,type='l',col=cbp[7],main='Edge dissolution (daily)',
     ylim=c(0,0.1),ylab='Dissolution rate',xlab='Day',cex.main=1.5,cex.lab=1.4,cex=0.8)
lines(dst_p2_day,col=cbp[7],type='p',cex=0.75)
#lines(dst_p2_day,type='p',col=cbp[7],cex=0.8)
#abline(h=mean(grt_p2_day),col=cbp[3],lty=3)
#abline(h=mean(dst_p2_day),col=cbp[7],lty=3)
#abline(h=1,lty=2)
dev.off()

#### (7) transitivity metrics ####

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

# produce man_grtrans.pdf (Figure S3, left panel)
pdf('real_data/data_plots_man/man_grtrans.pdf',width=7,height=5)
plot(as.numeric(colnames(gr_cn_summary)),
     gr_cn_summary[2,]/colSums(gr_cn_summary),
     type='p',pch=16,cex=(3+log(colSums(gr_cn_summary)))/6,col=cbp[3],
     main='Transitivity effects on grown edges',cex.main=1.5,
     xlab='Number of common neighbours',ylab='Relative freq. formed edge',cex.lab=1.4)
legend(x=17.5,y=0.25,pch=16,col=cbp[3],title='Sample size',
       pt.cex=(3+log(c(1e4,1e2,1)))/6,legend=c(TeX('$10^4$'),TeX('$10^2$'),TeX('$1$')))
#abline(h=0,lty=2)
#abline(h=1,lty=2)
dev.off()

# same transitivity analysis for dissolving edges
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

# produce man_dstrans.pdf (Figure S3, right panel)
pdf('real_data/data_plots_man/man_dstrans.pdf',width=7,height=5)
plot(as.numeric(colnames(ds_ncn_summary)),
     ds_ncn_summary[2,]/colSums(ds_ncn_summary),
     type='p',pch=16,cex=(3+log(colSums(ds_ncn_summary)))/6,col=cbp[7],
     main='Transitivity effects on dissolved edges',cex.main=1.5,
     xlab='Number of disjoint neighbours',ylab='Relative freq. dissolved edge',cex.lab=1.4)
legend(x=0,y=1,pch=16,col=cbp[7],title='Sample size',
       pt.cex=(3+log(c(1e4,1e2,1)))/6,legend=c(TeX('$10^4$'),TeX('$10^2$'),TeX('$1$')))
# abline(h=0,lty=2)
# abline(h=1,lty=2)
dev.off()

#### (8) transitivity metrics by day ####

# look at grown edge probability, plot against number of common neighbors
# look at dissolved edge probability, plot against previous number of
# common non-neighbors

gr_cn_day <- NULL
for(tt in 2:n_day){
  prev_cn_day <- X_day[,,tt-1]%*%X_day[,,tt-1]
  # growing edges
  grow_ind_day <- as.logical((1-X_day[,,tt-1])*X_day[,,tt])
  if(sum(grow_ind_day) > 0){
    gr_cn_day <- rbind(gr_cn_day,cbind(1,prev_cn_day[grow_ind_day]))
  }
  # non-growing edges
  nongrow_ind_day <- as.logical(hollowize((1-X_day[,,tt-1])*(1-X_day[,,tt])))
  if(sum(nongrow_ind_day)>0){
    gr_cn_day <- rbind(gr_cn_day,cbind(0,prev_cn_day[nongrow_ind_day]))
  }
}

gr_cn_summary_day <- table(gr_cn_day[,1],gr_cn_day[,2])
# rows by indicator of the growth of an edge
# columns by number of previous common neighbors

# produce man_grtrans_day.pdf (Figure S5, right panel)
pdf('real_data/data_plots_man/man_grtrans_day.pdf',width=7,height=5)
plot(as.numeric(colnames(gr_cn_summary_day)),
     gr_cn_summary_day[2,]/colSums(gr_cn_summary_day),
     type='p',pch=16,cex=(3+log(colSums(gr_cn_summary_day)))/6,col=cbp[3],
     main='Transitivity effects on grown edges (daily data)',cex.main=1.5,
     xlab='Number of common neighbours',ylab='Relative freq. formed edge',cex.lab=1.4,ylim=c(0,0.16))
legend(x=15,y=0.155,pch=16,col=cbp[3],title='Sample size',
       pt.cex=(3+log(c(1e4,1e2,1)))/6,legend=c(TeX('$10^4$'),TeX('$10^2$'),TeX('$1$')))
#abline(h=0,lty=2)
#abline(h=1,lty=2)
dev.off()

# same daily transitivity analysis for dissolving edges
ds_ncn_day <- NULL
for(tt in 2:n_day){
  prev_ncn_day <- X_day[,,tt-1] %*% hollowize(1-X_day[,,tt-1])
  # symmetrize
  prev_ncn_day <- .5*(prev_ncn_day + t(prev_ncn_day))
  # dissolving edges
  dissolve_ind_day <- as.logical(X_day[,,tt-1]*(1-X_day[,,tt]))
  if(sum(dissolve_ind_day) > 0){
    ds_ncn_day <- rbind(ds_ncn_day,cbind(1,prev_ncn_day[dissolve_ind_day]))
  }
  # non-growing edges
  nondissolve_ind_day <- as.logical(X_day[,,tt-1]*X_day[,,tt])
  if(sum(nondissolve_ind_day)>0){
    ds_ncn_day <- rbind(ds_ncn_day,cbind(0,prev_ncn_day[nondissolve_ind_day]))
  }
}

ds_ncn_summary_day <- table(ds_ncn_day[,1],ds_ncn_day[,2])
# rows by indicator of the growth of an edge
# columns by number of previous common neighbors

# produce man_dstrans_day.pdf (Figure S5, left panel)
pdf('real_data/data_plots_man/man_dstrans_day.pdf',width=7,height=5)
plot(as.numeric(colnames(ds_ncn_summary_day)),
     ds_ncn_summary_day[2,]/colSums(ds_ncn_summary_day),
     type='p',pch=16,cex=(3+log(colSums(ds_ncn_summary_day)))/6,col=cbp[7],
     main='Transitivity effects on dissolved edges (daily data)',cex.main=1.5,
     xlab='Number of disjoint neighbours',ylab='Relative freq. dissolved edge',cex.lab=1.4,ylim=c(0,1))
legend(x=0,y=0.4,pch=16,col=cbp[7],title='Sample size',
       pt.cex=(3+log(c(1e4,1e2,1)))/6,legend=c(TeX('$10^4$'),TeX('$10^2$'),TeX('$1$')))
# abline(h=0,lty=2)
# abline(h=1,lty=2)
dev.off()
