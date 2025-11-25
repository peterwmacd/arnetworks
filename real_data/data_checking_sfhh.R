# Preliminary analysis and plots for SFHH conference proximity data

# libraries
library(latex2exp)

# helpers
source('real_data/realdataFunctions.R')

# colorblind palete
cbp <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#### (1) importing data ####

# importing SFHH conference data (Genois & Barrat 2018), included in repository
# additional background/data terms and conditions at
# http://www.sociopatterns.org/datasets/sfhh-conference-data-set/
raw_data <- read.csv('real_data/data/tij_SFHH.dat_.gz',
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
dynamic_thresh <- 50 # keep nodes with at least 50 connections

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

#### (3) dynamic metrics ####

# metric to decide whether to keep, total dissolves and grows for a given node
grow <- X_raw[,,-1]*(1 - X_raw[,,-n])
diss <- X_raw[,,-n]*(1 - X_raw[,,-1])
dynamic_count <- apply(grow+diss,1,sum)

node_keep <- which(dynamic_count >= dynamic_thresh)
X <- X_raw[node_keep,node_keep,]
p <- dim(X)[1]

# save clean network snapshots
saveRDS(X,file='real_data/data/X_sfhh.rds')

# subset grow/dissolve arrays
grow_sub <- grow[node_keep,node_keep,]
diss_sub <- diss[node_keep,node_keep,]

# now edge density on the subnetwork
dens <- apply(X,3,mean)

# produces sfhh_dens.pdf (Figure S8, left panel)
pdf('real_data/data_plots_sfhh/sfhh_dens.pdf',width=7,height=5)
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

# produce sfhh_grds_p2.pdf (Figure S8, right panel)
pdf('real_data/data_plots_sfhh/sfhh_grds_p2.pdf',width=7,height=5)
plot(grt_p2,type='b',col=cbp[3],main='Dynamic activity',
     ylim=c(0,0.05),ylab='Dynamic activity - form/dissolve',
     xlab='20 minute interval',cex.main=1.5,cex.lab=1.4)
lines(dst_p2,type='b',col=cbp[7])
abline(h=mean(grt_p2),col=cbp[3],lty=3)
abline(h=mean(dst_p2),col=cbp[7],lty=3)
#abline(h=1,lty=2)
dev.off()

#### (4) transitivity metrics ####

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

# produce sfhh_grtrans.pdf (Figure S9, left panel)
pdf('real_data/data_plots_sfhh/sfhh_grtrans.pdf',width=7,height=5)
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

# produce sfhh_dstrans.pdf (Figure S9, right panel)
pdf('real_data/data_plots_sfhh/sfhh_dstrans.pdf',width=7,height=5)
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
