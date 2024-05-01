simuM = function(p, n, theta, eta,  a1=2, b1=2, a2=2, b2=2, initial = NA ) {
  # p: No. of nodes
  # n: No. of observed networks
  # a: coefficient of No. of common friends
  # b: coefficient of No. of non-(common friends) 
  
  #install.packages("Matrix")
  #library(Matrix)
  
  X = array(rep(0, p*p*n), c(p,p,n)) #  pxp adjacent matrix sequence of length n
  U = V =  array(0,dim = c(p,p,n-1))
  U_full = V_full =  array(0,dim = c(p,p,n))
  D = matrix(rep(0, n*p), nrow=n)  # Node degress
  
  ThetaM=theta%*%t(theta)
  EtaM=eta%*%t(eta)
  
  
  # Set initial adjacent matrix X[1,,] 
  zeroOne=c(0,1); zeroOne=as.vector(zeroOne)
  if (is.na(initial)){
    for(i in 1:(p-1)) for(j in (i+1):p){ 
      pb=ThetaM[i,j]/(EtaM[i,j]+ThetaM[i,j])
      X[i,j,1]=sample(zeroOne, 1, prob=c(1-pb, pb)); X[j,i,1]=X[i,j,1]}
    for(i in 1:p) D[1,i] = sum(X[i,,1])/(p-1) # Node degrees at time t=1
  }else{
    for(i in 1:(p-1)) for(j in (i+1):p){ 
      pb=initial
      X[i,j,1]=sample(zeroOne, 1, prob=c(1-pb, pb)); X[j,i,1]=X[i,j,1]}
    for(i in 1:p) D[1,i] = sum(X[i,,1])/(p-1) # Node degrees at time t=1
  }
  
  
  # Generate adjacent matrices for t=2, ... N
  for(t in 2:n){
    # Set transition probability matrix first, then generate X[t,,]
    for(i in 1:(p-1)) for(j in (i+1):p){ 
      tmp = U[i,j,t-1] = U[j,i,t-1] = sum(X[i,, t-1]*X[j,,t-1])/(p-1) # No. of common friends
      tmp1 = V[i,j,t-1] = V[j,i,t-1] = (sum(X[i,,t-1] + X[j,,t-1])- 2* X[i,j,t-1])/(p-1) - 2*tmp       # No. of uncommon friends
      if(X[i,j,t-1]==0) { 
        pb=ThetaM[i,j]*exp(a1*tmp)/(exp(a1*tmp)+exp(b1*tmp1)+1)
        X[i,j,t]=sample(zeroOne, 1, prob=c(1-pb, pb)); X[j,i,t]=X[i,j,t]
      } else {
        pb=EtaM[i,j]*exp(b2*tmp1)/(1+exp(a2*tmp)+exp(b2*tmp1))	# Prob of 1 to 0
        X[i,j,t]=sample(zeroOne, 1, prob=c(pb, 1-pb)); X[j,i,t]=X[i,j,t]
      }
    }
    # print(X[t,,])
    for(i in 1:p) D[t,i] = sum(X[i,,t])/(p-1) # Node degrees at time t
  }
  
  U_full[,,1:(n-1)] = U
  V_full[,,1:(n-1)] = V
  
  U_full[i,j,n] = U_full[j,i,n] = sum(X[i,, n]*X[j,,n])/(p-1) # No. of common friends
  V_full[i,j,n] = V_full[j,i,n] = ((sum(X[i,,n]+X[j,,n])-2*tmp) - 2* X[i,j,n])/(p-1) # No. of uncommon friends
  
  
  res = list()
  res$X = X
  res$D = D
  res$thetaM = ThetaM
  res$etaM = EtaM
  res$U = U
  res$V = V
  res$U_full = U_full
  res$V_full = V_full
  
  
  return(res)
}


simuM_dd = function(p, n, theta, eta, a , b) {
  # p: No. of nodes
  # n: No. of observed networks
  # a: coefficient of alpha
  # b: coefficient of beta
  
  #install.packages("Matrix")
  #library(Matrix)
  
  X = array(rep(0, p*p*n), c(p,p,n)) #  pxp adjacent matrix sequence of length n
  D = matrix(rep(0, n*p), nrow=n)  # Node degrees
  D_all = rep(0,n)
  
  ThetaM=theta%*%t(theta)
  EtaM=eta%*%t(eta)
  
   
  
  # Set initial adjacent matrix X[1,,] with the first 3 nodes as Harbour nodes
  zeroOne=c(0,1); zeroOne=as.vector(zeroOne)
  for(i in 1:(p-1)) for(j in (i+1):p){ 
    pb=ThetaM[i,j]/(EtaM[i,j]+ThetaM[i,j])
    X[i,j,1]=sample(zeroOne, 1, prob=c(1-pb, pb)); X[j,i,1]=X[i,j,1]}
  for(i in 1:p) D[1,i] = sum(X[i,,1]) # Node degrees at time t=1
  D_all[1] = sum(D[1,])/2 
  
  # Generate adjacent matrices for t=2, ... N
  for(t in 2:n){
    # Set transition probability matrix first, then generate X[t,,]
    for(i in 1:(p-1)) for(j in (i+1):p){ 
      tmp = a[1]* D_all[t-1] +a[2]*D[t-1,i] +a[3]*D[t-1,j]
      tmp1 = b[1]*( p*(p-1)/2 - D_all[t-1]) +b[2]*(p - 1 - D[t-1,i]) +b[3]*(p - 1 - D[t-1,j])
      if(X[i,j,t-1]==0) { 
        if (is.infinite(exp(tmp))){
          pb = 1
        }else{
          pb=ThetaM[i,j]* exp(tmp)/(1+exp(tmp))
        }
        X[i,j,t]=sample(zeroOne, 1, prob=c(1-pb, pb)); X[j,i,t]=X[i,j,t]
      } else {
        if (is.infinite(exp(tmp1))){
          pb = 1
        }else{
          pb= EtaM[i,j]*exp(tmp1)/(1+exp(tmp1))
        }
        	# Prob of 1 to 0
        X[i,j,t]=sample(zeroOne, 1, prob=c(pb, 1-pb)); X[j,i,t]=X[i,j,t]
      }
    }
    # print(X[t,,])
    for(i in 1:p) D[t,i] = sum(X[i,,t]) # Node degrees at time t
    D_all[t] = sum(D[t,])/2
  }
  

  
  
  res = list()
  res$X = X
  res$D = D
  
  
  return(res)
}






