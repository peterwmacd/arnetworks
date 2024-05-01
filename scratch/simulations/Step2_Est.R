library(doParallel)
library(Matrix)
library(MASS)
source("Step1_GenNet.R")
source("Functions.R")
library(lpSolve)
iter = 100
sim = 100
p = 100
n.seq = c(10,20,50,100,200)
burnin = 200
delta_init = 1
delta = 0.0001
paraset = matrix(c(
                   0.9,0.9,0.4,0.8,
                   0.7,0.9,0.8,0.8,
                   0.9,0.9,1.2,1.6,
                   0.7,0.9,0.8,2,
                   0.5,0.9,2,1.6,
                   0.3,0.9,0.8,0.4,
                   0.7,0.7,0.4,2,
                   0.5,0.7,2,2,
                   0.9,0.9,1.6,1.6,
                   0.3,0.3,0.8,1.2,
                   0.7,0.7,1.6,1.6,
                   0.1,0.7,0.8,0.8,
                   0.5,0.5,0.8,0.8
                   
) ,
ncol = 4, byrow = T
)



resTable = matrix(NA,dim(paraset)[1]*length(n.seq), 4+16  )

rownames(resTable) = as.character(rep(n.seq,dim(paraset)[1] ))
colnames(resTable) = c('theta','eta','a','b',  
                       "theta.mean","theta.sd","eta.mean","eta.sd","theta.abs.mean","theta.abs.sd","eta.abs.mean","eta.abs.sd",
                       "a.mean","a.sd","b.mean","b.sd", "a.abs.mean","a.abs.sd","b.abs.mean","b.abs.sd") 
ind = 1



for (ix1 in 1:(dim(paraset)[1])){
  
  theta = rep(paraset[ix1,1],p)
  eta = rep(paraset[ix1,2],p)
  a = paraset[ix1,3]
  b = paraset[ix1,4]
  
  for (n in n.seq){
    
    ############ Replicates
    cl = makeCluster(50)
    registerDoParallel(cl)
    
    out = foreach(try = 1:sim, .packages = c("Matrix","MASS",'lpSolve'), .errorhandling = "remove") %dopar%{
      set.seed(try+1000)
      data = simuM(p, n+burnin, theta, eta,  a, b, a, b)
      X = data$X[,, burnin+1: n]
      U = data$U[,, burnin+1: (n-1)]
      V = data$V[,, burnin+1: (n-1)]
      
      
      
      A1 = X[,,2:n]*(1 - X[,,2:n-1])
      B1 = (1 - X[,,2:n])*(1 - X[,,2:n-1])
      A2 = (1 - X[,,2:n])*( X[,,2:n-1])
      B2 = X[,,2:n]*(X[,,2:n-1])
      
      #Initialization:
      ab1 = c(a+delta_init,b+delta_init)
      thetaME = pmax(pmin(apply((1 +exp(ab1[1]*U) + exp(ab1[2]*V) )/exp(ab1[1]*U),c(1,2),min),theta+delta_init)-0.2,0)
      etaME =  pmax(pmin(apply((1 +exp(ab1[2]*V) + exp(ab1[1]*U) )/exp(ab1[2]*V),c(1,2),min),eta+delta_init)-0.2,0)
      thetaE = thetaEst(thetaME) 
      etaE =  thetaEst(etaME)
      
      tmp0 = optim(ab1, globalMLE_ab, gr = grr_globalMLE_ab, method = "L-BFGS-B", lower = c(0,0), A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, thetavec = thetaE, etavec = etaE)
      ab1 = tmp0$par
      fn1 = tmp0$value
      
      for (it in 1: iter){
        thetamax = apply((1 +exp(ab1[1]*U) + exp(ab1[2]*V) )/exp(ab1[1]*U),c(1,2),min)
        etamax =  apply((1 +exp(ab1[2]*V) + exp(ab1[1]*U) )/exp(ab1[2]*V),c(1,2),min)
        for (i in 1:(p-1)){
          for (j in (i+1):p){
            thetaij = optim(thetaME[i,j], localMLE.init, method = 'L-BFGS-B', Aij = A1[i,j,], Bij = B1[i,j,], Uij=U[i,j,], Vij=V[i,j,], ab = ab1, lower  = c(0), upper = c(thetamax[i,j]))$par
            thetaME[j,i] = thetaME[i,j] = thetaij
            
            etaij = optim(etaME[i,j], localMLE.init, method = 'L-BFGS-B', Aij = A2[i,j,], Bij = B2[i,j,], Uij=V[i,j,], Vij=U[i,j,], ab = c(ab1[2],ab1[1]), lower  = c(0), upper = c(etamax[i,j]))$par
            etaME[j,i] = etaME[i,j] = etaij
            
          }
        }
        
        thetaE = thetaEst(thetaME) 
        etaE = thetaEst(etaME)
        thetaME = outer(thetaE, thetaE)
        etaME = outer(etaE, etaE)
        
        
        tmp2 = optim(ab1, globalMLE_ab, gr = grr_globalMLE_ab, method = "L-BFGS-B", lower = c(0,0), A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, thetavec = thetaE, etavec = etaE)
        ab2 = tmp2$par
        fn2 = tmp2$value
        
        if(mean(abs(ab2 -ab1))<delta| fn1<=fn2) {
          break
        }else{
          ab1 = ab2
          fn1 = fn2
        }
      }
      
      
      ##### Initial value: ab1, thetaE and etaE
      ##### Now we apply the global MLE for (a,b) and local MLE for thetai and etai
      
      
      for (it in 1: iter){
        thetamax = apply((1 +exp(ab1[1]*U) + exp(ab1[2]*V) )/exp(ab1[1]*U),c(1,2),min)
        etamax =  apply((1 +exp(ab1[2]*V) + exp(ab1[1]*U) )/exp(ab1[2]*V),c(1,2),min)
        
        for (i in 1:p){
          thetaE[i] = optim(thetaE[i],  localMLE,  method = 'L-BFGS-B', Ai = A1[i,-i,] , Bi = B1[i,-i,], Ui = U[i,-i,], Vi=V[i,-i,], ab=ab1, thetavec_ic= thetaE[-i], lower  = c(0), upper = min(thetamax[i,-i]/thetaE[-i]))$par
          etaE[i] = optim(etaE[i],  localMLE,  method = 'L-BFGS-B', Ai = A2[i,-i,] , Bi = B2[i,-i,], Ui = V[i,-i,], Vi=U[i,-i,], ab=c(ab1[2],ab1[1]), thetavec_ic= etaE[-i], lower  = c(0), upper = min(etamax[i,-i]/etaE[-i]))$par
          }

        tmp2 = optim(ab1, globalMLE_ab, gr = grr_globalMLE_ab, method = "L-BFGS-B", lower = c(0,0), A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, thetavec = thetaE, etavec = etaE)
        ab2 = tmp2$par
        fn2 = tmp2$value
        
        if(mean(abs(ab2 -ab1))<delta| fn1<=fn2) {
          break
        }else{
          ab1 = ab2
          fn1 = fn2
        }
      }
      

      res = list()
      
      res$theta.err = mean((theta-thetaE)^2)
      res$eta.err = mean((eta-etaE)^2)
      res$theta.abs.err = mean(abs(theta-thetaE))
      res$eta.abs.err = mean(abs(eta-etaE))
      res$a.err = (ab1[1]-a)^2
      res$b.err = (ab1[2]-b)^2
      res$a.abs.err = abs(ab1[1]-a)
      res$b.abs.err = abs(ab1[2]-b)
      
      
      
      return(res)
    }
    
    save(out, file = paste0("est2_n", n, "_a", a, "_b", b, "_theta", theta[1], "_eta", eta[1], ".RData"))
    
    stopCluster(cl)
    
    resTable[ind,1:4] = c(theta[1],eta[1],a,b)
    
    for (ix2 in 0:7){
      resTable[ind,4+2*ix2+1] = mean(do.call(c,lapply(out, function(x) x[[ix2+1]])),na.rm = T)
      resTable[ind,4+2*ix2+2] = sd(do.call(c,lapply(out, function(x) x[[ix2+1]])),na.rm = T)
      
    }
    ind = ind+1
    
    save(resTable, file = paste0("resTable2_set1_",p,".RData"))
    write.csv(resTable, file = paste0("resTable2_set1_",p, ".csv"))
    
  }
  
  
}

save(resTable, file = paste0("resTable2_set1_",p, ".RData"))
write.csv(resTable, file = paste0("resTable2_set1_",p,".csv"))










