library(doParallel)
library(Matrix)
library(MASS)
source("Step1_GenNet.R")
source("Functions.R")
library(lpSolve)
iter = 100
sim = 100
p = 50
n.seq = c(100,200,400)
burnin = 200
delta = 0.05
delta_init = delta
gamma_seq = c(0.5)
r_seq = c(0.5,0.1)
paraset = matrix(c(
  0.7,0.7,1.2,1.6,
  0.7,0.8,0.8,2,
  0.6,0.8,2,1.6,
  0.8,0.7,0.8,1.6,
  0.6,0.7,0.4,0.8
) ,
ncol = 4, byrow = T
)


initset = matrix(c(
  0.5,0.5,
  0.5,0.7,
  0.5,0.9,
  0.7,0.5,
  0.7,0.7,
  0.7,0.9,
  0.9,0.5,
  0.9,0.7,
  0.9,0.9
) ,
ncol = 2, byrow = T
)



resTable = matrix(NA,dim(paraset)[1]*length(n.seq), 4+ (1+length(gamma_seq))*4  )
resTable_all = array(NA,dim = c(dim(paraset)[1]*length(n.seq), 4+ (1+length(gamma_seq))*4, dim(initset)[1] ))
rownames(resTable) = as.character(rep(n.seq,dim(paraset)[1] ))
colnames(resTable) = rep(c('theta','eta','a','b'), 1+  (1+length(gamma_seq))) 
ind = 1



for (ix1 in 1:(dim(paraset)[1])){
  
  theta = rep(paraset[ix1,1],p)
  eta = rep(paraset[ix1,2],p)
  a = paraset[ix1,3]
  b = paraset[ix1,4]
  
  for (n in n.seq){
    
    G = 2
    Gc = 2*p
    q = G + Gc
    S_g = p*(p-1)
    S_gc = p-1
    c_g2 = 0.01*max(q *log(n*S_g)/sqrt(n*S_g),q^{3/2} *(log(n*S_g))^{3/2}/sqrt(n)/S_g)
    c_gc2 = 0.01*max(q *log(n * S_gc)/sqrt(n * S_gc),q^{3/2} *(log(n*S_gc))^{3/2}/sqrt(n)/S_gc)
    
    delta_n_sqrt = sqrt(max(c_g2, c_gc2))
    # gamma \lessim delta_n_sqrt
    
    
    
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
      ab0 = ab1 = c(1,1)

      
      res0 = array(0,dim = c(4,1+length(gamma_seq),dim(initset)[1]))
      for (ix2 in 1:dim(initset)[1]){
        thetaE = theta0 = rep(initset[ix2,1],p)
        eta0 = etaE = rep(initset[ix2,2],p)
        
        tmp0 = optim(ab1, globalMLE_ab, gr = grr_globalMLE_ab, method = "L-BFGS-B", lower = c(0.01,0.01), A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, thetavec = thetaE, etavec = etaE)
        ab1 = tmp0$par
        
        ##### Initial value: ab1, thetaE and etaE
        ##### Now we apply the global MLE for (a,b) and local MLE for thetai and etai
        thetamax = apply((1 +exp(ab1[1]*U) + exp(ab1[2]*V) )/exp(ab1[1]*U),c(1,2),min)
        etamax =  apply((1 +exp(ab1[2]*V) + exp(ab1[1]*U) )/exp(ab1[2]*V),c(1,2),min)
        
        for (i in 1:p){
          tmp0 = optim(c(thetaE[i],etaE[i]), localMLE_thetaeta,  method = 'L-BFGS-B', A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, a = ab1[1], b = ab1[2], thetavec = thetaE, etavec = etaE, i = i, lower  = c(0.01,0.01), upper = c(min(thetamax[i,-i]/thetaE[-i]),min(etamax[i,-i]/etaE[-i])))$par
          thetaE[i] = tmp0[1]
          etaE[i] = tmp0[2]
        }
        
        
        
        
        
        
        
        
        res = vector(mode='list', length=(1+length(gamma_seq)))
        
        
        res[[1]]$theta.abs.err = mean(abs(theta-thetaE))
        res[[1]]$eta.abs.err = mean(abs(eta-etaE))
        res[[1]]$a.abs.err = abs(ab1[1]-a)
        res[[1]]$b.abs.err = abs(ab1[2]-b)
        
        
        ab10 = ab1
        thetaE0 = thetaE
        etaE0 = etaE
        
        
        
        ind_gr =1
        
        ### Global for a and b
        for (gamma_c in gamma_seq){
          gamma_theta = gamma_c * delta_n_sqrt
          gamma_a = 0.01 * delta_n_sqrt
          gamma_b = 0.01 * delta_n_sqrt
          gamma_eta = 0.5 * delta_n_sqrt
          ab2 = ab1 = ab10
          thetaR = thetaE = thetaE0
          etaR = etaE = etaE0
          
          for (r_tilda in r_seq){
            
            gArray = array(0, dim = c(2+p*2,2+p*2,p ))
            
            for (i in 1:p){
              for (j in (1:p)[-i]){
                A1ij = A1[i,j,]
                B1ij = B1[i,j,]
                A2ij = A2[i,j,]
                B2ij = B2[i,j,]
                Uij= U[i,j,]
                Vij= V[i,j,]
                thetai = thetaE[i]
                thetaj = thetaE[j]
                etai = etaE[i]
                etaj = etaE[j]
                gArray[1,1,i] = gArray[1,1,i] + dl2dada_ab(A1ij, B1ij, A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], thetai, thetaj, etai, etaj)
                gArray[1,2,i] = gArray[2,1,i] = gArray[1,2,i] + dl2dadb_ab(A1ij, B1ij, A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], thetai, thetaj, etai, etaj)
                gArray[2,2,i] = gArray[2,2,i] + dl2dbdb_ab(A1ij, B1ij, A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], thetai, thetaj, etai, etaj)
                
                gArray[1,2+i,i] =  gArray[2+i,1,i] = gArray[1,2+i,i] + dl2dadthetai(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], thetai, thetaj)
                gArray[1,2+j,i] =  gArray[2+j,1,i] = gArray[1,2+j,i] + dl2dadthetai(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], thetaj, thetai)
                gArray[2,2+i,i] =  gArray[2+i,2,i] = gArray[2,2+i,i] + dl2dbdthetai(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], thetai, thetaj)
                gArray[2,2+j,i] =  gArray[2+j,2,i] = gArray[2,2+j,i] + dl2dbdthetai(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], thetaj, thetai)
                
                gArray[1,2+p+i,i] = gArray[2+p+i,1,i] = gArray[1,2+p+i,i] + dl2dadetai_ab(A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], etai, etaj)
                gArray[1,2+p+j,i] = gArray[2+p+j,1,i] = gArray[1,2+p+j,i] + dl2dadetai_ab(A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], etaj, etai)
                gArray[2,2+p+i,i] = gArray[2+p+i,2,i] = gArray[2,2+p+i,i] + dl2dbdetai_ab(A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], etai, etaj)
                gArray[2,2+p+j,i] = gArray[2+p+j,2,i] = gArray[2,2+p+j,i] + dl2dbdetai_ab(A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], etaj, etai)
                
                gArray[2+i,2+i,i] = gArray[2+i,2+i,i] + dl2dthetaidthetai(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], thetai, thetaj)
                gArray[2+j,2+j,i] = gArray[2+j,2+j,i] + dl2dthetaidthetai(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], thetaj, thetaj)
                gArray[2+i,2+j,i] = gArray[2+j,2+i,i] =  gArray[2+i,2+j,i] + dl2dthetaidthetaj(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], thetai, thetaj)
                
                gArray[2+p+i,2+p+i,i] = gArray[2+p+i,2+p+i,i] + dl2dthetaidthetai(A2ij, B2ij, Vij, Uij, ab1[2], ab1[1], etai, etaj)
                gArray[2+p+j,2+p+j,i] =  gArray[2+p+j,2+p+j,i] +dl2dthetaidthetai(A2ij, B2ij, Vij, Uij, ab1[2], ab1[1], etaj, etai)
                gArray[2+p+i,2+p+j,i] = gArray[2+p+j,2+p+i,i] = gArray[2+p+i,2+p+j,i] +dl2dthetaidthetaj(A2ij, B2ij, Vij, Uij, ab1[2], ab1[1], etai, etaj)
                
              }
            }
            
            g_local = gArray/(p-1)
            g_global = apply(gArray, c(1,2),sum)/(p*(p-1))
            g_all = array(c(g_global, g_global, g_local, g_local),dim = c(2+2*p,2+2*p,2*p+2))
            
            El = diag(1, 2+2*p)
            Al = array(0, dim = c(2+2*p,2+2*p)) #order: a,b, theta1, \dots, thetap, eta1, \dots, etap
            
            Al[,1] = alSearch(g_all[,,1], El[,1], gamma_a)
            
            
            
            Al[,2] = alSearch(g_all[,,2], El[,2], gamma_b)
            
            
            for (i in 3:(2 + p)){
              Al[,i] = alSearch(g_all[,,i], El[,i], gamma_theta)
            }
            
            for (i in (p+3):(2 + 2*p)){
              Al[,i] = alSearch(g_all[,,i], El[,i], gamma_eta)
            }
            
            
            ab2[1] = optim(ab1[1],  globalMLE_refine_a,  method = 'L-BFGS-B', al = Al[,1], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, ab = ab1, thetavec = thetaE, etavec = etaE, lower = max(0.01, ab1[1]-r_tilda), upper = ab1[1]+r_tilda)$par
            ab2[2] = optim(ab1[2],  globalMLE_refine_b,  method = 'L-BFGS-B', al = Al[,2], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, ab = ab1, thetavec = thetaE, etavec = etaE, lower = max(0.01, ab1[2]-r_tilda), upper = ab1[2]+r_tilda)$par
            
            thetamax = apply((1 +exp(ab1[1]*U) + exp(ab1[2]*V) )/exp(ab1[1]*U),c(1,2),min)
            etamax =  apply((1 +exp(ab1[2]*V) + exp(ab1[1]*U) )/exp(ab1[2]*V),c(1,2),min)
            for (i in 1:p){
              thetaR[i] = optim(thetaE[i],  localMLE_refine_theta, grr_localMLE_refine_theta,  method = 'L-BFGS-B', al = Al[,2+i], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, ab = ab1, thetavec = thetaE, etavec = etaE, i = i, lower = max(0.01, thetaE[i]-r_tilda), upper = min(c(thetamax[i,-i]/thetaE[-i], thetaE[i]+r_tilda)))$par
              etaR[i] = optim(etaE[i],  localMLE_refine_eta, grr_localMLE_refine_eta,  method = 'L-BFGS-B', al = Al[,2+p+i], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, ab = ab1, thetavec = thetaE, etavec = etaE, i = i, lower = max(0.01, etaE[i]-r_tilda), upper = min(c(etamax[i,-i]/etaE[-i], etaE[i]+r_tilda)))$par
            }
            
            ab1 = ab2
            thetaE = thetaR
            etaE = etaR
          }
         
          res[[1+ind_gr]]$theta.abs.err = mean(abs(theta-thetaE))
          res[[1+ind_gr]]$eta.abs.err = mean(abs(eta-etaE))
          res[[1+ind_gr]]$a.abs.err = abs(ab1[1]-a)
          res[[1+ind_gr]]$b.abs.err = abs(ab1[2]-b)
          ind_gr = ind_gr+1 
          
          
          
        }
        
        res0[,,ix2] = matrix(c(sapply(res, function(x) x$theta.abs.err, simplify = "array"),
                             sapply(res, function(x) x$eta.abs.err, simplify = "array"),
                             sapply(res, function(x) x$a.abs.err, simplify = "array"),
                             sapply(res, function(x) x$b.abs.err, simplify = "array")), byrow = T ,nrow = 4)
        
      }
     
      
      return(res0)
      
      
      
    }
    
    save(out, file = paste0("est2_n", n, "_a", a, "_b", b, "_theta", theta[1], "_eta", eta[1], ".RData"))
    
    stopCluster(cl)
    
    for (ix2 in 1:dim(initset)[1]){
      resTable_all[ind,,ix2] = c(c(theta[1],eta[1],a,b),apply( matrix(do.call(c,lapply(out, function(x) x[,,ix2])),ncol = length(out)),1, function(x) mean(x,na.rm = T)))
      
    }

    resTable[ind,] = c(c(theta[1],eta[1],a,b),apply( matrix(do.call(c,lapply(out, function(x) x)),ncol = (dim(initset)[1])*length(out)),1, function(x) mean(x,na.rm = T)))
    
    ind = ind+1
    
    save(resTable_all, file = paste0("resArray_",p,gamma_seq[1],".RData"))
    save(resTable, file = paste0("resTable_",p,gamma_seq[1],".RData"))
    write.csv(resTable, file = paste0("resTable_",p,gamma_seq[1],".csv"))
    
  }
  
  
}











