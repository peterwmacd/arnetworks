logadj = function(x){
  x = ifelse(x<=0,0.0001,x)
  res = log(x)
  return(res)

}


fg_array = function(fg, global, stats){
  dm = dim(stats)
  p = dm[1]
  n = dm[3]
  res = array(fg(global, stats), dim = c(p,p,n))
  return(res)
}


### 1.1  part of global log likelihood function:  for alpha  or beta separately
logl_part = function(A, B, fg, global, stats, localMat){
  dm = dim(stats)
  p = dm[1]
  n = dm[3]

  diag(localMat) = 0
  tmp = fg_array(fg, global, stats)
  tmp = array(apply(tmp, c(3), function(x) localMat*x),dim = c(p,p,n))
  res =  sum(A*logadj(tmp) + B*logadj(1-tmp))
  return(res)
}

### 1.2  global log likelihood function
logl = function(A1, B1, A2, B2, fij, gAlphaVal, statsAlpha, xiMat, gij, gBetaVal, statsBeta, etaMat){
  res =  logl_part(A1, B1, fij, gAlphaVal, statsAlpha, xiMat) + logl_part(A2, B2, gij, gBetaVal, statsBeta, etaMat)
  return(res)
}

### 1.3  local log likelihood function for each (i,j)
logl_local = function(Aij, Bij, fg, global, statsij, localij){
  tmp = localij*fg(global, statsij)
  res =  sum(Aij*logadj(tmp) + Bij*logadj(1-tmp))
  return(res)
}

### 1.3  local log likelihood function for each (i,j)
logl_locali = function(A, B, fg, global, stats, localVec, i){
  dm = dim(stats)
  p = dm[1]
  n = dm[3]
  
  tmp = array(fg(global, stats[i,,,,drop = FALSE]),dim = c(p,n))
  tmp =  localVec[i]*localVec[-i]*tmp[-i,]
  res =  sum(A[i,-i,]*logadj(tmp) + B[i,-i,]*logadj(1-tmp))
  return(res)
}




### 2.1  Global log likelihood function
# shrGPrm: the number of shared global pamameters; if shrGPrm = 2, we have gAlphaVal[1] = gBetaVal[1]  and gAlphaVal[2] = gBetaVal[2].
# updateAlpha: whether the target is alpha-only global parameter
globalMLE = function(par, ix, A1, B1, A2, B2, fij, gAlphaVal, statsAlpha, xiMat, gij, gBetaVal, statsBeta, etaMat, shrGPrm, updateAlpha = TRUE){
  dm = dim(A1)
  p = dm[1]
  n = dm[3]

  da = length(gAlphaVal)
  db = length(gBetaVal)

  if ( da < shrGPrm | db < shrGPrm){
    warning("Please check the length of global parameters.")
  }

  res = 0
  if (ix <= shrGPrm){
    gAlphaVal[ix] = par
    gBetaVal[ix] = par
    res = - logl(A1, B1, A2, B2, fij, gAlphaVal, statsAlpha, xiMat, gij, gBetaVal, statsBeta, etaMat)
  }else if (updateAlpha){
    gAlphaVal[ix] = par
    res = - logl_part(A1, B1, fij, gAlphaVal, statsAlpha, xiMat)
  }else{
    gBetaVal[ix] = par
    res = - logl_part(A2, B2, gij, gBetaVal, statsBeta, etaMat)
  }

  res

}



### 2.2  Local log likelihood function for xi_ij and eta_ij
localMLE = function(par, Aij, Bij, fg, global, statsij){
  res = - logl_local(Aij, Bij, fg, global, statsij, par)
  res
}

locfr = function(logLoc, locM){
  # Given the values of xi matrix and eta matrix
  # Find the estimated values for xi and eta
  # Jiang et al (2024)
  logLocM = outer(logLoc,logLoc, "+")
  logLocM[lower.tri(logLocM,diag = T)] = 0
  sum(exp(logLocM))  - sum(logLoc * apply(locM,1,sum))
}

locEst = function(locMatE){
  # locMatE: symmetric matrix with diagonal elements being zero.
  diag(locMatE) = 0
  tmp = stats::optim(rep(1, dim(locMatE)[1]), locfr, method = 'BFGS', locM = locMatE)
  locE = exp(tmp$par)

  return(locE)
}

### 2.3  Local log likelihood function for xi_i and eta_i
localMLE_i = function(par, A, B, fg, global, stats, localVec, i){
  localVec[i] = par
  res = - logl_locali(A, B, fg, global, stats, localVec, i)
  res
}


globalMLE_group = function(par, A1, B1, A2, B2, fij, gAlphaVal, statsAlpha, xiMat, gij, gBetaVal, statsBeta, etaMat, shrGPrm, updateShare = T, updateAlpha = TRUE){
  dm = dim(A1)
  p = dm[1]
  n = dm[3]

  da = length(gAlphaVal)
  db = length(gBetaVal)

  if ( da < shrGPrm | db < shrGPrm){
    warning("Please check the length of global parameters.")
  }

  res = 0
  if (updateShare){
    gAlphaVal[1:shrGPrm] = par
    gBetaVal[1:shrGPrm]= par
    res = - logl(A1, B1, A2, B2, fij, gAlphaVal, statsAlpha, xiMat, gij, gBetaVal, statsBeta, etaMat)
  }else if (updateAlpha){
    gAlphaVal[(shrGPrm+1):da] = par
    res = - logl_part(A1, B1, fij, gAlphaVal, statsAlpha, xiMat)
  }else{
    gBetaVal[(shrGPrm+1):db] = par
    res = - logl_part(A2, B2, gij, gBetaVal, statsBeta, etaMat)
  }

  res

}
