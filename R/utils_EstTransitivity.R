# Utility functions for estimating transitivity model
# note that these functions have '_et' appended so that they are correctly
# referenced within estTransitivity

### 1 loglikelihood function for different (a, b) pair
logl_et = function(Aij, Bij, Uij, Vij, a, b, xii, xij){
  n = length(Aij)
  res = sum(Aij * logadj(xii * xij * exp(a* Uij)) +
              Bij * logadj(1+ exp(b* Vij) +  (1-xii * xij)* exp(a * Uij)) -
              (Aij + Bij)* logadj(1+ exp(b* Vij) +  exp(a * Uij)) )/(n-1)
  return (res)
}

### 1.2  Global loglikelihood function for different (a, b) pair
globalMLE_et = function(par, A, B, U, V, xivec){
  a = par[1]
  b = par[2]
  dm = dim(A)
  n = dm[3]
  p = dm[1]

  res = 0
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      Aij = A[i,j,]
      Bij = B[i,j,]
      Uij= U[i,j,]
      Vij= V[i,j,]
      xii = xivec[i]
      xij = xivec[j]
      res = res - logl_et(Aij, Bij, Uij, Vij, a, b, xii, xij)
    }
  }
  res
}

grr_globalMLE_et = function(par, A, B, U, V, xivec){
  a = par[1]
  b = par[2]
  dm = dim(A)
  n = dm[3]
  p = dm[1]


  vec1 = rep(0, 2)

  for (i in 1:(p-1)){
    for (j in (i+1):p){
      Aij = A[i,j,]
      Bij = B[i,j,]
      Uij= U[i,j,]
      Vij= V[i,j,]
      xii = xivec[i]
      xij = xivec[j]
      vec1[1] = vec1[1] - dlda_et(Aij, Bij, Uij, Vij, a, b, xii, xij)
      vec1[2] = vec1[2] - dldb_et(Aij, Bij, Uij, Vij, a, b, xii, xij)
    }
  }
  vec1

}


### 2. loglikelihood function when alpha and beta share the same pair of (a,b).
logl_ab_et = function(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj){
  return(logl_et(A1ij, B1ij, Uij, Vij, a, b, xii, xij) + logl_et(A2ij, B2ij, Vij, Uij, b, a, etai, etaj))
}

### 2.2  Global loglikelihood function for seperate a, b
globalMLE_ab_et = function(par, A1, B1, A2, B2, U, V, xivec, etavec){
  a = par[1]
  b = par[2]
  dm = dim(A1)
  n = dm[3]
  p = dm[1]

  res = 0
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      A1ij = A1[i,j,]
      B1ij = B1[i,j,]
      A2ij = A2[i,j,]
      B2ij = B2[i,j,]
      Uij= U[i,j,]
      Vij= V[i,j,]
      xii = xivec[i]
      xij = xivec[j]
      etai = etavec[i]
      etaj = etavec[j]
      res = res - logl_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai,etaj)
    }
  }
  res
}

grr_globalMLE_ab_et = function(par, A1, B1, A2, B2, U, V, xivec, etavec){
  a = par[1]
  b = par[2]
  dm = dim(A1)
  n = dm[3]
  p = dm[1]


  vec1 = rep(0, 2)

  for (i in 1:(p-1)){
    for (j in (i+1):p){
      A1ij = A1[i,j,]
      B1ij = B1[i,j,]
      A2ij = A2[i,j,]
      B2ij = B2[i,j,]
      Uij= U[i,j,]
      Vij= V[i,j,]
      xii = xivec[i]
      xij = xivec[j]
      etai = etavec[i]
      etaj = etavec[j]
      vec1[1] = vec1[1] - dlda_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai,etaj)
      vec1[2] = vec1[2] - dldb_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai,etaj)
    }
  }
  vec1

}



### 3  Local loglikelihood function for xi_ij and eta_ij
# Example for xi_ij
localMLE_init_et = function(par, Aij, Bij, Uij, Vij, ab){
  xiij = par
  a = ab[1]
  b = ab[2]

  n = length(Aij)
  res = -sum(Aij * logadj(xiij* exp(a* Uij)) +
               Bij * logadj(1+ exp(b* Vij) +  (1-xiij)* exp(a * Uij)) -
               (Aij + Bij)* logadj(1+ exp(b* Vij) +  exp(a * Uij)) )/(n-1)

  res
}




### 4  Local loglikelihood function for xi_i and eta_i
localMLE_et = function(par, Ai, Bi, Ui, Vi, ab, xivec_ic){
  xii = par
  a = ab[1]
  b = ab[2]
  dm = dim(Ai)
  n = dm[2]
  p = dm[1]

  res = 0

  for (j in 1:(p)){
    Aij =Ai[j,]
    Bij =Bi[j,]
    Uij=Ui[j,]
    Vij=Vi[j,]
    xij = xivec_ic[j]
    res = res - logl_et(Aij, Bij, Uij, Vij, a, b, xii, xij)
  }
  res
}


#######-------------------------- First direvatives of Loglikelihood function --------------
dlda_et = function(Aij, Bij, Uij, Vij, a, b, xii, xij){
  n = length(Aij)
  tmp = 1+ exp(b* Vij) +  (1-xii*xij)* exp(a * Uij)
  tmp =  ifelse(tmp==0,0.001,tmp)

  res = sum(Aij* Uij +
              Bij * (1- xii*xij)* exp(a* Uij) * Uij / (tmp) -
              (Aij + Bij)*exp(a* Uij) * Uij/(1+ exp(b* Vij) +  exp(a * Uij)) )/(n-1)

  return (res)

}

dldb_et = function(Aij, Bij, Uij, Vij, a, b, xii, xij){
  n = length(Aij)
  tmp = 1+ exp(b* Vij) +  (1-xii*xij)* exp(a * Uij)
  tmp =  ifelse(tmp==0,0.001,tmp)
  res = sum(Bij * exp(b* Vij) * Vij / (tmp) -
              (Aij + Bij)*exp(b* Vij) * Vij/(1+ exp(b* Vij) +  exp(a * Uij)) )/(n-1)

  return (res)

}

dlda_ab_et = function(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj){
  return(dlda_et(A1ij, B1ij, Uij, Vij, a, b, xii, xij) + dldb_et(A2ij, B2ij, Vij, Uij, b, a, etai, etaj))
}

dldb_ab_et = function(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj){
  return(dldb_et(A1ij, B1ij, Uij, Vij, a, b, xii, xij) + dlda_et(A2ij, B2ij, Vij, Uij, b, a, etai, etaj))
}


dldxii_et = function(Aij, Bij, Uij, Vij, a, b, xii, xij){
  n = length(Aij)
  tmp = 1+ exp(b* Vij) +  (1-xii*xij)* exp(a * Uij)
  tmp =  ifelse(tmp==0,0.001,tmp)
  res = sum( Aij/xii  -
               Bij * exp(a* Uij) * xij/(tmp))/(n-1)

  return (res)

}

###--------------Second derivatives of Loglikelihood function------------------###

dl2dada_et = function(Aij, Bij, Uij, Vij, a, b, xii, xij){
  n = length(Aij)

  tmp = 1+ exp(b* Vij) +  (1-xii*xij)* exp(a * Uij)
  tmp =  ifelse(tmp==0,0.001,tmp)

  tmp1 = (1-xii * xij)* exp(a* Uij) * Uij / (tmp)

  tmp2 = exp(a* Uij) * Uij/(1+ exp(b* Vij) +  exp(a * Uij))

  res = sum(Bij * (Uij * tmp1 -tmp1^2) -
              (Aij + Bij)*(Uij * tmp2 - tmp2^2))/(n-1)

  return (res)

}

dl2dbdb_et = function(Aij, Bij, Uij, Vij, a, b, xii, xij){
  n = length(Aij)

  tmp = 1+ exp(b* Vij) +  (1-xii*xij)* exp(a * Uij)
  tmp =  ifelse(tmp==0,0.001,tmp)

  tmp1 = exp(b* Vij) *Vij / (tmp)

  tmp2 = exp(b* Vij) *Vij/(1+ exp(b* Vij) +  exp(a * Uij))

  res = sum(Bij * (Vij * tmp1 -tmp1^2) -
              (Aij + Bij)*(Vij * tmp2 - tmp2^2) )/(n-1)

  return (res)

}


dl2dada_ab_et = function(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj){
  return(dl2dada_et(A1ij, B1ij,Uij, Vij, a, b, xii, xij) + dl2dbdb_et(A2ij, B2ij,Vij, Uij, b, a, etai, etaj))
}

dl2dbdb_ab_et = function(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj){
  return(dl2dbdb_et(A1ij, B1ij,Uij, Vij, a, b, xii, xij) + dl2dada_et(A2ij, B2ij,Vij, Uij, b, a, etai, etaj))
}



dl2dadb_et = function(Aij, Bij, Uij, Vij, a, b, xii, xij){
  n = length(Aij)
  tmp = 1+ exp(b* Vij) +  (1-xii*xij)* exp(a * Uij)
  tmp =  ifelse(tmp==0,0.001,tmp)

  res = sum((Aij + Bij)* exp(a* Uij) * Uij * exp(b* Vij) * Vij/ (1+ exp(b* Vij) +  exp(a * Uij))^2 -
              Bij * exp(a* Uij) * Uij * exp(b* Vij) * Vij *(1-xii * xij)/(tmp)^2 )/(n-1)

  return (res)

}

dl2dadb_ab_et = function(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj){
  return(dl2dadb_et(A1ij, B1ij,Uij, Vij, a, b, xii, xij) + dl2dadb_et(A2ij, B2ij,Vij, Uij, b, a, etai, etaj))
}



dl2dadxii_et = function(Aij, Bij, Uij, Vij, a, b, xii, xij){
  n = length(Aij)
  tmp = 1+ exp(b* Vij) +  (1-xii*xij)* exp(a * Uij)
  tmp =  ifelse(tmp==0,0.001,tmp)

  res = sum( Bij * exp(a* Uij)^2 * Uij  *(1-xii * xij)/(tmp)^2 *xij -
               Bij * exp(a* Uij) * Uij /(tmp)*xij)/(n-1)

  return (res)

}



dl2dbdxii_et = function(Aij, Bij, Uij, Vij, a, b, xii, xij){
  n = length(Aij)
  tmp = 1+ exp(b* Vij) +  (1-xii*xij)* exp(a * Uij)
  tmp =  ifelse(tmp==0,0.001,tmp)

  res = sum( Bij * exp(a* Uij) * Vij * exp(b* Vij) /(tmp)^2 * xij)/(n-1)

  return (res)

}


dl2dadetai_ab_et = function(Aij, Bij, Uij, Vij, a, b, etai, etaj){
  return (dl2dbdxii_et(Aij, Bij, Vij, Uij, b, a, etai, etaj))
}

dl2dbdetai_ab_et = function(Aij, Bij, Uij, Vij, a, b, etai, etaj){
  return (dl2dadxii_et(Aij, Bij, Vij, Uij, b, a, etai, etaj))
}


dl2dxiidxii_et = function(Aij, Bij, Uij, Vij, a, b, xii, xij){
  n = length(Aij)
  tmp = 1+ exp(b* Vij) +  (1-xii*xij)* exp(a * Uij)
  tmp =  ifelse(tmp==0,0.001,tmp)

  res = sum( - Aij * xii^{-2}-
               Bij *  (exp(a * Uij))^2 * xij^2 /(tmp)^2)/(n-1)

  return (res)

}



dl2dxiidxij_et = function(Aij, Bij, Uij, Vij, a, b, xii, xij){
  n = length(Aij)
  tmp = 1+ exp(b* Vij) +  (1-xii*xij)* exp(a * Uij)
  tmp =  ifelse(tmp==0,0.001,tmp)

  res = sum(- Bij * exp(a * Uij)/(tmp)
            - Bij* xii * xij * (exp(a * Uij))^2 /(tmp)^2)/(n-1)

  return (res)

}

alSearch_et = function(gMat, el, gamma){
  r = length(el)
  objective.in <- rep(1,2*r)

  const.mat <- rbind(cbind(gMat,-gMat),cbind(-gMat,gMat),-diag(2*r))
  const.dir <- "<="
  const.rhs <- c(rep(1,r)*gamma+el,rep(1,r)*gamma-el,rep(0,2*r))
  res <- lpSolve::lp(direction = "min",objective.in = objective.in, const.mat = const.mat,
            const.dir = const.dir, const.rhs = const.rhs)
  res.dir <- res$solution
  a_l <- res.dir[1:r] - res.dir[-(1:r)]
  return(a_l)
}

localMLE_refine_xi_et = function(par, al, A1, B1, A2, B2, U, V, ab, xivec, etavec, i){
  xii = par
  p <- dim(A1)[1]
  a = ab[1]
  b = ab[2]

  gArray = array(0,dim = c(2*p+2, p-1))
  jix = 1
  for (j in (1:p)[-i]){
    A1ij = A1[i,j,]
    B1ij = B1[i,j,]
    A2ij = A2[i,j,]
    B2ij = B2[i,j,]
    Uij= U[i,j,]
    Vij= V[i,j,]
    xij = xivec[j]
    etai = etavec[i]
    etaj = etavec[j]

    gArray[1,jix] = dlda_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj)
    gArray[2,jix] = dldb_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj)

    gArray[2+i,jix] =  dldxii_et(A1ij, B1ij, Uij, Vij, a, b, xii, xij)
    gArray[2+j,jix] =  dldxii_et(A1ij, B1ij, Uij, Vij, a, b, xij, xii)

    gArray[2+p+i,jix] = dldxii_et(A2ij, B2ij, Vij, Uij, b, a, etai, etaj)
    gArray[2+p+j,jix] = dldxii_et(A2ij, B2ij, Vij, Uij, b, a, etaj, etai)

    jix = jix + 1

  }


  gvec = apply(gArray, c(1), mean)
  sum(al*gvec)^2

}

grr_localMLE_refine_xi_et = function(par, al, A1, B1, A2, B2, U, V, ab, xivec, etavec, i){
  xii = par
  p <- dim(A1)[1]
  a = ab[1]
  b = ab[2]

  gArray = array(0,dim = c(2*p+2, p-1))
  jix = 1
  for (j in (1:p)[-i]){
    A1ij = A1[i,j,]
    B1ij = B1[i,j,]
    A2ij = A2[i,j,]
    B2ij = B2[i,j,]
    Uij= U[i,j,]
    Vij= V[i,j,]
    xij = xivec[j]
    etai = etavec[i]
    etaj = etavec[j]

    gArray[1,jix] = dlda_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj)
    gArray[2,jix] = dldb_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj)

    gArray[2+i,jix] =  dldxii_et(A1ij, B1ij, Uij, Vij, a, b, xii, xij)
    gArray[2+j,jix] =  dldxii_et(A1ij, B1ij, Uij, Vij, a, b, xij, xii)

    gArray[2+p+i,jix] = dldxii_et(A2ij, B2ij, Vij, Uij, b, a, etai, etaj)
    gArray[2+p+j,jix] = dldxii_et(A2ij, B2ij, Vij, Uij, b, a, etaj, etai)

    jix = jix + 1

  }


  gvec = apply(gArray, c(1), mean)


  grrArray = array(0,dim = c(2*p+2, p-1))
  jix = 1
  for (j in (1:p)[-i]){
    A1ij = A1[i,j,]
    B1ij = B1[i,j,]
    A2ij = A2[i,j,]
    B2ij = B2[i,j,]
    Uij= U[i,j,]
    Vij= V[i,j,]
    xij = xivec[j]
    etai = etavec[i]
    etaj = etavec[j]

    grrArray[1,jix] = dl2dadxii_et(A1ij, B1ij, Uij, Vij, a, b, xii, xij)
    grrArray[2,jix] = dl2dbdxii_et(A1ij, B1ij, Uij, Vij, a, b, xii, xij)

    grrArray[2+i,jix] =  dl2dxiidxii_et(A1ij, B1ij, Uij, Vij, a, b, xii, xij)
    grrArray[2+j,jix] =  dl2dxiidxij_et(A1ij, B1ij, Uij, Vij, a, b, xij, xij)

    jix = jix + 1

  }

  grrvec = apply(grrArray, c(1), mean)
  2*sum(al*gvec)*sum(al*grrvec)

}

localMLE_refine_eta_et = function(par, al, A1, B1, A2, B2, U, V, ab, xivec, etavec, i){
  etai = par
  p <- dim(A1)[1]
  a = ab[1]
  b = ab[2]

  gArray = array(0,dim = c(2*p+2, p-1))
  jix = 1
  for (j in (1:p)[-i]){
    A1ij = A1[i,j,]
    B1ij = B1[i,j,]
    A2ij = A2[i,j,]
    B2ij = B2[i,j,]
    Uij= U[i,j,]
    Vij= V[i,j,]
    etaj = etavec[j]
    xii = xivec[i]
    xij = xivec[j]

    gArray[1,jix] = dlda_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj)
    gArray[2,jix] = dldb_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj)

    gArray[2+i,jix] =  dldxii_et(A1ij, B1ij, Uij, Vij, a, b, xii, xij)
    gArray[2+j,jix] =  dldxii_et(A1ij, B1ij, Uij, Vij, a, b, xij, xii)

    gArray[2+p+i,jix] = dldxii_et(A2ij, B2ij, Vij, Uij, b, a, etai, etaj)
    gArray[2+p+j,jix] = dldxii_et(A2ij, B2ij, Vij, Uij, b, a, etaj, etai)

    jix = jix + 1

  }


  gvec = apply(gArray, c(1), mean)
  sum(al*gvec)^2

}

grr_localMLE_refine_eta_et = function(par, al, A1, B1, A2, B2, U, V, ab, xivec, etavec, i){
  etai = par
  p <- dim(A1)[1]
  a = ab[1]
  b = ab[2]

  gArray = array(0,dim = c(2*p+2, p-1))
  jix = 1
  for (j in (1:p)[-i]){
    A1ij = A1[i,j,]
    B1ij = B1[i,j,]
    A2ij = A2[i,j,]
    B2ij = B2[i,j,]
    Uij= U[i,j,]
    Vij= V[i,j,]
    etaj = etavec[j]
    xii = xivec[i]
    xij = xivec[j]

    gArray[1,jix] = dlda_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj)
    gArray[2,jix] = dldb_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj)

    gArray[2+i,jix] =  dldxii_et(A1ij, B1ij, Uij, Vij, a, b, xii, xij)
    gArray[2+j,jix] =  dldxii_et(A1ij, B1ij, Uij, Vij, a, b, xij, xii)

    gArray[2+p+i,jix] = dldxii_et(A2ij, B2ij, Vij, Uij, b, a, etai, etaj)
    gArray[2+p+j,jix] = dldxii_et(A2ij, B2ij, Vij, Uij, b, a, etaj, etai)

    jix = jix + 1

  }


  gvec = apply(gArray, c(1), mean)


  grrArray = array(0,dim = c(2*p+2, p-1))
  jix = 1
  for (j in (1:p)[-i]){
    A1ij = A1[i,j,]
    B1ij = B1[i,j,]
    A2ij = A2[i,j,]
    B2ij = B2[i,j,]
    Uij= U[i,j,]
    Vij= V[i,j,]
    xij = xivec[j]
    etai = etavec[i]
    etaj = etavec[j]

    grrArray[1,jix] = dl2dadetai_ab_et(A2ij, B2ij, Uij, Vij, a, b, etai, etaj)
    grrArray[2,jix] = dl2dbdetai_ab_et(A2ij, B2ij, Uij, Vij, a, b, etai, etaj)


    grrArray[2+p+i,jix] = dl2dxiidxii_et(A2ij, B2ij, Vij, Uij, b, a, etai, etaj)
    grrArray[2+p+j,jix] = dl2dxiidxij_et(A2ij, B2ij, Vij, Uij, b, a, etai, etaj)


    jix = jix + 1

  }

  grrvec = apply(grrArray, c(1), mean)

  2*sum(al*gvec)*sum(al*grrvec)

}

globalMLE_refine_a_et = function(par, al, A1, B1, A2, B2, U, V, ab, xivec, etavec){
  p <- dim(A1)[1]
  a = par
  b = ab[2]

  gArray = array(0,dim = c(2*p+2,p))
  for (i in 1:p){
    for (j in (1:p)[-i]){
      A1ij = A1[i,j,]
      B1ij = B1[i,j,]
      A2ij = A2[i,j,]
      B2ij = B2[i,j,]
      Uij= U[i,j,]
      Vij= V[i,j,]
      xii = xivec[i]
      xij = xivec[j]
      etai = etavec[i]
      etaj = etavec[j]

      gArray[1,i] = gArray[1,i]+ dlda_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj)
      gArray[2,i] = gArray[2,i] +dldb_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj)

      gArray[2+i,i] =  gArray[2+i,i] + dldxii_et(A1ij, B1ij, Uij, Vij, a, b, xii, xij)
      gArray[2+j,i] =  gArray[2+j,i] + dldxii_et(A1ij, B1ij, Uij, Vij, a, b, xij, xii)

      gArray[2+p+i,i] = gArray[2+p+i,i] + dldxii_et(A2ij, B2ij, Vij, Uij, b, a, etai, etaj)
      gArray[2+p+j,i] = gArray[2+p+j,i] + dldxii_et(A2ij, B2ij, Vij, Uij, b, a, etaj, etai)

    }
  }

  gvec = apply(gArray, c(1), sum)/(p*(p-1))
  sum(al*gvec)^2

}

globalMLE_refine_b_et = function(par, al, A1, B1, A2, B2, U, V, ab, xivec, etavec){
  p <- dim(A1)[1]
  b = par
  a = ab[1]

  gArray = array(0,dim = c(2*p+2,p))
  for (i in 1:p){
    for (j in (1:p)[-i]){
      A1ij = A1[i,j,]
      B1ij = B1[i,j,]
      A2ij = A2[i,j,]
      B2ij = B2[i,j,]
      Uij= U[i,j,]
      Vij= V[i,j,]
      xii = xivec[i]
      xij = xivec[j]
      etai = etavec[i]
      etaj = etavec[j]

      gArray[1,i] = gArray[1,i]+ dlda_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj)
      gArray[2,i] = gArray[2,i] +dldb_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, a, b, xii, xij, etai, etaj)

      gArray[2+i,i] =  gArray[2+i,i] + dldxii_et(A1ij, B1ij, Uij, Vij, a, b, xii, xij)
      gArray[2+j,i] =  gArray[2+j,i] + dldxii_et(A1ij, B1ij, Uij, Vij, a, b, xij, xii)

      gArray[2+p+i,i] = gArray[2+p+i,i] + dldxii_et(A2ij, B2ij, Vij, Uij, b, a, etai, etaj)
      gArray[2+p+j,i] = gArray[2+p+j,i] + dldxii_et(A2ij, B2ij, Vij, Uij, b, a, etaj, etai)

    }
  }

  gvec = apply(gArray, c(1), sum)/(p*(p-1))
  sum(al*gvec)^2

}



#######-------------------------- Inference-related functions --------------
dgammada = function(Xij, Uij, Vij, a, b, xii, xij, etai, etaj){
  tmp = (1+ exp(a*Uij) +exp(b*Vij))^2
  return(
    ((1-Xij)*xii*xij*Uij*exp(a*Uij)*(1+exp(b*Vij))
     +Xij*etai*etaj*exp(b*Vij)*exp(a*Uij)*Uij)/tmp
  )
}


dgammadb = function(Xij, Uij, Vij, a, b, xii, xij, etai, etaj){
  tmp = (1+ exp(a*Uij) +exp(b*Vij))^2
  return(
    -(Xij*etai*etaj*Vij*exp(b*Vij)*(1+exp(a*Uij))
      +(1-Xij)*xii*xij*exp(b*Vij)*exp(a*Uij)*Vij)/tmp
  )
}

dgammadxii = function(Xij, Uij, Vij, a, b, xij){
  tmp = (1+ exp(a*Uij) +exp(b*Vij))
  return(
    (1-Xij)*xij*exp(a*Uij)/tmp
  )
}

dgammadetai = function(Xij, Uij, Vij, a, b, etaj){
  tmp = (1+ exp(a*Uij) +exp(b*Vij))
  return(
    -Xij*etaj*exp(b*Vij)/tmp
  )
}

gamma_root = function(Xij, Uij, Vij, a, b, xii, xij, etai, etaj){
  tmp = (1+ exp(a*Uij) +exp(b*Vij))
  alpha = xii*xij * exp(a*Uij)/tmp
  beta = etai*etaj*exp(b*Vij)/tmp
  alpha = pmin(pmax(alpha, 0.000001), 1-0.000001)
  beta  = pmin(pmax(beta, 0.000001), 1-0.000001)
  gamma = alpha +Xij*(1-alpha - beta)
  return(
    (gamma*(1-gamma))^{1/2}
  )
}


gt2 = function(Xij, Uij, Vij, a, b, xii, xij, etai, etaj){
  n = length(Xij)
  tmp = (1+ exp(a*Uij) +exp(b*Vij))
  alpha = xii*xij * exp(a*Uij)/tmp
  beta = etai*etaj*exp(b*Vij)/tmp
  gamma = alpha +Xij[-n]*(1-alpha - beta)
  return(
    (Xij[-1]-gamma)/(gamma*(1-gamma))
  )
}


### Simplified version of estNet for internal use within estTransitivity.R
estNet_trs = function(X, statsAlpha, statsBeta, globInitAlpha, globInitBeta, shrGPrm,
                  initXi = NULL, initEta = NULL, tol = 0.01, maxIter = 100){
  #Preparation:
  p = dim(X)[1]
  n = dim(X)[3]

  A1 = X[,,2:n]*(1 - X[,,2:n-1])
  B1 = (1 - X[,,2:n])*(1 - X[,,2:n-1])
  A2 = (1 - X[,,2:n])*( X[,,2:n-1])
  B2 = X[,,2:n]*(X[,,2:n-1])

  # check dimensions for initial Xi,Eta
  if((!is.null(initXi) & !(length(initXi)==p))){
    stop('Incorrect dimension for initXi')
  }
  if((!is.null(initEta) & !(length(initEta)==p))){
    stop('Incorrect dimension for initEta')
  }

  # Check the first three dimensions of statsAlpha and statsBeta
  expectedDims = c(p, p, n - 1)

  if (!all(dim(statsAlpha)[1:3] == expectedDims)) {
    stop("The first 3 dimensions of statsAlpha do not match the required dimensions of (p * p * (n - 1)).")
  }

  if (!all(dim(statsBeta)[1:3] == expectedDims)) {
    stop("The first 3 dimensions of statsBeta do not match the required dimensions of (p * p * (n - 1)).")
  }

  ### 0.1  fij function: alpha_ij  = xi_i * xi_j * f_ij
  fij = function(global, stats){
    return (exp(global[1]*stats[,,,1])/(1 + exp(global[1]*stats[,,,1]) + exp(global[2]*stats[,,,2])))
  }

  ### 0.2  gij function: beta_ij  = eta_i * eta_j * g_ij
  gij = function(global, stats){
    return ( exp(global[2]*stats[,,,2])/(1 + exp(global[1]*stats[,,,1]) + exp(global[2]*stats[,,,2])))
  }

  if (is.null(initXi)){
    xiE = rep(1,p)
  }else{
    xiE = initXi
  }

  if (is.null(initEta)){
    etaE = rep(1,p)
  }else{
    etaE = initEta
  }

  xiME = outer(xiE, xiE)
  etaME = outer(etaE, etaE)
  # Check if loglikelihood is finite
  fn = - logl(A1, B1, A2, B2, fij, globInitAlpha, statsAlpha, xiME, gij, globInitBeta, statsBeta, etaME)
  if (!is.finite(fn)) {
    stop("Please ensure initial values lead to a finite loglikelihood.")
  }


  #Initialization:
  da = length(globInitAlpha)
  db = length(globInitBeta)

  gAlphaVal.E0 = gAlphaVal.E = globInitAlpha
  gBetaVal.E0 = gBetaVal.E = globInitBeta


  if ( da < shrGPrm | db < shrGPrm){
    stop("Please check the length of global parameters.")
  }

  if (length(dim(statsAlpha))==3){
    statsAlpha = array(statsAlpha, dim = c(p,p,n-1,1))
  }
  if (length(dim(statsBeta))==3){
    statsBeta = array(statsBeta, dim = c(p,p,n-1,1))
  }



    # Updates shared global parameters simultaneously.
    # Updates unique global parameters specific to alpha simultaneously..
    # Similarly, updates unique global parameters exclusive to beta simultaneously..
    if (shrGPrm == 0){
      for (ix in 1:da){
        tmp0 = stats::optim(gAlphaVal.E0, globalMLE_group, method = "L-BFGS-B", lower = rep(0.01,da), A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateShare = F,updateAlpha = TRUE)
        gAlphaVal.E0 =  tmp0$par
      }

      fn1 = - logl_part(A1, B1, fij, gAlphaVal.E0, statsAlpha, xiME)

      for (it in 1: maxIter){
        xiMax = apply(1/fg_array(fij, gAlphaVal.E0, statsAlpha) , c(1,2), min)
        for (i in 1:(p-1)){
          for (j in (i+1):p){
            xiij = stats::optim(xiME[i,j], localMLE, method = 'L-BFGS-B', Aij = A1[i, j,], Bij = B1[i,j,], fg = fij, global = gAlphaVal.E0, stats = statsAlpha[i,j, ,,drop = FALSE], lower  = c(0.01), upper = c(xiMax[i,j]))$par
            xiME[j,i] = xiME[i,j] = xiij
          }
        }

        xiE = locEst(xiME)
        xiME = tcrossprod(xiE)

        for (ix in 1:da){
          tmp2 = stats::optim(gAlphaVal.E0, globalMLE_group, method = "L-BFGS-B", lower = rep(0.01,da), A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateShare= F,updateAlpha = TRUE)
          gAlphaVal.E = tmp2$par
        }

        fn2 = - logl_part(A1, B1, fij, gAlphaVal.E, statsAlpha, xiME)


        if(mean(abs(gAlphaVal.E -gAlphaVal.E0))<tol| fn1<=fn2) {
          break
        }else{
          gAlphaVal.E0 = gAlphaVal.E
          fn1 = fn2
        }
      }


      for (ix in 1:db){
        tmp0 = stats::optim(gBetaVal.E0, globalMLE_group, method = "L-BFGS-B", lower = rep(0.01,db), A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateShare = F,updateAlpha = F)
        gBetaVal.E0 =  tmp0$par
      }

      fn1 = - logl_part(A2, B2, gij, gBetaVal.E0, statsBeta, etaME)

      for (it in 1: maxIter){
        etaMax = apply(1/fg_array(gij, gBetaVal.E0, statsBeta) , c(1,2),min)
        for (i in 1:(p-1)){
          for (j in (i+1):p){
            etaij = stats::optim(etaME[i,j], localMLE, method = 'L-BFGS-B', Aij = A2[i,j,], Bij = B2[i,j,], fg = gij, global = gBetaVal.E0, stats = statsBeta[i,j, ,,drop = FALSE], lower  = c(0.01), upper = c(etaMax[i,j]))$par
            etaME[j,i] = etaME[i,j] = etaij
          }
        }

        etaE = locEst(etaME)
        etaME = tcrossprod(etaE)

        for (ix in 1:db){
          tmp2 = stats::optim(gBetaVal.E0, globalMLE_group, method = "L-BFGS-B", lower = rep(0.01,db),  A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateShare = F, updateAlpha = F)
          gBetaVal.E = tmp2$par
        }

        fn2 = - logl_part(A2, B2, gij, gBetaVal.E, statsBeta, etaME)


        if(mean(abs(gBetaVal.E -gBetaVal.E0))<tol| fn1<=fn2) {
          break
        }else{
          gBetaVal.E0 = gBetaVal.E
          fn1 = fn2
        }
      }


    }


    if ( shrGPrm>0){
      tmp0 = stats::optim(gAlphaVal.E0[1:shrGPrm], globalMLE_group, method = "L-BFGS-B", lower = rep(0.01,shrGPrm), A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateShare = T,updateAlpha = TRUE)
      gAlphaVal.E0[1:shrGPrm] = gBetaVal.E0[1:shrGPrm] = tmp0$par
    }

    if (da > shrGPrm){
      tmp0 = stats::optim(gAlphaVal.E0[(shrGPrm + 1):da], globalMLE_group, method = "L-BFGS-B", lower = rep(0.01, da-shrGPrm), A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateShare = F, updateAlpha = TRUE)
      gAlphaVal.E0[(shrGPrm + 1):da] =  tmp0$par
    }

    if (db > shrGPrm){
      tmp0 = stats::optim(gBetaVal.E0[(shrGPrm + 1):db], globalMLE_group, method = "L-BFGS-B", lower = rep(0.01, db-shrGPrm), A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateShare = F, updateAlpha = F)
      gBetaVal.E0[(shrGPrm + 1):db] =  tmp0$par
    }

    fn1 = - logl(A1, B1, A2, B2, fij, gAlphaVal.E0, statsAlpha, xiME, gij, gBetaVal.E0, statsBeta, etaME)


    for (it in 1: maxIter){
      xiMax = apply(1/fg_array(fij, gAlphaVal.E0, statsAlpha) , c(1,2),min)
      etaMax = apply(1/fg_array(gij, gBetaVal.E0, statsBeta) , c(1,2),min)
      for (i in 1:(p-1)){
        for (j in (i+1):p){
          xiij = stats::optim(xiME[i,j], localMLE, method = 'L-BFGS-B', Aij = A1[i, j,], Bij = B1[i,j,], fg = fij, global = gAlphaVal.E0, stats = statsAlpha[i,j, ,,drop = FALSE], lower  = c(0.01), upper = c(xiMax[i,j]))$par
          xiME[j,i] = xiME[i,j] = xiij

          etaij = stats::optim(etaME[i,j], localMLE, method = 'L-BFGS-B', Aij = A2[i,j,], Bij = B2[i,j,], fg = gij, global = gBetaVal.E0, stats = statsBeta[i,j, ,,drop = FALSE], lower  = c(0.01), upper = c(etaMax[i,j]))$par
          etaME[j,i] = etaME[i,j] = etaij
        }
      }

      xiE = locEst(xiME)
      etaE = locEst(etaME)
      xiME = tcrossprod(xiE)
      etaME = tcrossprod(etaE)

      if ( shrGPrm>0){
        tmp0 = stats::optim(gAlphaVal.E0[1:shrGPrm], globalMLE_group, method = "L-BFGS-B", lower = rep(0.01,shrGPrm), A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateShare = T,updateAlpha = TRUE)
        gAlphaVal.E[1:shrGPrm] = gBetaVal.E[1:shrGPrm] = tmp0$par
      }

      if (da > shrGPrm){
        tmp0 = stats::optim(gAlphaVal.E0[(shrGPrm + 1):da], globalMLE_group, method = "L-BFGS-B", lower = rep(0.01, da-shrGPrm), A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateShare = F, updateAlpha = TRUE)
        gAlphaVal.E[(shrGPrm + 1):da] =  tmp0$par
      }

      if (db > shrGPrm){
        tmp0 = stats::optim(gBetaVal.E0[(shrGPrm + 1):db], globalMLE_group, method = "L-BFGS-B", lower = rep(0.01, db-shrGPrm), A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateShare = F, updateAlpha = F)
        gBetaVal.E[(shrGPrm + 1):db] =  tmp0$par
      }

      fn2 = - logl(A1, B1, A2, B2, fij, gAlphaVal.E, statsAlpha, xiME, gij, gBetaVal.E, statsBeta, etaME)


      if(mean(abs(gAlphaVal.E - gAlphaVal.E0)) + mean(abs(gBetaVal.E - gBetaVal.E0))<tol| fn1<=fn2) {
        break
      }else{
        gAlphaVal.E0 = gAlphaVal.E
        gBetaVal.E0 = gBetaVal.E
        fn1 = fn2
      }
    }



  res = list()
  res$gAlphaVal = gAlphaVal.E0
  res$gBetaVal = gBetaVal.E0
  res$xi = xiE
  res$eta = etaE

  return(res)

}
