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


