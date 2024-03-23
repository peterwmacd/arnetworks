# wrapper function for estimation with the transitivity model


estim_transitivity <- function(X,
                               # optimization parameters
                               iter=10,delta=1e-4,verbose=FALSE,
                               # initializer tuning (w/ default values)
                               ab_init=c(0,0),
                               thetamax_init=1,etamax_init=1){
  # dimensions
  p <- dim(X)[1]
  n <- dim(X)[3]
  # transitivity statistics
  tstat <- transitivity_stats(X)
  U <- tstat$U
  V <- tstat$V
  # dynamic network counts
  A1 = X[,,2:n]*(1 - X[,,2:n-1])
  B1 = (1 - X[,,2:n])*(1 - X[,,2:n-1])
  A2 = (1 - X[,,2:n])*( X[,,2:n-1])
  B2 = X[,,2:n]*(X[,,2:n-1])
  # rough maximum for a,b: want exp(bV) and exp(aU) < Inf; ok if < 400
  ab_max <- (400*(p-1)) / max(c(U,V))
  #Initialization:
  ab1 <- ab_init
  thetaME = pmax(pmin(apply((1 +exp(ab1[1]*U) + exp(ab1[2]*V) )/exp(ab1[1]*U),c(1,2),min),thetamax_init)-0.2,0)
  etaME =  pmax(pmin(apply((1 +exp(ab1[2]*V) + exp(ab1[1]*U) )/exp(ab1[2]*V),c(1,2),min),etamax_init)-0.2,0)
  thetaE = thetaEst(thetaME)
  etaE =  thetaEst(etaME)
  if(verbose){
    print(stats::quantile(thetaE))
    print(stats::quantile(etaE))
  }

  tmp0 = stats::optim(ab1, globalMLE_ab, gr = grr_globalMLE_ab, method = "L-BFGS-B",
                      lower = c(0,0), upper=rep(ab_max,2),
                      A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, thetavec = thetaE, etavec = etaE)
  ab1 = tmp0$par
  fn1 = tmp0$value
  if(verbose){
    print(ab1)
    print('done rough initialization')
  }

  for (it in 1:iter){
    thetamax = apply((1 +exp(ab1[1]*U) + exp(ab1[2]*V) )/exp(ab1[1]*U),c(1,2),min)
    etamax =  apply((1 +exp(ab1[2]*V) + exp(ab1[1]*U) )/exp(ab1[2]*V),c(1,2),min)
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        thetaij = stats::optim(thetaME[i,j], localMLE.init, method = 'L-BFGS-B', Aij = A1[i,j,], Bij = B1[i,j,], Uij=U[i,j,], Vij=V[i,j,], ab = ab1, lower  = c(0), upper = c(thetamax[i,j]))$par
        thetaME[j,i] = thetaME[i,j] = thetaij

        etaij = stats::optim(etaME[i,j], localMLE.init, method = 'L-BFGS-B', Aij = A2[i,j,], Bij = B2[i,j,], Uij=V[i,j,], Vij=U[i,j,], ab = c(ab1[2],ab1[1]), lower  = c(0), upper = c(etamax[i,j]))$par
        etaME[j,i] = etaME[i,j] = etaij

      }
    }

    thetaE = thetaEst(thetaME)
    etaE = thetaEst(etaME)
    if(verbose){
      print(stats::quantile(thetaE))
      print(stats::quantile(etaE))
    }
    thetaME = outer(thetaE, thetaE)
    etaME = outer(etaE, etaE)


    tmp2 = stats::optim(ab1, globalMLE_ab, gr = grr_globalMLE_ab, method = "L-BFGS-B",
                        lower = c(0,0), upper=rep(ab_max,2),
                        A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, thetavec = thetaE, etavec = etaE)
    ab2 = tmp2$par
    fn2 = tmp2$value
    if(verbose){
      print(ab2)
    }

    if(mean(abs(ab2 -ab1))<delta | fn1<=fn2) {
      break
    }else{
      ab1 = ab2
      fn1 = fn2
      if(verbose){
        print(paste0('done initialization iter ',it))
      }
    }
  }
  if(verbose){
    print('done initialization')
  }

  ##### Initial value: ab1, thetaE and etaE
  ##### Now we apply the global MLE for (a,b) and local MLE for thetai and etai

  # shrink ab_max if ab1 << ab_max
  ab_max_est <- min(10*max(ab1),ab_max)
  for (it in 1: iter){
    thetamax = apply((1 +exp(ab1[1]*U) + exp(ab1[2]*V) )/exp(ab1[1]*U),c(1,2),min)
    etamax =  apply((1 +exp(ab1[2]*V) + exp(ab1[1]*U) )/exp(ab1[2]*V),c(1,2),min)

    for (i in 1:p){
      thetaE[i] = stats::optim(thetaE[i],  localMLE,  method = 'L-BFGS-B', Ai = A1[i,-i,] , Bi = B1[i,-i,], Ui = U[i,-i,], Vi=V[i,-i,], ab=ab1, thetavec_ic= thetaE[-i], lower  = c(delta), upper = min(thetamax[i,-i]/thetaE[-i]))$par
      etaE[i] = stats::optim(etaE[i],  localMLE,  method = 'L-BFGS-B', Ai = A2[i,-i,] , Bi = B2[i,-i,], Ui = V[i,-i,], Vi=U[i,-i,], ab=c(ab1[2],ab1[1]), thetavec_ic= etaE[-i], lower  = c(delta), upper = min(etamax[i,-i]/etaE[-i]))$par
    }
    if(verbose){
      print(stats::quantile(thetaE))
      print(stats::quantile(etaE))
    }

    tmp2 = stats::optim(ab1, globalMLE_ab, gr = grr_globalMLE_ab, method = "L-BFGS-B",
                        lower = c(delta,delta), upper=rep(ab_max_est,2),
                        A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, thetavec = thetaE, etavec = etaE)
    ab2 = tmp2$par
    fn2 = tmp2$value
    if(verbose){
      print(ab2)
    }

    if(mean(abs(ab2 -ab1))<delta| fn1<=fn2) {
      break
    }else{
      ab1 = ab2
      fn1 = fn2
      if(verbose){
        print(paste0('done local estimation iter',it))
      }
    }
  }
  if(verbose){
    print('done local estimation')
  }

  # refinement

  # tuning parameters
  r_seq <- c(0.5,0.1)
  # rate parameters
  G = 2
  Gc = 2*p
  q = G + Gc
  S_g = p*(p-1)
  S_gc = p-1
  c_g2 = 0.01*max(q *log(n*S_g)/sqrt(n*S_g),q^{3/2} *(log(n*S_g))^{3/2}/sqrt(n)/S_g)
  c_gc2 = 0.01*max(q *log(n * S_gc)/sqrt(n * S_gc),q^{3/2} *(log(n*S_gc))^{3/2}/sqrt(n)/S_gc)
  delta_n_sqrt = sqrt(max(c_g2, c_gc2))
  gamma_theta = 0.5 * delta_n_sqrt
  gamma_a = 0.01 * delta_n_sqrt
  gamma_b = 0.01 * delta_n_sqrt
  gamma_eta = 0.5 * delta_n_sqrt
  # initialize
  ab2 = ab1
  thetaR = thetaE
  etaR = etaE
  # loop refinement
  refine <- 1
  for (r_tilda in r_seq){
    # initialize gradient array
    gArray = array(0, dim = c(2+p*2,2+p*2,p ))
    # populate gradient array
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
    # split gradient array
    g_local = gArray/(p-1)
    g_global = apply(gArray, c(1,2),sum)/(p*(p-1))
    g_all = array(c(g_global, g_global, g_local, g_local),dim = c(2+2*p,2+2*p,2*p+2))
    # space for linear projections
    El = diag(1, 2+2*p)
    Al = array(0, dim = c(2+2*p,2+2*p)) #order: a,b, theta1, \dots, thetap, eta1, \dots, etap
    # search for linear projection
    Al[,1] = alSearch(g_all[,,1], El[,1], gamma_a)
    Al[,2] = alSearch(g_all[,,2], El[,2], gamma_b)
    for (i in 3:(2 + p)){
      Al[,i] = alSearch(g_all[,,i], El[,i], gamma_theta)
    }
    for (i in (p+3):(2 + 2*p)){
      Al[,i] = alSearch(g_all[,,i], El[,i], gamma_eta)
    }
    # refinement for a,b
    ab2[1] = stats::optim(ab1[1],  globalMLE_refine_a,  method = 'L-BFGS-B', al = Al[,1], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, ab = ab1, thetavec = thetaE, etavec = etaE, lower = max(0.01, ab1[1]-(p-1)*r_tilda), upper = ab1[1]+(p-1)*r_tilda)$par
    ab2[2] = stats::optim(ab1[2],  globalMLE_refine_b,  method = 'L-BFGS-B', al = Al[,2], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, ab = ab1, thetavec = thetaE, etavec = etaE, lower = max(0.01, ab1[2]-(p-1)*r_tilda), upper = ab1[2]+(p-1)*r_tilda)$par
    if(verbose){
      print(ab2)
      print(paste0('Done a,b refinement ',refine))
    }
    # refinement for theta,eta
    thetamax = apply((1 +exp(ab1[1]*U) + exp(ab1[2]*V) )/exp(ab1[1]*U),c(1,2),min)
    etamax =  apply((1 +exp(ab1[2]*V) + exp(ab1[1]*U) )/exp(ab1[2]*V),c(1,2),min)
    for (i in 1:p){
      thetaR[i] = stats::optim(thetaE[i],  localMLE_refine_theta, grr_localMLE_refine_theta,  method = 'L-BFGS-B', al = Al[,2+i], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, ab = ab1, thetavec = thetaE, etavec = etaE, i = i, lower = max(0.01, thetaE[i]-r_tilda), upper = min(c(thetamax[i,-i]/thetaE[-i], thetaE[i]+r_tilda)))$par
      etaR[i] = stats::optim(etaE[i],  localMLE_refine_eta, grr_localMLE_refine_eta,  method = 'L-BFGS-B', al = Al[,2+p+i], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, ab = ab1, thetavec = thetaE, etavec = etaE, i = i, lower = max(0.01, etaE[i]-r_tilda), upper = min(c(etamax[i,-i]/etaE[-i], etaE[i]+r_tilda)))$par
    }
    if(verbose){
      print(stats::quantile(thetaR))
      print(stats::quantile(etaR))
      print(paste0('Done theta,eta refinement ',refine))
    }
    # update
    ab1 = ab2
    thetaE = thetaR
    etaE = etaR
    refine <- refine+1
  }
  # store estimators
  res <- list()
  res$theta <- thetaE
  res$eta <- etaE
  res$a <- ab1[1]
  res$b <- ab1[2]
  return(res)
}
