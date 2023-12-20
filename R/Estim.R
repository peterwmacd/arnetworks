# wrapper function for estimation with the transitivity model

estim_transitivity <- function(X,
                               # optimization parameters
                               iter=100,delta=1e-4,verbose=FALSE,
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
  ab_max <- 400 / max(c(U,V))
  #Initialization:
  ab1 <- ab_init
  thetaME = pmax(pmin(apply((1 +exp(ab1[1]*U) + exp(ab1[2]*V) )/exp(ab1[1]*U),c(1,2),min),thetamax_init)-0.2,0)
  etaME =  pmax(pmin(apply((1 +exp(ab1[2]*V) + exp(ab1[1]*U) )/exp(ab1[2]*V),c(1,2),min),etamax_init)-0.2,0)
  thetaE = thetaEst(thetaME)
  etaE =  thetaEst(etaME)
  if(verbose){
    print(quantile(thetaE))
    print(quantile(etaE))
  }

  tmp0 = stats::optim(ab1, globalMLE_ab, gr = grr_globalMLE_ab, method = "L-BFGS-B",
                      lower = c(0,0), upper=rep(ab_max,2),
                      A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = U, V = V, thetavec = thetaE, etavec = etaE)
  ab1 = tmp0$par
  fn1 = tmp0$value
  if(verbose){
    print(ab1)
  }
  print('done rough initialization')

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
      print(quantile(thetaE))
      print(quantile(etaE))
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
      print(quantile(thetaE))
      print(quantile(etaE))
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
  # store estimators
  res <- list()
  res$theta <- thetaE
  res$eta <- etaE
  res$a <- ab1[1]
  res$b <- ab1[2]
  return(res)
}
