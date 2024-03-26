#' Transitivity statistics for Dynamic Networks
#'
#' This function calculates statistics which summarize the higher-order structure
#' of a dynamic network object: normalized number of common and disjoint neighbours
#' between each pair of nodes in each snapshot:
#' \deqn{U^{t-1}_{i,j} = \frac{1}{p-2} \sum_k X_{i,k}^{t-1} X_{k,j}^{t-1}}
#' \deqn{V^{t-1}_{i,j} = \frac{1}{p-2} \sum_k \{ X_{i,k}^{t-1}(1 - X_{k,j}^{t-1}) +  (1 - X_{i,k}^{t-1})X_{k,j}^{t-1} \}}
#'
#' @param X A \eqn{p \times p \times n} array of network's adjacency matrices over time.
#'
#' @return A list containing the statistics:
#' \itemize{
#'   \item \code{U}: A \eqn{p \times p \times n} array of the normalized number of common neighbours
#'   \item \code{V}: A \eqn{p \times p \times n} array of the normalized number of disjoint neighbours
#' }
#'
#' @examples
#' p = 30; n = 20
#' xi = rep(0.7, p); eta = rep(0.8, p)
#' a = 30; b = 15
#'
#' # Simulate data using simulate_transitivity function
#' simulated_data = simulate_transitivity(p, n, xi, eta, a, b)
#' X = simulated_data$X
#'
#' # Calculate transitivity statistics
#' statsGlobal = statsTransitivity(X)
#'
#'
#' @export
statsTransitivity <- function(X){
  # dimensions
  p <- dim(X)[1]
  n <- dim(X)[3]
  # initialize U,V
  U <- V <- array(0,dim = c(p,p,n))
  # populate U,V
  for(t in 1:n){
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        # Common and uncommon friends calculation
        common <- sum(X[i, , t] * X[j, , t]) / (p - 2)
        uncommon <- (sum(X[i, , t] + X[j, , t]) - 2 * X[i, j, t]) / (p - 2) - 2 * common

        U[i, j, t] <- U[j, i, t] <- common
        V[i, j, t] <- V[j, i, t] <- uncommon
      }
    }
  }
  return(list(U=U,V=V))
}

#' Estimation for Autoregressive Networks with Transitivity
#'
#' This function estimates the parameters of an autoregressive network model with transitivity
#' effects. The method performs local and global likelihood estimation, then applies a
#' refinement procedure to each parameter estimate.
#'
#' This autoregressive network model contains a total of \eqn{2 + 2p} parameters.
#' The global parameters \eqn{a} and \eqn{b} are shared between the
#' formation and dissolution models, and control the influence of transitivity on the edge variables through the number of common neighbours and
#' disjoint neighbours in the previous snapshot.
#' The local parameters \eqn{\{\xi_i\}} and \eqn{\{\eta_i\}} control the node-specific propensity to form or dissolve edges, respectively.
#'
#' Edge probabilities are modeled as
#' \deqn{\alpha_{i,j}^{t-1} = \frac{\xi_i \xi_j \exp(a U_{i,j}^{t-1})}{1 + \exp(a U_{i,j}^{t-1}) + \exp(b V_{i,j}^{t-1})}}
#' \deqn{\beta_{i,j}^{t-1} = \frac{\eta_i \eta_j \exp(b V_{i,j}^{t-1})}{1 + \exp(a U_{i,j}^{t-1}) + \exp(b V_{i,j}^{t-1})}}
#' where
#' \deqn{U^{t-1}_{i,j} = \frac{1}{p-2} \sum_k X_{i,k}^{t-1} X_{k,j}^{t-1}}
#' is the normalized number of common neighbours of nodes \eqn{i} and \eqn{j} in snapshot \eqn{t-1}, and
#' \deqn{V^{t-1}_{i,j} = \frac{1}{p-2} \sum_k \{ X_{i,k}^{t-1}(1 - X_{k,j}^{t-1}) +  (1 - X_{i,k}^{t-1})X_{k,j}^{t-1} \}}
#' is the normalized number of disjoint neighbours of nodes \eqn{i} and \eqn{j} in snapshot \eqn{t-1}.
#'
#' Estimation Procedure:
#' \enumerate{
#'   \item \strong{Initialization}: automatic initialization of global parameters.
#'   \item \strong{Local Optimization}: With global parameters set to their initial values, local
#'    log-likelihood functions are maximized through a Newton method to update local parameter estimates.
#'   \item \strong{Global Optimization}: With local parameters estimated in the previous stage, a global
#'    log-likelihood function is maximized through a Newton method to
#'    estimate global parameters.
#'   \item \strong{Refinement}: For a prespecified number of iterations, global and local parameters estimates
#'   are refined by maximizing a linear restriction of the local or global likelihood functions.
#' }
#'
#'
#' @param X A \eqn{p \times p \times n} array of network's adjacency matrices over time.
#' @param U A \eqn{p \times p \times n} or \eqn{p \times p \times (n-1)} array of the normalized number of
#' common neighbour statistics. If a \eqn{p \times p \times n} array is provided, the \code{n}th slice is
#' ignored in the likelihood calculations. If not provided, it will be calculated internally from \code{X}.
#' @param V A \eqn{p \times p \times n} or \eqn{p \times p \times (n-1)} array of the normalized number of
#' disjoint neighbour statistics. If a \eqn{p \times p \times n} array is provided, the \code{n}th slice is
#' ignored in the likelihood calculations. If not provided, it will be calculated internally from \code{X}.
#' @param rSeqGlob The sequence of control parameters for the refinement of global parameters. Defaults
#' to \code{c(0.5,0.1)}. Note that \code{rSeqGlob} and \code{rSeqLoc} must have the same length.
#' @param rSeqLoc The sequence of control parameters for the refinement of local parameters. Defaults
#' to \code{c(0.5,0.1)}.
#' @param verbose An indicator, if \code{TRUE} there will be printed console output providing optimization
#' progress. defaults to \code{FALSE}.
#'
#' @return A list containing the estimated parameters:
#' \itemize{
#'   \item \code{gVal}: Estimated global parameters \eqn{a, b}.
#'   \item \code{xi}: Estimated values of \eqn{\xi_1, \dots, \xi_p}.
#'   \item \code{eta}: Estimated values of \eqn{\eta_1, \dots, \eta_p}.
#' }
#'
#' @examples
#' p = 30; n = 20
#' xi = rep(0.7, p); eta = rep(0.8, p)
#' a = 30; b = 15
#'
#' # Simulate data using simulate_transitivity function
#' simulated_data = simulate_transitivity(p, n, xi, eta, a, b)
#' X = simulated_data$X
#' U = simulated_data$U
#' V = simulated_data$V
#'
#' result = estTransitivity(X, U, V)
#'
#' @export
estTransitivity <- function(X,U=NULL,V=NULL,
                            rSeqGlob=c(0.5,0.1),rSeqLoc=c(0.5,0.1),
                            verbose=FALSE){
  # dimensions
  p <- dim(X)[1]
  n <- dim(X)[3]

  # additional parameter checking
  if(!(length(rSeqGlob)==length(rSeqLoc))){
    stop('Global and local refinement control sequences must be the same length')
  }

  # transitivity statistics
  statsGlobal <- NULL
  # U: common neighbours
  if(is.null(U)){
    statsGlobal <- statsTransitivity(X)
    # likelihood uses clipped U without the final snapshot
    Uc <- statsGlobal$U[,,-n]
  }
  else{
    if(all(dim(U)==c(p,p,n))){
      Uc <- U[,,-n]
    }
    else if(all(dim(U)==c(p,p,n))){
      Uc <- U
    }
    else{
      stop('Incorrect dimensions for U')
    }
  }
  # V: uncommon neighbours
  if(is.null(V)){
    if(is.null(statsGlobal)){
      statsGlobal <- statsTransitivity(X)
    }
    # likelihood uses clipped U without the final snapshot
    Vc <- statsGlobal$V[,,-n]
  }
  else{
    if(all(dim(V)==c(p,p,n))){
      Vc <- V[,,-n]
    }
    else if(all(dim(V)==c(p,p,n))){
      Vc <- V
    }
    else{
      stop('Incorrect dimensions for V')
    }
  }

  # dynamic network counts
  A1 = X[,,2:n]*(1 - X[,,2:n-1])
  B1 = (1 - X[,,2:n])*(1 - X[,,2:n-1])
  A2 = (1 - X[,,2:n])*( X[,,2:n-1])
  B2 = X[,,2:n]*(X[,,2:n-1])

  # additional control for optimization of a,b: want exp(bV) and exp(aU) < Inf
  ab_max <- (400*(p-1)) / max(c(Uc,Vc))
  #ximax_init <- 1
  #etamax_init <- 1
  maxIter <- 5
  tol <- 1e-4

  #Rough Initialization:
  # ab1 <- globInit
  # xiME = pmax(pmin(apply((1 +exp(ab1[1]*Uc) + exp(ab1[2]*Vc) )/exp(ab1[1]*Uc),c(1,2),min),ximax_init)-0.2,0)
  # etaME =  pmax(pmin(apply((1 +exp(ab1[2]*Vc) + exp(ab1[1]*Uc) )/exp(ab1[2]*Vc),c(1,2),min),etamax_init)-0.2,0)
  # xiE = thetaEst_et(xiME)
  # etaE =  thetaEst_et(etaME)
  # if(verbose){
  #   cat('Current xi quartiles:',stats::quantile(xiE),
  #       '\nCurrent eta quartiles:',stats::quantile(etaE),'\n')
  # }

  tmp0 = stats::optim(c(tol,tol), globalMLE_ab_et, gr = grr_globalMLE_ab_et, method = "L-BFGS-B",
                      lower = c(0,0), upper=rep(ab_max,2),
                      A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, thetavec = rep(1,p), etavec = rep(1,p))
  ab1 = tmp0$par
  fn1 = tmp0$value
  # if(verbose){
  #   cat('Current (a,b):',ab1,'\nDone rough initialization\n')
  # }

  for (it in 1:maxIter){
    ximax = apply((1 +exp(ab1[1]*Uc) + exp(ab1[2]*Vc) )/exp(ab1[1]*Uc),c(1,2),min)
    etamax =  apply((1 +exp(ab1[2]*Vc) + exp(ab1[1]*Uc) )/exp(ab1[2]*Vc),c(1,2),min)
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        xiij = stats::optim(xiME[i,j], localMLE_init_et, method = 'L-BFGS-B', Aij = A1[i,j,], Bij = B1[i,j,], Uij=Uc[i,j,], Vij=Vc[i,j,], ab = ab1, lower  = c(0), upper = c(ximax[i,j]))$par
        xiME[j,i] = xiME[i,j] = xiij

        etaij = stats::optim(etaME[i,j], localMLE_init_et, method = 'L-BFGS-B', Aij = A2[i,j,], Bij = B2[i,j,], Uij=Vc[i,j,], Vij=Uc[i,j,], ab = c(ab1[2],ab1[1]), lower  = c(0), upper = c(etamax[i,j]))$par
        etaME[j,i] = etaME[i,j] = etaij

      }
    }

    xiE = thetaEst_et(xiME)
    etaE = thetaEst_et(etaME)
    # if(verbose){
    #   cat('Current xi quartiles:',stats::quantile(xiE),
    #       '\nCurrent eta quartiles:',stats::quantile(etaE),'\n')
    # }
    xiME = outer(xiE, xiE)
    etaME = outer(etaE, etaE)


    tmp2 = stats::optim(ab1, globalMLE_ab_et, gr = grr_globalMLE_ab_et, method = "L-BFGS-B",
                        lower = c(0,0), upper=rep(ab_max,2),
                        A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, thetavec = xiE, etavec = etaE)
    ab2 = tmp2$par
    fn2 = tmp2$value
    # if(verbose){
    #   cat('Current (a,b):',ab1,'\nDone initialization iter',it,'\n')
    # }

    if(mean(abs(ab2-ab1))<tol | fn1<=fn2) {
      break
    }else{
      ab1 = ab2
      fn1 = fn2
    }
  }
  if(verbose){
    cat('Done initialization\n')
  }

  ##### Initial value: ab1, xiE and etaE
  ##### Now we apply ONE ITERATION OF the global MLE for (a,b) and local MLE for xii and etai

  # additional control for estimation of a,b
  ab_max_est <- min(10*max(ab1),ab_max)
  ximax = apply((1 +exp(ab1[1]*Uc) + exp(ab1[2]*Vc) )/exp(ab1[1]*Uc),c(1,2),min)
  etamax =  apply((1 +exp(ab1[2]*Vc) + exp(ab1[1]*Uc) )/exp(ab1[2]*Vc),c(1,2),min)

  # local MLE
  for (i in 1:p){
    xiE[i] = stats::optim(xiE[i],  localMLE_et,  method = 'L-BFGS-B', Ai = A1[i,-i,] , Bi = B1[i,-i,], Ui = Uc[i,-i,], Vi=Vc[i,-i,], ab=ab1, thetavec_ic= xiE[-i], lower  = c(tol), upper = min(ximax[i,-i]/xiE[-i]))$par
    etaE[i] = stats::optim(etaE[i],  localMLE_et,  method = 'L-BFGS-B', Ai = A2[i,-i,] , Bi = B2[i,-i,], Ui = Vc[i,-i,], Vi=Uc[i,-i,], ab=c(ab1[2],ab1[1]), thetavec_ic= etaE[-i], lower  = c(tol), upper = min(etamax[i,-i]/etaE[-i]))$par
  }
  if(verbose){
    cat('Current xi quartiles:',stats::quantile(xiE),
        '\nCurrent eta quartiles:',stats::quantile(etaE),'\nDone local optimization\n')
  }

  # global MLE
  tmp2 = stats::optim(ab1, globalMLE_ab_et, gr = grr_globalMLE_ab_et, method = "L-BFGS-B",
                      lower = c(tol,tol), upper=rep(ab_max_est,2),
                      A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, thetavec = xiE, etavec = etaE)
  ab1 = tmp2$par
  fn1 = tmp2$value
  if(verbose){
    cat('Current (a,b):',ab1,'\nDone global optimization\n')
  }

  # Linear refinement
  # additional rate parameters
  G = 2
  Gc = 2*p
  q = G + Gc
  S_g = p*(p-1)
  S_gc = p-1
  c_g2 = 0.01*max(q *log(n*S_g)/sqrt(n*S_g),q^{3/2} *(log(n*S_g))^{3/2}/sqrt(n)/S_g)
  c_gc2 = 0.01*max(q *log(n * S_gc)/sqrt(n * S_gc),q^{3/2} *(log(n*S_gc))^{3/2}/sqrt(n)/S_gc)
  delta_n_sqrt = sqrt(max(c_g2, c_gc2))
  gamma_xi = 0.5 * delta_n_sqrt
  gamma_a = 0.01 * delta_n_sqrt
  gamma_b = 0.01 * delta_n_sqrt
  gamma_eta = 0.5 * delta_n_sqrt
  # initialize
  ab2 = ab1
  xiR = xiE
  etaR = etaE
  # loop refinement
  for (rr in 1:length(rSeqGlob)){
    # set refinement control parameters
    rtildeGlob <- (p-2)*rSeqGlob[rr]
    rtildeLoc <- rSeqLoc[rr]
    # initialize gradient array
    gArray = array(0, dim = c(2+p*2,2+p*2,p ))
    # populate gradient array
    for (i in 1:p){
      for (j in (1:p)[-i]){
        A1ij = A1[i,j,]
        B1ij = B1[i,j,]
        A2ij = A2[i,j,]
        B2ij = B2[i,j,]
        Uij= Uc[i,j,]
        Vij= Vc[i,j,]
        xii = xiE[i]
        xij = xiE[j]
        etai = etaE[i]
        etaj = etaE[j]
        gArray[1,1,i] = gArray[1,1,i] + dl2dada_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], xii, xij, etai, etaj)
        gArray[1,2,i] = gArray[2,1,i] = gArray[1,2,i] + dl2dadb_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], xii, xij, etai, etaj)
        gArray[2,2,i] = gArray[2,2,i] + dl2dbdb_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], xii, xij, etai, etaj)

        gArray[1,2+i,i] =  gArray[2+i,1,i] = gArray[1,2+i,i] + dl2dadthetai_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xii, xij)
        gArray[1,2+j,i] =  gArray[2+j,1,i] = gArray[1,2+j,i] + dl2dadthetai_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xij, xii)
        gArray[2,2+i,i] =  gArray[2+i,2,i] = gArray[2,2+i,i] + dl2dbdthetai_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xii, xij)
        gArray[2,2+j,i] =  gArray[2+j,2,i] = gArray[2,2+j,i] + dl2dbdthetai_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xij, xii)

        gArray[1,2+p+i,i] = gArray[2+p+i,1,i] = gArray[1,2+p+i,i] + dl2dadetai_ab_et(A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], etai, etaj)
        gArray[1,2+p+j,i] = gArray[2+p+j,1,i] = gArray[1,2+p+j,i] + dl2dadetai_ab_et(A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], etaj, etai)
        gArray[2,2+p+i,i] = gArray[2+p+i,2,i] = gArray[2,2+p+i,i] + dl2dbdetai_ab_et(A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], etai, etaj)
        gArray[2,2+p+j,i] = gArray[2+p+j,2,i] = gArray[2,2+p+j,i] + dl2dbdetai_ab_et(A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], etaj, etai)

        gArray[2+i,2+i,i] = gArray[2+i,2+i,i] + dl2dthetaidthetai_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xii, xij)
        gArray[2+j,2+j,i] = gArray[2+j,2+j,i] + dl2dthetaidthetai_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xij, xij)
        gArray[2+i,2+j,i] = gArray[2+j,2+i,i] =  gArray[2+i,2+j,i] + dl2dthetaidthetaj_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xii, xij)

        gArray[2+p+i,2+p+i,i] = gArray[2+p+i,2+p+i,i] + dl2dthetaidthetai_et(A2ij, B2ij, Vij, Uij, ab1[2], ab1[1], etai, etaj)
        gArray[2+p+j,2+p+j,i] =  gArray[2+p+j,2+p+j,i] +dl2dthetaidthetai_et(A2ij, B2ij, Vij, Uij, ab1[2], ab1[1], etaj, etai)
        gArray[2+p+i,2+p+j,i] = gArray[2+p+j,2+p+i,i] = gArray[2+p+i,2+p+j,i] +dl2dthetaidthetaj_et(A2ij, B2ij, Vij, Uij, ab1[2], ab1[1], etai, etaj)

      }
    }
    # split gradient array
    g_local = gArray/(p-1)
    g_global = apply(gArray, c(1,2),sum)/(p*(p-1))
    g_all = array(c(g_global, g_global, g_local, g_local),dim = c(2+2*p,2+2*p,2*p+2))
    # space for linear projections
    El = diag(1, 2+2*p)
    Al = array(0, dim = c(2+2*p,2+2*p)) #order: a,b, xi1, \dots, xip, eta1, \dots, etap
    # search for linear projection
    Al[,1] = alSearch_et(g_all[,,1], El[,1], gamma_a)
    Al[,2] = alSearch_et(g_all[,,2], El[,2], gamma_b)
    for (i in 3:(2 + p)){
      Al[,i] = alSearch_et(g_all[,,i], El[,i], gamma_xi)
    }
    for (i in (p+3):(2 + 2*p)){
      Al[,i] = alSearch_et(g_all[,,i], El[,i], gamma_eta)
    }
    # refinement for a,b
    ab2[1] = stats::optim(ab1[1],  globalMLE_refine_a_et,  method = 'L-BFGS-B', al = Al[,1], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, ab = ab1, thetavec = xiE, etavec = etaE, lower = max(0.01, ab1[1]-rtildeGlob), upper = ab1[1]+rtildeGlob)$par
    ab2[2] = stats::optim(ab1[2],  globalMLE_refine_b_et,  method = 'L-BFGS-B', al = Al[,2], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, ab = ab1, thetavec = xiE, etavec = etaE, lower = max(0.01, ab1[2]-rtildeGlob), upper = ab1[2]+rtildeGlob)$par
    if(verbose){
      cat('Current (a,b):',ab2,'\nDone global parameter refinement',rr,'of',length(rSeqGlob),'\n')
    }
    # refinement for xi,eta
    ximax = apply((1 +exp(ab1[1]*Uc) + exp(ab1[2]*Vc) )/exp(ab1[1]*Uc),c(1,2),min)
    etamax =  apply((1 +exp(ab1[2]*Vc) + exp(ab1[1]*Uc) )/exp(ab1[2]*Vc),c(1,2),min)
    for (i in 1:p){
      xiR[i] = stats::optim(xiE[i],  localMLE_refine_theta_et, grr_localMLE_refine_theta_et,  method = 'L-BFGS-B', al = Al[,2+i], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, ab = ab1, thetavec = xiE, etavec = etaE, i = i, lower = max(0.01, xiE[i]-rtildeLoc), upper = min(c(ximax[i,-i]/xiE[-i], xiE[i]+rtildeLoc)))$par
      etaR[i] = stats::optim(etaE[i],  localMLE_refine_eta_et, grr_localMLE_refine_eta_et,  method = 'L-BFGS-B', al = Al[,2+p+i], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, ab = ab1, thetavec = xiE, etavec = etaE, i = i, lower = max(0.01, etaE[i]-rtildeLoc), upper = min(c(etamax[i,-i]/etaE[-i], etaE[i]+rtildeLoc)))$par
    }
    if(verbose){
      cat('Current xi quartiles:',stats::quantile(xiE),
          '\nCurrent eta quartiles:',stats::quantile(etaE),
          '\nDone local parameter refinement',rr,'of',length(rSeqGlob),'\n')
    }
    # update
    ab1 = ab2
    xiE = xiR
    etaE = etaR
  }
  # store estimators
  res <- list()
  #res$gAlphaVal <- ab1
  #res$gBetaVal <- ab1
  res$gVal <- ab1
  res$xi <- xiE
  res$eta <- etaE

  return(res)
}
