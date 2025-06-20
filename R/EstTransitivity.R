#' Transitivity Statistics for Dynamic Networks
#'
#' This function calculates statistics which summarize the higher-order structure
#' of a dynamic network object: normalized number of common and disjoint neighbours
#' between each pair of nodes in each snapshot:
#' \deqn{U_{i,j}^{t} = \frac{1}{p-2} \sum_{k \neq i,j} X_{i,k}^{t}X_{j,k}^{t}}
#' \deqn{V_{i,j}^{t} =  \frac{1}{p-2} \sum_{k \neq i,j} \{ X_{i,k}^{t-1}(1-X_{j,k}^{t}) + (1-X_{i,k}^{t-1})X_{j,k}^{t}\}}
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
#' # Simulate data using simulateTransitivity function
#' simulated_data = simulateTransitivity(p, n, xi, eta, a, b)
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

#' Estimation/Inference for Autoregressive Networks with Transitivity
#'
#' This function estimates the parameters of an autoregressive network model with transitivity
#' effects, as detailed in "Autoregressive Networks with Dependent Edges". The method performs local and global likelihood estimation,  followed by a
#' projection-based refinement step for each parameter estimate.
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
#' \deqn{U_{i,j}^{t} = \frac{1}{p-2} \sum_{k \neq i,j} X_{i,k}^{t}X_{j,k}^{t}}
#' is the normalized number of common neighbours of nodes \eqn{i} and \eqn{j} in snapshot \eqn{t-1}, and
#' \deqn{V_{i,j}^{t} =  \frac{1}{p-2} \sum_{k \neq i,j} \{ X_{i,k}^{t-1}(1-X_{j,k}^{t}) + (1-X_{i,k}^{t-1})X_{j,k}^{t}\}}
#' is the normalized number of disjoint neighbours of nodes \eqn{i} and \eqn{j} in snapshot \eqn{t-1}.
#'
#' Estimation Procedure:
#' \enumerate{
#'   \item \strong{Initialization}: Set initial values for local parameters.
#'   \item \strong{Global Optimization}: Utilize initialized local parameters within a global
#'    log-likelihood function, optimizing it through a quasi-Newton method to
#'    estimate global parameters.
#'   \item \strong{Local Optimization}: With global parameters estimated in the previous stage, local
#'    log-likelihood functions are maximized through a Newton method to update local parameter estimates.
#'    \item \strong{IMoM Estimation}: Apply an efficient Iterative Method-of-Moments (IMoM) estimation to reduce the sensitivity to initial values.
#'   \item \strong{Refinement}: Incorporate a projection-based refinement process to improve the
#'    previous maximum likelihood estimates (MLEs). The tuning parameters \eqn{\tau} are selected internally as documented in Chang et al. (2024).
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
#' @param initXi Initial values for \eqn{\xi_1, \dots, \xi_p}. It will be calculated internally if NULL.
#' @param initEta Initial values for \eqn{\eta_1, \dots, \eta_p}. It will be calculated internally if NULL.
#' @param tauSeq_a A numeric vector specifying candidate values for the tuning parameter \eqn{\tau_a}, used in the projection-based refinement of the global parameter \eqn{a}.
#' The selected value is multiplied by \eqn{\Delta_n^{1/2}} to form the constraint threshold in the projection step (see Equation 5.10). If a vector of length one is provided, the corresponding value is used directly without selection.
#' See also:
#' \deqn{
#' \hat{\varphi}_l = \arg \min_{\mathbf{u} \in \mathbb{R}^q} \|\mathbf{u}\|_1 \quad \text{s.t.} \quad \left\| \left( \frac{1}{n - m} \sum_{t=m+1}^{n} \nabla_\theta g_t^{(l)}(\tilde{\boldsymbol{\theta}}) \right)^\top \mathbf{u} - \mathbf{e}_l \right\|_\infty \leq \tau \Delta_n^{1/2},
#' }
#' where \eqn{\tau} is chosen from the \code{tauSeq_*} vector.
#' @param tauSeq_b A numeric vector specifying candidate values for the tuning parameter \eqn{\tau_b}, used in the projection-based refinement of the global parameter \eqn{b}.
#' @param tauSeq_xi A numeric vector specifying candidate values for \eqn{\tau_{\xi}}, used to refine the local formation parameters \eqn{\xi_1, \dots, \xi_p}. Selection is performed individually for each node \eqn{i \in \{1, \dots, p\}}.
#' @param tauSeq_eta A numeric vector specifying candidate values for \eqn{\tau_{\eta}}, used to refine the local dissolution parameters \eqn{\eta_1, \dots, \eta_p}. Selection is performed individually for each node \eqn{i \in \{1, \dots, p\}}.
#' @param rSeqGlob A numeric vector of length two specifying the local search radii used in the two-stage refinement of the global parameters \eqn{a} and \eqn{b}. Each value defines the size of the neighborhood around the initial estimate within which the final refinement is performed. Defaults to \code{c(10, 2)}.
#' @param rSeqLoc A numeric vector of length two specifying the local search radii used in the two-stage refinement of the local parameters \eqn{\xi_1, \dots, \xi_p} and \eqn{\eta_1, \dots, \eta_p}. Defaults to \code{c(0.2, 0.05)}.
#' @param verbose Logical. If \code{TRUE}, there will be console output providing optimization
#' progress. Defaults to \code{FALSE}.
#' @param doInference Logical. If \code{TRUE}, the function performs inference by computing the asymptotic standard errors for all estimated parameters. Defaults to \code{TRUE}.
#'
#' @return A list containing the estimated parameters and optional inference results:
#' \itemize{
#'   \item \code{gVal_init}: Initial global parameter estimates \eqn{(a, b)}.
#'   \item \code{xi_init}: Initial local formation parameter estimates \eqn{(\xi_1, \dots, \xi_p)}.
#'   \item \code{eta_init}: Initial local dissolution parameter estimates \eqn{(\eta_1, \dots, \eta_p)}.
#'   \item \code{gVal}: Refined global parameter estimates \eqn{(a, b)}.
#'   \item \code{xi}: Refined local formation parameter estimates \eqn{(\xi_1, \dots, \xi_p)}.
#'   \item \code{eta}: Refined local dissolution parameter estimates \eqn{(\eta_1, \dots, \eta_p)}.
#'   \item \code{se_estimates}: A numeric vector of length \eqn{2 + 2p} containing the asymptotic standard errors for the final parameter estimates \eqn{\hat{\theta}}. The entries correspond to:
#'   \eqn{\hat{a}, \hat{b}, \hat{\xi}_1, \dots, \hat{\xi}_p, \hat{\eta}_1, \dots, \hat{\eta}_p}.
#'   Returned only if \code{doInference = TRUE}.
#' }

#'
#' @examples
#' p = 20; n = 30
#' xi = rep(0.8, p); eta = rep(0.9, p)
#' a = 25; b = 15
#' initXi = runif(p,0,1); initEta =runif(p,0,1)
#'
#' # Simulate data using simulateTransitivity function
#' simulated_data = simulateTransitivity(p, n, xi, eta, a, b)
#' X = simulated_data$X
#' U = simulated_data$U
#' V = simulated_data$V
#'
#' result = estTransitivity(X, U, V, initXi, initEta, tauSeq_a = seq(0.3, 0.4, 2), tauSeq_b = seq(0.3, 0.4, 2), tauSeq_xi = seq(0.05, 0.06, 2), tauSeq_eta = seq(0.05, 0.06, 2))
#'@references
#' Chang, J., Fang, Q., Kolaczyk, E. D., MacDonald, P. W., & Yao, Q. (2024). Autoregressive Networks with Dependent Edges.
#' @export
estTransitivity <- function(X,U=NULL,V=NULL,
                            initXi=NULL,initEta=NULL,
                            tauSeq_a = seq(0.3,0.4,length = 5),
                            tauSeq_b = seq(0.3,0.4,length = 5),
                            tauSeq_xi = seq(0.05, 0.06,length = 5),
                            tauSeq_eta = seq(0.05, 0.06,length = 5),
                            rSeqGlob=c(10,2), rSeqLoc=c(0.2,0.05),
                            verbose=FALSE, doInference = TRUE){
  # dimensions
  p <- dim(X)[1]
  n <- dim(X)[3]

  # additional parameter checking
  if(!(length(rSeqGlob)==length(rSeqLoc))){
    stop('Global and local refinement control sequences must be the same length')
  }
  if((!is.null(initXi) & !(length(initXi)==p))){
    stop('Incorrect dimension for initXi')
  }
  if((!is.null(initEta) & !(length(initEta)==p))){
    stop('Incorrect dimension for initEta')
  }

  # transitivity statistics
  statsGlobal <- NULL
  # U: common neighbours
  if(is.null(U)){
    statsGlobal <- statsTransitivity(X)
    # likelihood uses clipped U without the final snapshot
    Uc <- statsGlobal$U[,,-n]
  }else{
    if(all(dim(U)==c(p,p,n))){
      Uc <- U[,,-n]
    }
    else if(all(dim(U)==c(p,p,n-1))){
      Uc <- U
    }else{
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
  }else{
    if(all(dim(V)==c(p,p,n))){
      Vc <- V[,,-n]
    }
    else if(all(dim(V)==c(p,p,n-1))){
      Vc <- V
    }else{
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
  maxIter <- 10
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

  # Initialization:
  if (is.null(initXi) | is.null(initEta)){
    xiE = rep(1,p)  # If NA, create a vector of 1's of length p
    xiME = matrix(1,p,p)
    etaE = rep(1,p)  # If NA, create a vector of 1's of length p
    etaME = matrix(1,p,p)

    # initialize all xis, etas as 1
    tmp0 = stats::optim(c(tol,tol), globalMLE_ab_et, gr = grr_globalMLE_ab_et, method = "L-BFGS-B",
                        lower = c(0,0), upper=rep(ab_max,2),
                        A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, xivec = xiE, etavec = etaE)
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
          xiij = stats::optim(xiME[i,j], localMLE_init_et, method = 'L-BFGS-B', Aij = A1[i,j,], Bij = B1[i,j,], Uij=Uc[i,j,], Vij=Vc[i,j,], ab = ab1, lower  = c(0.0001), upper = c(ximax[i,j]))$par
          xiME[j,i] = xiME[i,j] = xiij

          etaij = stats::optim(etaME[i,j], localMLE_init_et, method = 'L-BFGS-B', Aij = A2[i,j,], Bij = B2[i,j,], Uij=Vc[i,j,], Vij=Uc[i,j,], ab = c(ab1[2],ab1[1]), lower  = c(0.0001), upper = c(etamax[i,j]))$par
          etaME[j,i] = etaME[i,j] = etaij

        }
      }

      xiE = locEst(xiME)
      etaE = locEst(etaME)
      xiME = outer(xiE,xiE)
      etaME = outer(etaE,etaE)

      tmp2 = stats::optim(ab1, globalMLE_ab_et, gr = grr_globalMLE_ab_et, method = "L-BFGS-B",
                          lower = c(0.01,0.01), upper=rep(ab_max,2),
                          A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, xivec = xiE, etavec = etaE)
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
    # additional control for estimation of a,b
    ab_max_est <- min(10*max(ab1),ab_max)

    for (it in 1:maxIter){
      ximax = apply((1 +exp(ab1[1]*Uc) + exp(ab1[2]*Vc) )/exp(ab1[1]*Uc),c(1,2),min)
      etamax =  apply((1 +exp(ab1[2]*Vc) + exp(ab1[1]*Uc) )/exp(ab1[2]*Vc),c(1,2),min)

      for (i in 1:p){
        xiE[i] = stats::optim(xiE[i],  localMLE_et,  method = 'L-BFGS-B', Ai = A1[i,-i,] , Bi = B1[i,-i,], Ui = Uc[i,-i,], Vi=Vc[i,-i,], ab=ab1, xivec_ic= xiE[-i], lower  = c(0.01), upper = min(ximax[i,-i]/xiE[-i]))$par
        etaE[i] = stats::optim(etaE[i],  localMLE_et,  method = 'L-BFGS-B', Ai = A2[i,-i,] , Bi = B2[i,-i,], Ui = Vc[i,-i,], Vi=Uc[i,-i,], ab=c(ab1[2],ab1[1]), xivec_ic= etaE[-i], lower  = c(0.01), upper = min(etamax[i,-i]/etaE[-i]))$par
      }


      tmp2 = stats::optim(ab1, globalMLE_ab_et, gr = grr_globalMLE_ab_et, method = "L-BFGS-B",
                          lower = c(tol,tol), upper=rep(ab_max_est,2),
                          A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, xivec = xiE, etavec = etaE)
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
  }else if (length(initXi) != p){
    stop("Error: Length of initXi does not match p.")
  }else if (length(initEta) != p){
    stop("Error: Length of initEta does not match p.")
  }else{
    xiE = initXi
    etaE = initEta
    ab1 = c(1,1)
    ab_max_est = ab_max
  }


  if(verbose){
    cat('Done initialization\n')
  }



  ##### Initial value: ab1, xiE and etaE
  ##### Now we apply ONE ITERATION OF the global MLE for (a,b) and local MLE for xii and etai #####
  # global MLE
  ab1  = stats::optim(ab1, globalMLE_ab_et, gr = grr_globalMLE_ab_et, method = "L-BFGS-B",
                      lower = c(tol,tol), upper=rep(ab_max_est,2),
                      A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, xivec = xiE, etavec = etaE)$par
  if(verbose){
    cat('Current (a,b):',ab1,'\nDone global optimization\n')
  }

  ximax = apply((1 +exp(ab1[1]*Uc) + exp(ab1[2]*Vc) )/exp(ab1[1]*Uc),c(1,2),min)
  etamax =  apply((1 +exp(ab1[2]*Vc) + exp(ab1[1]*Uc) )/exp(ab1[2]*Vc),c(1,2),min)

  # local MLE
  for (i in 1:p){
    xiE[i] = stats::optim(xiE[i],  localMLE_et,  method = 'L-BFGS-B', Ai = A1[i,-i,] , Bi = B1[i,-i,], Ui = Uc[i,-i,], Vi=Vc[i,-i,], ab=ab1, xivec_ic= xiE[-i], lower  = c(0.01), upper = min(ximax[i,-i]/xiE[-i]))$par
    etaE[i] = stats::optim(etaE[i],  localMLE_et,  method = 'L-BFGS-B', Ai = A2[i,-i,] , Bi = B2[i,-i,], Ui = Vc[i,-i,], Vi=Uc[i,-i,], ab=c(ab1[2],ab1[1]), xivec_ic= etaE[-i], lower  = c(0.01), upper = min(etamax[i,-i]/etaE[-i]))$par
  }
  if(verbose){
    cat('Current xi quartiles:',stats::quantile(xiE),
        '\nCurrent eta quartiles:',stats::quantile(etaE),'\nDone local optimization\n')
  }


  res <- list()
  res$gVal_init <- ab1
  res$xi_init <- xiE
  res$eta_init <- etaE


  ##### Now we apply IMoM Estimation #####
  stats_alpha = array(c(Uc, Vc), dim = c(p,p,n-1,2))

  IMoM = estNet_trs(X, stats_alpha, stats_alpha, ab1, ab1, shrGPrm = 2, initXi = xiE, initEta = etaE, maxIter = 100)
  rm(stats_alpha)
  xiE = IMoM$xi
  etaE = IMoM$eta
  ab1 = IMoM$gAlphaVal


  ##### Now we apply Refinement Estimation #####
  G = 2
  Gc = 2*p
  q = G + Gc
  S_g = p*(p-1)
  S_gc = p-1
  c_g2 = 0.01*max(q *log(n*S_g)/sqrt(n*S_g),q^{3/2} *(log(n*S_g))^{3/2}/sqrt(n)/S_g)
  c_gc2 = 0.01*max(q *log(n * S_gc)/sqrt(n * S_gc),q^{3/2} *(log(n*S_gc))^{3/2}/sqrt(n)/S_gc)
  delta_n_sqrt = sqrt(max(c_g2, c_gc2))

  # initialize
  ab2 = ab_tmp = ab1
  xiR = xi_tmp =  xiE
  etaR = eta_tmp = etaE

  rtildeGlob <- rSeqGlob[1]
  rtildeLoc <- rSeqLoc[1]


  # initialize shared objects
  g_global <- matrix(0, q, q)
  Al <- matrix(0, q, q)
  El <- diag(1, q)
  tau_xi_seq = tau_eta_seq = rep(0, p)
  tau_a = tau_b = 0

  for (i in 1:p){
    g_i <- matrix(0, q, q)
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
      g_i[1,1] = g_i[1,1] + dl2dada_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], xii, xij, etai, etaj)
      g_i[1,2] = g_i[2,1] = g_i[1,2] + dl2dadb_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], xii, xij, etai, etaj)
      g_i[2,2] = g_i[2,2] + dl2dbdb_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], xii, xij, etai, etaj)

      g_i[1,2+i] =  g_i[2+i,1] = g_i[1,2+i] + dl2dadxii_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xii, xij)
      g_i[1,2+j] =  g_i[2+j,1] = g_i[1,2+j] + dl2dadxii_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xij, xii)
      g_i[2,2+i] =  g_i[2+i,2] = g_i[2,2+i] + dl2dbdxii_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xii, xij)
      g_i[2,2+j] =  g_i[2+j,2] = g_i[2,2+j] + dl2dbdxii_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xij, xii)

      g_i[1,2+p+i] = g_i[2+p+i,1] = g_i[1,2+p+i] + dl2dadetai_ab_et(A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], etai, etaj)
      g_i[1,2+p+j] = g_i[2+p+j,1] = g_i[1,2+p+j] + dl2dadetai_ab_et(A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], etaj, etai)
      g_i[2,2+p+i] = g_i[2+p+i,2] = g_i[2,2+p+i] + dl2dbdetai_ab_et(A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], etai, etaj)
      g_i[2,2+p+j] = g_i[2+p+j,2] = g_i[2,2+p+j] + dl2dbdetai_ab_et(A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], etaj, etai)

      g_i[2+i,2+i] = g_i[2+i,2+i] + dl2dxiidxii_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xii, xij)
      g_i[2+j,2+j] = g_i[2+j,2+j] + dl2dxiidxii_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xij, xii)
      g_i[2+i,2+j] = g_i[2+j,2+i] =  g_i[2+i,2+j] + dl2dxiidxij_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xii, xij)

      g_i[2+p+i,2+p+i] = g_i[2+p+i,2+p+i] + dl2dxiidxii_et(A2ij, B2ij, Vij, Uij, ab1[2], ab1[1], etai, etaj)
      g_i[2+p+j,2+p+j] =  g_i[2+p+j,2+p+j] +dl2dxiidxii_et(A2ij, B2ij, Vij, Uij, ab1[2], ab1[1], etaj, etai)
      g_i[2+p+i,2+p+j] = g_i[2+p+j,2+p+i] = g_i[2+p+i,2+p+j] +dl2dxiidxij_et(A2ij, B2ij, Vij, Uij, ab1[2], ab1[1], etai, etaj)
    }

    g_i <- g_i / (p - 1)
    g_global <- g_global + g_i

    # Select tuning parameter tau for each local parameters:
    # Starting with xi
    var_xi = rep(NA,length(tauSeq_xi))
    var_eta = rep(NA,length(tauSeq_eta))
    Al_xi = matrix(0,q,length(tauSeq_xi))
    Al_eta = matrix(0,q,length(tauSeq_eta))

    if (length(tauSeq_xi) == 1){
      Al[, 2 + i] <- alSearch_et(g_i, El[, 2 + i], tauSeq_xi * delta_n_sqrt)
    }else{
      for (ixi in 1:length(tauSeq_xi)){
        tau_xi = tauSeq_xi[ixi] * delta_n_sqrt
        xi_tmp = xiE
        g_i_tmp = rep(0,q)
        tauvec = array( 0, dim = c(q, n-1))
        var_2 = rep(0,n-1)

        Al_xi[,ixi] = Al[, 2 + i] <- alSearch_et(g_i, El[, 2 + i], tau_xi)
        # refinement for xi
        ximax = apply((1 +exp(ab1[1]*Uc) + exp(ab1[2]*Vc) )/exp(ab1[1]*Uc),c(1,2),min)

        xi_tmp[i] = stats::optim(xiE[i],  localMLE_refine_xi_et, grr_localMLE_refine_xi_et,  method = 'L-BFGS-B', al = Al[,2+i], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, ab = ab1, xivec = xiE, etavec = etaE, i = i, lower = max(0.01, xiE[i]-rtildeLoc), upper = min(c(ximax[i,-i]/xiE[-i], xiE[i]+rtildeLoc)))$par

        # check the variance of \check \xi_l
        for (j in (1:p)[-i]){
          A1ij = A1[i,j,]
          B1ij = B1[i,j,]
          A2ij = A2[i,j,]
          B2ij = B2[i,j,]
          Uij= Uc[i,j,]
          Vij= Vc[i,j,]
          xii = xi_tmp[i]
          xij = xi_tmp[j]
          etai = etaE[i]
          etaj = etaE[j]

          g_i_tmp[1] = g_i_tmp[1] + dl2dadxii_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xii, xij)
          g_i_tmp[2] = g_i_tmp[2] + dl2dbdxii_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xii, xij)
          g_i_tmp[2+i] = g_i_tmp[2+i] + dl2dxiidxii_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xii, xij)
          g_i_tmp[2+j] = g_i_tmp[2+j] + dl2dxiidxij_et(A1ij, B1ij, Uij, Vij, ab1[1], ab1[2], xii, xij)

          Xij_full = X[i,j,]
          Xij = Xij_full[-n]
          Uij = Uc[i,j,]
          Vij = Vc[i,j,]

          tmp = gt2(Xij_full, Uij, Vij, ab1[1], ab1[2], xii, xij, etai, etaj)
          tauvec[1,] = dgammada(Xij, Uij, Vij, ab1[1], ab1[2], xii, xij, etai, etaj)
          tauvec[2,] = dgammadb(Xij, Uij, Vij, ab1[1], ab1[2], xii, xij, etai, etaj)
          tauvec[2+i,] = dgammadxii(Xij, Uij, Vij, ab1[1], ab1[2], xij)
          tauvec[2+j,] = dgammadxii(Xij, Uij, Vij, ab1[1], ab1[2], xii)
          tauvec[2+p+i,] = dgammadetai(Xij, Uij, Vij, ab1[1], ab1[2], etaj)
          tauvec[2+p+j,] = dgammadetai(Xij, Uij, Vij, ab1[1], ab1[2], etai)

          var_2 = var_2 + (Al[,2+i]%*%tauvec)*tmp/(p-1)
        }

        var_xi[ixi] = mean(var_2^2)/(sum(Al[, 2 + i]*g_i_tmp))^2

      }

      tau_xi_seq[i] = tauSeq_xi[which.min(round(sqrt(var_xi),3))]
      Al[, 2 + i] <- alSearch_et(g_i, El[, 2 + i], tau_xi_seq[i]* delta_n_sqrt)
    }

    # for eta
    if (length(tauSeq_eta) == 1){
      Al[, 2 + p + i] <- alSearch_et(g_i, El[, 2 + p + i], tauSeq_eta* delta_n_sqrt)
    }else{
      for (ieta in 1:length(tauSeq_eta)){
        tau_eta = tauSeq_eta[ieta] * delta_n_sqrt
        eta_tmp = etaE
        g_i_tmp = rep(0,q)
        tauvec = array( 0, dim = c(q, n-1))
        var_2 = rep(0,n-1)

        Al_eta[,ieta] = Al[, 2 + p + i] <- alSearch_et(g_i, El[, 2 + p + i], tau_eta)
        # refinement for eta
        etamax =  apply((1 +exp(ab1[2]*Vc) + exp(ab1[1]*Uc) )/exp(ab1[2]*Vc),c(1,2),min)
        eta_tmp[i] = stats::optim(etaE[i],  localMLE_refine_eta_et, grr_localMLE_refine_eta_et,  method = 'L-BFGS-B', al = Al[,2+p+i], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, ab = ab1, xivec = xiE, etavec = etaE, i = i, lower = max(0.01, etaE[i]-rtildeLoc), upper = min(c(etamax[i,-i]/etaE[-i], etaE[i]+rtildeLoc)))$par

        # check the variance of \check \eta_l
        for (j in (1:p)[-i]){
          A1ij = A1[i,j,]
          B1ij = B1[i,j,]
          A2ij = A2[i,j,]
          B2ij = B2[i,j,]
          Uij= Uc[i,j,]
          Vij= Vc[i,j,]
          xii = xiE[i]
          xij = xiE[j]
          etai = eta_tmp[i]
          etaj = eta_tmp[j]

          g_i_tmp[1] = g_i_tmp[1] + dl2dadetai_ab_et(A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], etai, etaj)
          g_i_tmp[2] = g_i_tmp[2] +  dl2dbdetai_ab_et(A2ij, B2ij, Uij, Vij, ab1[1], ab1[2], etai, etaj)
          g_i_tmp[2+p+i] = g_i_tmp[2+p+i] + dl2dxiidxii_et(A2ij, B2ij, Vij, Uij, ab1[2], ab1[1], etai, etaj)
          g_i_tmp[2+p+j] = g_i_tmp[2+p+j] + dl2dxiidxij_et(A2ij, B2ij, Vij, Uij, ab1[2], ab1[1], etai, etaj)

          Xij_full = X[i,j,]
          Xij = Xij_full[-n]
          Uij = Uc[i,j,]
          Vij = Vc[i,j,]

          tmp = gt2(Xij_full, Uij, Vij, ab1[1], ab1[2], xii, xij, etai, etaj)
          tauvec[1,] = dgammada(Xij, Uij, Vij, ab1[1], ab1[2], xii, xij, etai, etaj)
          tauvec[2,] = dgammadb(Xij, Uij, Vij, ab1[1], ab1[2], xii, xij, etai, etaj)
          tauvec[2+i,] = dgammadxii(Xij, Uij, Vij, ab1[1], ab1[2], xij)
          tauvec[2+j,] = dgammadxii(Xij, Uij, Vij, ab1[1], ab1[2], xii)
          tauvec[2+p+i,] = dgammadetai(Xij, Uij, Vij, ab1[1], ab1[2], etaj)
          tauvec[2+p+j,] = dgammadetai(Xij, Uij, Vij, ab1[1], ab1[2], etai)

          var_2 = var_2 + (Al[,2+p+i]%*%tauvec)*tmp/(p-1)
        }
        var_eta[ieta] = mean(var_2^2)/(sum(Al[, 2 +p+ i]*g_i_tmp))^2
      }

      tau_eta_seq[i] = tauSeq_eta[which.min(round(sqrt(var_eta),3))]
      Al[, 2 + p + i] <- alSearch_et(g_i, El[, 2 + p + i], tau_eta_seq[i]* delta_n_sqrt)

    }
  }






  # Final gradient matrix
  g_global <- g_global/p
  var_a = rep(NA,length(tauSeq_a))
  var_b = rep(NA,length(tauSeq_b))
  Al_a = matrix(0,q,length(tauSeq_a))
  Al_b = matrix(0,q,length(tauSeq_b))


  if (length(tauSeq_b) == 1){
    Al[,2] = alSearch_et(g_global, El[,2], tauSeq_b*delta_n_sqrt)
  } else {
    for (ib in 1:length(tauSeq_b)){
      tau_b = tauSeq_b[ib] * delta_n_sqrt
      ab_tmp = ab1
      g_i_tmp = rep(0,q)
      tauvec = array( 0, dim = c(q, n-1))
      var_2 = rep(0,n-1)

      # refinement for b
      Al_b[,ib]= Al[,2] = alSearch_et(g_global, El[,2], tau_b)
      ab_tmp[2] = stats::optim(ab1[2],  globalMLE_refine_b_et,  method = 'L-BFGS-B', al = Al[,2], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, ab = ab1, xivec = xiE, etavec = etaE, lower = max(0.01, ab1[2]-rtildeGlob), upper = ab1[2]+rtildeGlob)$par

      # check the variance of \check b
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

          g_i_tmp[1] = g_i_tmp[1] + dl2dadb_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, ab_tmp[1], ab_tmp[2], xii, xij, etai, etaj)
          g_i_tmp[2] = g_i_tmp[2] + dl2dbdb_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, ab_tmp[1], ab_tmp[2], xii, xij, etai, etaj)
          g_i_tmp[2+i] = g_i_tmp[2+i] + dl2dbdxii_et(A1ij, B1ij, Uij, Vij, ab_tmp[1], ab_tmp[2], xii, xij)
          g_i_tmp[2+j] = g_i_tmp[2+j] + dl2dbdxii_et(A1ij, B1ij, Uij, Vij, ab_tmp[1], ab_tmp[2], xij, xii)
          g_i_tmp[2+p+i] = g_i_tmp[2+p+i] + dl2dbdetai_ab_et(A2ij, B2ij, Uij, Vij, ab_tmp[1], ab_tmp[2], etai, etaj)
          g_i_tmp[2+p+j] = g_i_tmp[2+p+j] + dl2dbdetai_ab_et(A2ij, B2ij, Uij, Vij, ab_tmp[1], ab_tmp[2], etaj, etai)

          Xij_full = X[i,j,]
          Xij = Xij_full[-n]
          Uij = Uc[i,j,]
          Vij = Vc[i,j,]

          tmp = gt2(Xij_full, Uij, Vij, ab_tmp[1], ab_tmp[2], xii, xij, etai, etaj)
          tauvec[1,] = dgammada(Xij, Uij, Vij, ab_tmp[1], ab_tmp[2], xii, xij, etai, etaj)
          tauvec[2,] = dgammadb(Xij, Uij, Vij, ab_tmp[1], ab_tmp[2], xii, xij, etai, etaj)
          tauvec[2+i,] = dgammadxii(Xij, Uij, Vij, ab_tmp[1], ab_tmp[2], xij)
          tauvec[2+j,] = dgammadxii(Xij, Uij, Vij, ab_tmp[1], ab_tmp[2], xii)
          tauvec[2+p+i,] = dgammadetai(Xij, Uij, Vij, ab_tmp[1], ab_tmp[2], etaj)
          tauvec[2+p+j,] = dgammadetai(Xij, Uij, Vij, ab_tmp[1], ab_tmp[2], etai)

          var_2 = var_2 + (Al[,2]%*%tauvec)*tmp/((p - 1)*p)
        }
      }
      g_i_tmp <-g_i_tmp /((p - 1)*p)
      var_b[ib] = mean(var_2^2)/(sum(Al[, 2]*g_i_tmp))^2
    }

    tau_b = tauSeq_b[which.min(round(sqrt(var_b)/100,2))]
    Al[,2] = alSearch_et(g_global, El[,2], tau_b*delta_n_sqrt)

  }

  if(length(tauSeq_a) == 1){
    Al[,1] = alSearch_et(g_global, El[,1], tauSeq_a* delta_n_sqrt)
  }else{
    for (ia in 1:length(tauSeq_a)){
      tau_a = tauSeq_a[ia] * delta_n_sqrt
      ab_tmp = ab1
      g_i_tmp = rep(0,q)
      tauvec = array( 0, dim = c(q, n-1))
      var_2 = rep(0,n-1)

      # refinement for a
      Al_a[,ia] = Al[,1] = alSearch_et(g_global, El[,1], tau_a)
      ab_tmp[1] = stats::optim(ab1[1],  globalMLE_refine_a_et,  method = 'L-BFGS-B', al = Al[,1], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, ab = ab1, xivec = xiE, etavec = etaE, lower = max(0.01, ab1[1]-rtildeGlob), upper = ab1[1]+rtildeGlob)$par

      # check the variance of \check a
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
          g_i_tmp[1] = g_i_tmp[1] + dl2dada_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, ab_tmp[1], ab_tmp[2], xii, xij, etai, etaj)
          g_i_tmp[2] = g_i_tmp[2] + dl2dadb_ab_et(A1ij, B1ij, A2ij, B2ij, Uij, Vij, ab_tmp[1], ab_tmp[2], xii, xij, etai, etaj)
          g_i_tmp[2+i] = g_i_tmp[2+i] + dl2dadxii_et(A1ij, B1ij, Uij, Vij, ab_tmp[1], ab_tmp[2], xii, xij)
          g_i_tmp[2+j] = g_i_tmp[2+j] + dl2dadxii_et(A1ij, B1ij, Uij, Vij, ab_tmp[1], ab_tmp[2], xij, xii)
          g_i_tmp[2+p+i] = g_i_tmp[2+p+i] + dl2dadetai_ab_et(A2ij, B2ij, Uij, Vij, ab_tmp[1], ab_tmp[2], etai, etaj)
          g_i_tmp[2+p+j] = g_i_tmp[2+p+j] + dl2dadetai_ab_et(A2ij, B2ij, Uij, Vij, ab_tmp[1], ab_tmp[2], etaj, etai)

          Xij_full = X[i,j,]
          Xij = Xij_full[-n]
          Uij = Uc[i,j,]
          Vij = Vc[i,j,]

          tmp = gt2(Xij_full, Uij, Vij, ab_tmp[1], ab_tmp[2], xii, xij, etai, etaj)
          tauvec[1,] = dgammada(Xij, Uij, Vij, ab_tmp[1], ab_tmp[2], xii, xij, etai, etaj)
          tauvec[2,] = dgammadb(Xij, Uij, Vij, ab_tmp[1], ab_tmp[2], xii, xij, etai, etaj)
          tauvec[2+i,] = dgammadxii(Xij, Uij, Vij, ab_tmp[1], ab_tmp[2], xij)
          tauvec[2+j,] = dgammadxii(Xij, Uij, Vij, ab_tmp[1], ab_tmp[2], xii)
          tauvec[2+p+i,] = dgammadetai(Xij, Uij, Vij, ab_tmp[1], ab_tmp[2], etaj)
          tauvec[2+p+j,] = dgammadetai(Xij, Uij, Vij, ab_tmp[1], ab_tmp[2], etai)

          var_2 = var_2 + (Al[,1]%*%tauvec)*tmp/((p - 1)*p)
        }
      }

      g_i_tmp <-g_i_tmp /((p - 1)*p)
      var_a[ia] = mean(var_2^2)/(sum(Al[, 1]*g_i_tmp))^2
    }

    tau_a = tauSeq_a[which.min(round(sqrt(var_a)/100,2))]
    Al[,1] = alSearch_et(g_global, El[,1], tau_a* delta_n_sqrt)
  }




  # With the finalized Al, we perform the refinement estimation
  for (rr in 1:length(rSeqGlob)){
    # set refinement control parameters
    rtildeGlob <- rSeqGlob[rr]
    rtildeLoc <- rSeqLoc[rr]

    # refinement for a,b
    ab2[1] = stats::optim(ab1[1],  globalMLE_refine_a_et,  method = 'L-BFGS-B', al = Al[,1], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, ab = ab1, xivec = xiE, etavec = etaE, lower = max(0.01, ab1[1]-rtildeGlob), upper = ab1[1]+rtildeGlob)$par
    ab2[2] = stats::optim(ab1[2],  globalMLE_refine_b_et,  method = 'L-BFGS-B', al = Al[,2], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, ab = ab1, xivec = xiE, etavec = etaE, lower = max(0.01, ab1[2]-rtildeGlob/2), upper = ab1[2]+rtildeGlob)$par
    if(verbose){
      cat('Current (a,b):',ab2,'\nDone global parameter refinement',rr,'of',length(rSeqGlob),'\n')
    }
    # refinement for xi,eta
    ximax = apply((1 +exp(ab1[1]*Uc) + exp(ab1[2]*Vc) )/exp(ab1[1]*Uc),c(1,2),min)
    etamax =  apply((1 +exp(ab1[2]*Vc) + exp(ab1[1]*Uc) )/exp(ab1[2]*Vc),c(1,2),min)

    for (i in 1:p){
      xiR[i] = stats::optim(xiE[i],  localMLE_refine_xi_et, grr_localMLE_refine_xi_et,  method = 'L-BFGS-B', al = Al[,2+i], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, ab = ab1, xivec = xiE, etavec = etaE, i = i, lower = max(0.01, xiE[i]-rtildeLoc), upper = min(c(ximax[i,-i]/xiE[-i], xiE[i]+rtildeLoc)))$par
      etaR[i] = stats::optim(etaE[i],  localMLE_refine_eta_et, grr_localMLE_refine_eta_et,  method = 'L-BFGS-B', al = Al[,2+p+i], A1 = A1, B1 = B1, A2 = A2, B2 = B2, U = Uc, V = Vc, ab = ab1, xivec = xiE, etavec = etaE, i = i, lower = max(0.01, etaE[i]-rtildeLoc), upper = min(c(etamax[i,-i]/etaE[-i], etaE[i]+rtildeLoc)))$par
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
  res$gVal <- ab1
  res$xi <- xiE
  res$eta <- etaE


  ########Inference #######
  if (doInference){
    coef_inf = rep(0,2*p+2)
    ## Global parameter: a and b
    ## for a,b
    zeta_sum_a = zeta_sum_b = 0
    # for local parameter
    zeta_sum = matrix(0, nrow = 2, ncol = p)
    gammavec = array( 0, dim = c(2*p+2, n-1))

    for (i in 1:p){
      for (j in (1:p)[-i]){
        Xij = X[i,j,-n]
        Uij= Uc[i,j,]
        Vij= Vc[i,j,]
        xii = xiE[i]
        xij = xiE[j]
        etai = etaE[i]
        etaj = etaE[j]
        tmp = gamma_root(Xij, Uij, Vij, ab1[1], ab1[2], xii, xij, etai, etaj)
        gammavec[1,] = dgammada(Xij, Uij, Vij, ab1[1], ab1[2], xii, xij, etai, etaj)/tmp
        gammavec[2,] = dgammadb(Xij, Uij, Vij, ab1[1], ab1[2], xii, xij, etai, etaj)/tmp
        gammavec[2+i,] = dgammadxii(Xij, Uij, Vij, ab1[1], ab1[2], xij)/tmp
        gammavec[2+j,] = dgammadxii(Xij, Uij, Vij, ab1[1], ab1[2], xii)/tmp
        gammavec[2+p+i,] = dgammadetai(Xij, Uij, Vij, ab1[1], ab1[2], etaj)/tmp
        gammavec[2+p+j,] = dgammadetai(Xij, Uij, Vij, ab1[1], ab1[2], etai)/tmp
        zeta_sum_a = zeta_sum_a + mean((Al[,1]%*%gammavec)^2)
        zeta_sum_b = zeta_sum_b + mean((Al[,2]%*%gammavec)^2)
        zeta_sum[1,i] = zeta_sum[1,i] + mean((Al[,2+i]%*%gammavec)^2)
        zeta_sum[2,i] = zeta_sum[2,i] + mean((Al[,2+p+i]%*%gammavec)^2)
      }
    }
    zeta_sum_a = zeta_sum_a/(p*(p-1))
    coef_inf[1] = sqrt(n*(p*(p-1)/2)/zeta_sum_a)

    zeta_sum_b = zeta_sum_b/(p*(p-1))
    coef_inf[2] =sqrt(n*(p*(p-1)/2)/zeta_sum_b)

    zeta_sum = zeta_sum/(p-1)
    tmp = sqrt(n*(p-1)/zeta_sum)
    coef_inf[1:p+2] = tmp[1,]
    coef_inf[1:p+2+p] = tmp[2,]

    # store estimators
    res$se_estimates = 1/coef_inf
  }


  return(res)
}
