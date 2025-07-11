#' Estimation for General Autoregressive Networks
#'
#' This function estimates the parameters of a user-specified autoregressive network model
#' using an Iterative Method-of-Moments (IMoM) estimation
#' procedure. The method alternates
#' between optimizing global and local parameters, capturing edge formation and
#' dissolution dynamics in dynamic networks.
#'
#' The network dynamics are modeled through autoregressive equations for the probabilities of edge formation and dissolution between nodes i and j at time t-1:
#' \deqn{\alpha_{i,j}^{t-1} = \xi_i \xi_j f(\mathbf{X_{t-1}},...,\mathbf{X_{t-m_f}}; \mathbf{\theta}_f),}{"alpha_{i,j}^{t-1} = xi_i * xi_j * f(X_{t-1},...,X_{t-m_f}; theta_f),"}
#' \deqn{\beta_{i,j}^{t-1} = \eta_i \eta_j g(\mathbf{X_{t-1}},...,\mathbf{X_{t-m_g}}; \mathbf{\theta}_g),}{"beta_{i,j}^{t-1} = eta_i * eta_j * g(X_{t-1},...,X_{t-m_g}; theta_g),"}
#' where \eqn{\xi_i, \xi_j, \eta_i, \eta_j} are node-specific parameters, and \eqn{f}, \eqn{g} represent functions defining global edge dynamics, parameterized by \eqn{\mathbf{\theta}_f} and \eqn{\mathbf{\theta}_g}, respectively.
#'
#' User-Defined Functions:
#' The `fij` and `gij` functions, modeling global edge formation and dissolution dynamics, must adhere to specific input/output rules:
#' \itemize{
#'   \item \strong{First Argument}: A vector of global parameters, with the first \code{shrGPrm} parameters being shared across both functions.
#'   \item \strong{Second Argument}: A \eqn{p \times p \times (n-1) \times d} array of sufficient statistics from observed networks, where \eqn{d} is the number of sufficient statistics.
#'   \item \strong{Output}: A \eqn{p \times p \times (n-1)} array representing the values \eqn{f(\mathbf{X_{t-1}},...,\mathbf{X_{t-m_f}}; \mathbf{\theta}_f)} or \eqn{g(\mathbf{X_{t-1}},...,\mathbf{X_{t-m_g}}; \mathbf{\theta}_g)} per node pair across time.
#' }
#'
#' Estimation Procedure:
#' \enumerate{
#'   \item \strong{Initialization}: Sets initial values for local parameters.
#'   \item \strong{Global Optimization}: Utilizes initialized local parameters within a global
#'    log-likelihood function, optimizing it through a quasi-Newton method to
#'    estimate global parameters.
#'   \item \strong{Local Optimization}: With newly estimated global parameters, local
#'    log-likelihood functions are maximized through a quasi-Newton method to update local parameter estimates.
#'   \item \strong{Iteration and Convergence}: The process iterates until convergence criteria are met, typically when the mean absolute difference between successive global parameter estimates drops below \code{tol},
#'    or when the global log-likelihood decreases, signaling optimal estimates.
#' }
#'
#'
#' @param X A \eqn{p \times p \times n} array of network's adjacency matrices over time.
#' @param fij Function for global edge formation behaviour, dependent on the global parameters \eqn{\mathbf{\theta}_f} and sufficient statistics built upon past obervations. See also ‘Details’.
#' @param gij Function for global edge dissolution behaviour, dependent on the global parameters \eqn{\mathbf{\theta}_g} and sufficient statistics built upon past obervations. See also ‘Details’.
#' @param statsAlpha Sufficient statistics input for \code{fij}.
#' @param statsBeta Sufficient statistics input for \code{gij}.
#' @param globInitAlpha Initial values for global parameters in \code{fij}.
#' @param globInitBeta Initial values for global parameters in \code{gij}.
#' @param shrGPrm Number of shared global parameters between \code{fij} and
#'   \code{gij}.
#' @param initXi Initial values for \eqn{\xi_1, \dots, \xi_p}. Defaults to a vector of ones if NULL.
#' @param initEta Initial values for \eqn{\eta_1, \dots, \eta_p}. Defaults to a vector of ones if NULL.
#' @param updateMethod Method for updating global parameters ("sequential" or "batch").
#' @param tol Tolerance for the optimization process.
#' @param maxIter Maximum iterations allowed for the optimization algorithm.
#'
#' @return A list containing the estimated parameters:
#' \itemize{
#'   \item \code{gAlphaVal}: Estimated global parameters in \code{fij}.
#'   \item \code{gBetaVal}: Estimated global parameters in \code{gij}.
#'   \item \code{xi}: Estimated values of \eqn{\xi_1, \dots, \xi_p}.
#'   \item \code{eta}: Estimated values of \eqn{\eta_1, \dots, \eta_p}.
#' }
#'
#' @examples
#' # Example 1: transitivity model.
#' p = 30; n = 20
#' xi = rep(0.7, p); eta = rep(0.8, p)
#' a = 30; b = 15
#'
#' # Simulate data using simulateTransitivity function
#' simulated_data = simulateTransitivity(p, n, xi, eta, a, b)
#' X = simulated_data$X
#' U = simulated_data$U[, , 1:(n - 1)]
#' V = simulated_data$V[, , 1:(n - 1)]
#'
#' statsAlpha = statsBeta = array(c(U, V), dim = c(p,p,n-1,2))
#'
#' # Initialize global parameters for fij and gij
#' globInitAlpha = globInitBeta = rep(100, 2)
#'
#' # Define edge formation and dissolution functions based on transitivity model
#' fij = function(global, stats) {
#'   exp(global[1] * stats[, , , 1]) /
#'   (1 + exp(global[1] * stats[, , , 1]) + exp(global[2] * stats[, , , 2]))
#' }
#' gij = function(global, stats) {
#'   exp(global[2] * stats[, , , 2]) /
#'   (1 + exp(global[1] * stats[, , , 1]) + exp(global[2] * stats[, , , 2]))
#' }
#'
#' result = estNet(X, fij, gij, statsAlpha, statsBeta, globInitAlpha, globInitBeta,
#'                 shrGPrm = 2, updateMethod = "sequential", maxIter=2)
#'
#' result = estNet(X, fij, gij, statsAlpha, statsBeta, globInitAlpha, globInitBeta,
#'                 shrGPrm = 2, updateMethod = "batch", maxIter=2)
#' # NOTE: small maxIter used to limit example runtime, default value is 100
#'
#' # Example 2: persistence model.
#' p = 30; n = 20
#' xi = runif(p,0.5,0.9); eta = runif(p,0.5,0.9)
#' a = 0.5; b = 0.5
#'
#' # Simulate data using simulatePersistence function
#' simulated_data = simulatePersistence(p, n, xi, eta, a, b)
#' X = simulated_data$X
#' statsAlpha =  1-X[,,4:n-2]+ ( 1-X[,,4:n-2])*( 1-X[,,4:n-3])
#' statsBeta =  X[,,4:n-2]+ ( X[,,4:n-2])*( X[,,4:n-3])
#' X = X[,, 3: n]
#'
#' # Initialize global parameters for fij and gij
#' globInitAlpha = globInitBeta = 1
#'
#' # Define edge formation and dissolution functions based on persistence model
#' fij = function(global, stats) {
#'   return (exp(-1 -global*stats))
#' }
#' gij <- function(global, stats) {
#'   return (exp(-1 -global*stats))
#' }
#'
#' result <- estNet(X, fij, gij, statsAlpha, statsBeta, globInitAlpha, globInitBeta,
#'                  shrGPrm = 0, maxIter = 2)
#' # NOTE: small maxIter used to limit example runtime, default value is 100
#' # NOTE: The provided examples are for illustration purpose.
#' #       Larger sample, e.g. p = 50 and n = 100, is recommended for better estimation performance.
#'
#' @export
estNet = function(X, fij, gij, statsAlpha, statsBeta, globInitAlpha, globInitBeta, shrGPrm,
                  initXi = NULL, initEta = NULL, updateMethod="sequential", tol = 0.01, maxIter = 100){
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

  # Validate fij and gij function structure
  if (!is.function(fij) || !is.function(gij)) {
    stop("fij and gij must be functions.")
  }

  # Test fij function structure
  tryCatch({
    testOutput <- fij(globInitAlpha, statsAlpha)
    if (!(is.array(testOutput) && all(dim(testOutput) == c(p, p, n - 1)))) {
      stop("The output of fij does not meet the expected criteria: it should be an array with dimensions (p * p * (n - 1)).")
    }

  }, error = function(e) {
    stop("Error in testing fij function structure: ", e$message)
  })
  # Test gij function structure
  tryCatch({
    testOutput <- gij(globInitBeta, statsBeta)
    if (!(is.array(testOutput) && all(dim(testOutput) == c(p, p, n - 1)))) {
      stop("The output of gij does not meet the expected criteria: it should be an array with dimensions (p * p * (n - 1)).")
    }

  }, error = function(e) {
    stop("Error in testing gij function structure: ", e$message)
  })


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


  # Iteratively updates global parameters and local parameters in sequence, adjusting one parameter at a time to optimize the model.
  if(updateMethod=="sequential"){
    if (shrGPrm == 0){
      for (ix in 1:da){
        tmp0 = stats::optim(gAlphaVal.E0[ix], globalMLE, method = "L-BFGS-B", lower = c(0.01), ix = ix, A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateAlpha = TRUE)
        gAlphaVal.E0[ix] =  tmp0$par
      }

      fn1 = - logl_part(A1, B1, fij, gAlphaVal.E0, statsAlpha, xiME)

      for (it in 1: maxIter){
        xiMax = apply(1/fg_array(fij, gAlphaVal.E0, statsAlpha) , c(1,2), min)
        for (i in 1:(p-1)){
          for (j in (i+1):p){
            xiij = stats::optim(xiME[i,j], localMLE, method = 'L-BFGS-B', Aij = A1[i, j,], Bij = B1[i,j,], fg = fij, global = gAlphaVal.E0, stats = statsAlpha[i,j, , ,drop = FALSE], lower  = c(0.01), upper = c(xiMax[i,j]))$par
            xiME[j,i] = xiME[i,j] = xiij
          }
        }

        xiE = locEst(xiME)
        xiME = tcrossprod(xiE)

        for (ix in 1:da){
          tmp2 = stats::optim(gAlphaVal.E0[ix], globalMLE, method = "L-BFGS-B", lower = c(0.01), ix = ix, A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateAlpha = TRUE)
          gAlphaVal.E[ix] = tmp2$par
        }

        fn2 = - logl_part(A1, B1, fij, gAlphaVal.E, statsAlpha, xiME)


        if(mean(abs(gAlphaVal.E - gAlphaVal.E0))<tol| fn1<=fn2) {
          break
        }else{
          gAlphaVal.E0 = gAlphaVal.E
          fn1 = fn2
        }
      }

      for (it in 1: maxIter){
        xiMax = apply(1/fg_array(fij, gAlphaVal.E0, statsAlpha) , c(1,2), min)
        for (i in 1:p){
          xiE[i] = stats::optim(xiE[i],  localMLE_i,  method = 'L-BFGS-B', A = A1, B = B1, fg = fij, global = gAlphaVal.E0, stats = statsAlpha, localVec = xiE, i = i, lower  = c(0.01), upper = min(xiMax[i,-i]/xiE[-i]))$par
        }
        xiME = outer(xiE,xiE)
        for (ix in 1:da){
          tmp2 = stats::optim(gAlphaVal.E0[ix], globalMLE, method = "L-BFGS-B", lower = c(0.01), ix = ix, A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateAlpha = TRUE)
          gAlphaVal.E[ix] = tmp2$par
        }
        fn2 = - logl_part(A1, B1, fij, gAlphaVal.E, statsAlpha, xiME)

        if(mean(abs(gAlphaVal.E - gAlphaVal.E0))<tol| fn1<=fn2) {
          break
        }else{
          gAlphaVal.E0 = gAlphaVal.E
          fn1 = fn2
        }
      }


      for (ix in 1:db){
        tmp0 = stats::optim(gBetaVal.E0[ix], globalMLE, method = "L-BFGS-B", lower = c(0.01), ix = ix, A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateAlpha = F)
        gBetaVal.E0[ix] =  tmp0$par
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
          tmp2 = stats::optim(gBetaVal.E0[ix], globalMLE, method = "L-BFGS-B", lower = c(0.01), ix = ix, A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateAlpha = F)
          gBetaVal.E[ix] = tmp2$par
        }

        fn2 = - logl_part(A2, B2, gij, gBetaVal.E, statsBeta, etaME)


        if(mean(abs(gBetaVal.E - gBetaVal.E0))<tol| fn1<=fn2) {
          break
        }else{
          gBetaVal.E0 = gBetaVal.E
          fn1 = fn2
        }
      }

      for (it in 1: maxIter){
        etaMax = apply(1/fg_array(gij, gBetaVal.E0, statsBeta) , c(1,2),min)

        for (i in 1:p){
          etaE[i] = stats::optim(etaE[i],  localMLE_i,  method = 'L-BFGS-B', A = A2, B = B2, fg = gij, global = gBetaVal.E0, stats = statsBeta, localVec = etaE, i = i, lower  = c(0.01), upper = min(etaMax[i,-i]/etaE[-i]))$par
        }
        etaME = outer(etaE,etaE)

        for (ix in 1:db){
          tmp2 = stats::optim(gBetaVal.E0[ix], globalMLE, method = "L-BFGS-B", lower = c(0.01), ix = ix, A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateAlpha = F)
          gBetaVal.E[ix] = tmp2$par
        }

        fn2 = - logl_part(A2, B2, gij, gBetaVal.E, statsBeta, etaME)


        if(mean(abs(gBetaVal.E - gBetaVal.E0))<tol| fn1<=fn2) {
          break
        }else{
          gBetaVal.E0 = gBetaVal.E
          fn1 = fn2
        }
      }

    }


    if ( shrGPrm>0){
      for (ix in 1:shrGPrm){
        tmp0 = stats::optim(gAlphaVal.E0[ix], globalMLE, method = "L-BFGS-B", lower = c(0.01), ix = ix, A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateAlpha = TRUE)
        gAlphaVal.E0[ix] = gBetaVal.E0[ix] = tmp0$par
      }
    }

    if (da > shrGPrm){
      for (ix in (shrGPrm + 1):length(gAlphaVal.E0)){
        tmp0 = stats::optim(gAlphaVal.E0[ix], globalMLE, method = "L-BFGS-B", lower = c(0.01), ix = ix, A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateAlpha = TRUE)
        gAlphaVal.E0[ix] =  tmp0$par
      }
    }

    if (db > shrGPrm){
      for (ix in (shrGPrm + 1):length(gBetaVal.E0)){
        tmp0 = stats::optim(gBetaVal.E0[ix], globalMLE, method = "L-BFGS-B", lower = c(0.01), ix = ix, A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateAlpha = FALSE)
        gBetaVal.E0[ix] =  tmp0$par
      }
    }

    fn1 = - logl(A1, B1, A2, B2, fij, gAlphaVal.E0, statsAlpha, xiME, gij, gBetaVal.E0, statsBeta, etaME)

    for (it in 1: maxIter){
      xiMax = apply(1/fg_array(fij, gAlphaVal.E0, statsAlpha) , c(1,2), min)
      etaMax = apply(1/fg_array(gij, gBetaVal.E0, statsBeta) , c(1,2), min)
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
        for (ix in 1:shrGPrm){
          tmp2 = stats::optim(gAlphaVal.E0[ix], globalMLE, method = "L-BFGS-B", lower = c(0.01), ix = ix, A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateAlpha = TRUE)
          gAlphaVal.E[ix] = gBetaVal.E[ix] = tmp2$par
        }
      }

      if (da > shrGPrm){
        for (ix in (shrGPrm + 1):length(gAlphaVal.E0)){
          tmp2 = stats::optim(gAlphaVal.E0[ix], globalMLE, method = "L-BFGS-B", lower = c(0.01), ix = ix, A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateAlpha = TRUE)
          gAlphaVal.E[ix] =  tmp2$par
        }
      }

      if (db > shrGPrm){
        for (ix in (shrGPrm + 1):length(gBetaVal.E0)){
          tmp2 = stats::optim(gBetaVal.E0[ix], globalMLE, method = "L-BFGS-B", lower = c(0.01), ix = ix, A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateAlpha = F)
          gBetaVal.E[ix] =  tmp2$par
        }
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

    for (it in 1: maxIter){
      xiMax = apply(1/fg_array(fij, gAlphaVal.E0, statsAlpha) , c(1,2), min)
      etaMax = apply(1/fg_array(gij, gBetaVal.E0, statsBeta) , c(1,2), min)
      for (i in 1:p){
        xiE[i] = stats::optim(xiE[i],  localMLE_i,  method = 'L-BFGS-B', A = A1, B = B1, fg = fij, global = gAlphaVal.E0, stats = statsAlpha, localVec = xiE, i = i, lower  = c(0.01), upper = min(xiMax[i,-i]/xiE[-i]))$par
        etaE[i] = stats::optim(etaE[i],  localMLE_i,  method = 'L-BFGS-B', A = A2, B = B2, fg = gij, global = gBetaVal.E0, stats = statsBeta, localVec = etaE, i = i, lower  = c(0.01), upper = min(etaMax[i,-i]/etaE[-i]))$par
      }


      xiME = tcrossprod(xiE)
      etaME = tcrossprod(etaE)

      if ( shrGPrm>0){
        for (ix in 1:shrGPrm){
          tmp2 = stats::optim(gAlphaVal.E0[ix], globalMLE, method = "L-BFGS-B", lower = c(0.01), ix = ix, A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateAlpha = TRUE)
          gAlphaVal.E[ix] = gBetaVal.E[ix] = tmp2$par
        }
      }

      if (da > shrGPrm){
        for (ix in (shrGPrm + 1):length(gAlphaVal.E0)){
          tmp2 = stats::optim(gAlphaVal.E0[ix], globalMLE, method = "L-BFGS-B", lower = c(0.01), ix = ix, A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateAlpha = TRUE)
          gAlphaVal.E[ix] =  tmp2$par
        }
      }

      if (db > shrGPrm){
        for (ix in (shrGPrm + 1):length(gBetaVal.E0)){
          tmp2 = stats::optim(gBetaVal.E0[ix], globalMLE, method = "L-BFGS-B", lower = c(0.01), ix = ix, A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateAlpha = F)
          gBetaVal.E[ix] =  tmp2$par
        }
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
  }else if (updateMethod=="batch"){
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

      for (it in 1: maxIter){
        xiMax = apply(1/fg_array(fij, gAlphaVal.E0, statsAlpha) , c(1,2), min)
        for (i in 1:p){
          xiE[i] = stats::optim(xiE[i],  localMLE_i,  method = 'L-BFGS-B', A = A1, B = B1, fg = fij, global = gAlphaVal.E0, stats = statsAlpha, localVec = xiE, i = i, lower  = c(0.01), upper = min(xiMax[i,-i]/xiE[-i]))$par
        }
        xiME = outer(xiE,xiE)
        for (ix in 1:da){
          tmp2 = stats::optim(gAlphaVal.E0, globalMLE_group, method = "L-BFGS-B", lower = rep(0.01,da), A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateShare= F,updateAlpha = TRUE)
          gAlphaVal.E = tmp2$par
        }
        fn2 = - logl_part(A1, B1, fij, gAlphaVal.E, statsAlpha, xiME)

        if(mean(abs(gAlphaVal.E - gAlphaVal.E0))<tol| fn1<=fn2) {
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

      for (it in 1: maxIter){
        etaMax = apply(1/fg_array(gij, gBetaVal.E0, statsBeta) , c(1,2),min)

        for (i in 1:p){
          etaE[i] = stats::optim(etaE[i],  localMLE_i,  method = 'L-BFGS-B', A = A2, B = B2, fg = gij, global = gBetaVal.E0, stats = statsBeta, localVec = etaE, i = i, lower  = c(0.01), upper = min(etaMax[i,-i]/etaE[-i]))$par
        }
        etaME = outer(etaE,etaE)

        for (ix in 1:db){
          tmp2 = stats::optim(gBetaVal.E0, globalMLE_group, method = "L-BFGS-B", lower = rep(0.01,db),  A1 = A1, B1 = B1, A2 = A2, B2 = B2, fij = fij, gAlphaVal = gAlphaVal.E0, statsAlpha = statsAlpha, xiMat = xiME, gij = gij, gBetaVal = gBetaVal.E0, statsBeta = statsBeta, etaMat = etaME, shrGPrm = shrGPrm, updateShare = F, updateAlpha = F)
          gBetaVal.E = tmp2$par
        }

        fn2 = - logl_part(A2, B2, gij, gBetaVal.E, statsBeta, etaME)


        if(mean(abs(gBetaVal.E - gBetaVal.E0))<tol| fn1<=fn2) {
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

    for (it in 1: maxIter){
      xiMax = apply(1/fg_array(fij, gAlphaVal.E0, statsAlpha) , c(1,2), min)
      etaMax = apply(1/fg_array(gij, gBetaVal.E0, statsBeta) , c(1,2), min)
      for (i in 1:p){
        xiE[i] = stats::optim(xiE[i],  localMLE_i,  method = 'L-BFGS-B', A = A1, B = B1, fg = fij, global = gAlphaVal.E0, stats = statsAlpha, localVec = xiE, i = i, lower  = c(0.01), upper = min(xiMax[i,-i]/xiE[-i]))$par
        etaE[i] = stats::optim(etaE[i],  localMLE_i,  method = 'L-BFGS-B', A = A2, B = B2, fg = gij, global = gBetaVal.E0, stats = statsBeta, localVec = etaE, i = i, lower  = c(0.01), upper = min(etaMax[i,-i]/etaE[-i]))$par
      }

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
  }else {
    stop("Invalid update method. Choose 'sequential' or 'batch'.")
  }

  res = list()
  res$gAlphaVal = gAlphaVal.E0
  res$gBetaVal = gBetaVal.E0
  res$xi = xiE
  res$eta = etaE

  return(res)

}


