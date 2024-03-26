#' Prediction for general autoregressive network models
#'
#' This function predicts the edge probabilities of the next snapshot of a sequence of
#' dynamic networks based on a model fit using
#' fit using \code{estNet}.
#'
#' @param estimates A list of parameter estimates from \code{estNet}.
#' @param Xnew A \eqn{p \times p} matrix of the network's current adjacency matrix from which to
#' perform one step prediction
#' @param statsAlphaNew A \eqn{p \times p \times 1 \times d} array, current snapshot sufficient statistics
#' input for \code{fij}.
#' @param statsBetaNew A \eqn{p \times p \times 1 \times d} array, current snapshot sufficient statistics
#' input for \code{gij}.
#' @param fij Function for global edge formation behaviour, dependent on the global parameters \eqn{\mathbf{\theta}_f} and sufficient statistics. See ‘Details’ for \code{estNet} documentation.
#' @param gij Function for global edge dissolution behaviour, dependent on the global parameters \eqn{\mathbf{\theta}_g} and sufficient statistics. See ‘Details’ for \code{estNet} documentation.
#'
#' @return A \eqn{p \times p} matrix, the predicted edge probabilites for the next network snapshot.
#'
#' @examples
#' # Example with transitivity model.
#' p = 50; n = 100
#' xi = rep(0.7, p); eta = rep(0.8, p)
#' a = 30; b = 15
#'
#' # Simulate data using simulate_transitivity function
#' simulated_data = simulate_transitivity(p, n, xi, eta, a, b)
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
#'   exp(global[1] * stats[, , , 1,drop=FALSE]) /
#'   (1 + exp(global[1] * stats[, , , 1,drop=FALSE]) + exp(global[2] * stats[, , , 2,drop=FALSE]))
#' }
#' gij = function(global, stats) {
#'   exp(global[2] * stats[, , , 2,drop=FALSE]) /
#'   (1 + exp(global[1] * stats[, , , 1,drop=FALSE]) + exp(global[2] * stats[, , , 2,drop=FALSE]))
#' }
#'
#' result = estNet(X, fij, gij, statsAlpha, statsBeta, globInitAlpha, globInitBeta,
#'                 shrGPrm = 2)
#'
#' Xnew = X[,,n]
#' statsAlphaNew = statsBetaNew = array(c(U[,,n],V[,,n]),c(p,p,1,2))
#' pred = predictNet(result,Xnew,statsAlphaNew,statsBetaNew,fij,gij)
#'
#' @export
predictNet <- function(estimates,
                       Xnew,statsAlphaNew,statsBetaNew,
                       fij,gij){
  # calculate alphas for next snapshot
  Alpha <- tcrossprod(estimates$xi)*fij(estimates$gAlphaVal,statsAlphaNew)[,,1]
  # calculate betas for next snapshot
  Beta <- tcrossprod(estimates$eta)*gij(estimates$gBetaVal,statsBetaNew)[,,1]
  # calculate probability of next edge
  Pred <- Alpha + Xnew*(1 - Alpha - Beta)
  Pred <- Pred - diag(Pred)
  return(Pred)
}

# prediction function for transitivity
predictTransitivity <- function(estimates,
                                Xnew,
                                Unew=NULL,Vnew=NULL,
                                n_step=1){
  # dimensions
  p <- dim(X_prev)[1]
  # initialize current observation
  Xcurr <- Xnew

  ##### to do from here
  # populate current U/V statistics
  U <- V <- matrix(0,p,p)
  # populate U,V
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      tmp <- U[i,j] <- U[j,i] <- sum(X_curr[i,]*X_curr[j,])/(p-1)
      V[i,j] <- V[j,i] <- (sum(X_curr[i,]+X_curr[j,])-2*X_curr[i,j])/(p-1) - 2*tmp
    }
  }
  stats_curr <- list(U=U,V=V)
  # theta/eta degree parameters
  Ttheta <- tcrossprod(fit$theta)
  Eeta <- tcrossprod(fit$eta)
  # initialize array
  X_out <- array(NA,c(p,p,n_out))
  for(t in 1:n_out){
    # predict
    eaU <- exp(fit$a * stats_curr$U)
    ebV <- exp(fit$b * stats_curr$V)
    alpha <- (Ttheta * eaU) / (1 + eaU + ebV)
    beta <- (Eeta * ebV) / (1 + eaU + ebV)
    X_out[,,t] <- alpha + X_curr*(1 - alpha - beta)
    if(t < n_out){
      # update X
      X_curr <- X_out[,,t]
      # update U,V
      U <- V <- matrix(0,p,p)
      # populate U,V
      for(i in 1:(p-1)){
        for(j in (i+1):p){
          tmp <- U[i,j] <- U[j,i] <- sum(X_curr[i,]*X_curr[j,])/(p-1)
          V[i,j] <- V[j,i] <- (sum(X_curr[i,]+X_curr[j,])-2*X_curr[i,j])/(p-1) - 2*tmp
        }
      }
      stats_curr <- list(U=U,V=V)
    }
  }
  if(n_out==1){
    return(X_out[,,1])
  }
  else{
    return(X_out)
  }
}
