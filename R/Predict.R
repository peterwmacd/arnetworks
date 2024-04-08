#' Prediction for General Autoregressive Networks
#'
#' This function predicts the edge probabilities of the next snapshot of a sequence of
#' dynamic networks based on a model
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
#'   p <- dim(stats)[1]
#'   n <- dim(stats)[3]
#'   temp <- exp(global[1] * stats[, , , 1]) /
#'   (1 + exp(global[1] * stats[, , , 1]) + exp(global[2] * stats[, , , 2]))
#'   return(array(temp,c(p,p,n)))
#' }
#' gij = function(global, stats) {
#'   p <- dim(stats)[1]
#'   n <- dim(stats)[3]
#'   temp <- exp(global[2] * stats[, , , 2]) /
#'   (1 + exp(global[1] * stats[, , , 1]) + exp(global[2] * stats[, , , 2]))
#'   return(array(temp,c(p,p,n)))
#' }
#'
#' # parameters in the format of estNet
#' result = list(gAlphaVal=c(25,5),gBetaVal=c(25,5),
#'               xi=rep(0.8,p),eta=rep(0.7,p))
#'
#' Xnew = X[,,n]
#' statsAlphaNew = statsBetaNew = array(c(simulated_data$U[,,n],simulated_data$V[,,n]),c(p,p,1,2))
#' pred = predictNet(result,Xnew,statsAlphaNew,statsBetaNew,fij,gij)
#'
#' @export
predictNet <- function(estimates,
                       Xnew,statsAlphaNew,statsBetaNew,
                       fij,gij){
  # Define the expected names
  expected_names <- c('gAlphaVal', 'gBetaVal', 'xi', 'eta')
  
  # Check if all names in 'estimates' match the expected names
  if (!all(names(estimates) %in% expected_names)) {
    # Create an error message that lists the expected names
    expected_names_str <- paste("'", expected_names, "'", sep = "", collapse = ", ")
    stop(sprintf("Estimate names do not match expected names. Expected names are: %s.", expected_names_str))
  }
  # calculate alphas for next snapshot
  Alpha <- tcrossprod(estimates$xi)*fij(estimates$gAlphaVal,statsAlphaNew)[,,1]
  # calculate betas for next snapshot
  Beta <- tcrossprod(estimates$eta)*gij(estimates$gBetaVal,statsBetaNew)[,,1]
  # calculate probability of next edge
  Pred <- Alpha + Xnew*(1 - Alpha - Beta)
  Pred <- Pred - diag(diag(Pred))
  return(Pred)
}

#' Prediction for Autoregressive Networks with Transitivity
#'
#' This function recursively predicts the edge probabilities of future snapshots of a sequence of
#' dynamic networks based on a model
#' fit using \code{estTransitivity}.
#'
#' @param estimates A list of parameter estimates from \code{estTransitivity}.
#' @param Xnew A \eqn{p \times p} matrix of the network's current adjacency matrix from which to
#' perform prediction
#' @param Unew A \eqn{p \times p \times 1} array or \eqn{p \times p } matrix of the normalized number of
#' common neighbour statistics for the current snapshot. If not provided, it will be calculated internally from \code{Xnew}.
#' @param Vnew A \eqn{p \times p \times 1} array or \eqn{p \times p } matrix of the normalized number of
#' disjoint neighbour statistics for the current snapshot. If not provided, it will be calculated internally from \code{Xnew}.
#' @param nStep A positive integer, the number of steps to recursively predict with the model.
#'
#' @return A \eqn{p \times p} matrix (if \code{nStep=1}) or a \eqn{p \times p \times} \code{nStep} array,
#' the predicted edge probabilites for the next \code{nStep} network snapshots.
#'
#' @examples
#' p = 30; n = 20
#' xi = rep(0.7, p); eta = rep(0.8, p)
#' a = 30; b = 15
#'
#' # Simulate data using simulateTransitivity function
#' simulated_data = simulateTransitivity(p, n, xi, eta, a, b)
#' X = simulated_data$X
#' U = simulated_data$U
#' V = simulated_data$V
#'
#' # Model parameters in format of estTransitivity
#' result = list(gVal=c(25,5), xi=rep(0.8, p), eta=rep(0.7, p))
#'
#' # Predict the next 2 snapshots using predictTransitivity
#' pred = predictTransitivity(result,X[,,n],nStep=2)
#'
#' @export
predictTransitivity <- function(estimates,
                                Xnew,
                                Unew=NULL,Vnew=NULL,
                                nStep=1){
  # dimensions
  p <- dim(Xnew)[1]
  if(nStep < 1){
    stop('nStep must be a positive integer')
  }
  # initialize current observation
  Xcurr <- Xnew
  # initialize statistics
  if(is.null(Unew)){
    statsCurr <- statsTransitivity(array(Xnew,c(p,p,1)))
    # likelihood uses clipped U without the final snapshot
    Ucurr <- statsCurr$U[,,1]
  }
  else{
    if(all(dim(Unew)==c(p,p,1))){
      Ucurr <- Unew[,,1]
    }
    else if(all(dim(Unew)==c(p,p))){
      Ucurr <- Unew
    }
    else{
      stop('Incorrect dimensions for U')
    }
  }
  # V: uncommon neighbours
  if(is.null(Vnew)){
    if(is.null(statsCurr)){
      statsCurr <- statsTransitivity(array(Xnew,c(p,p,1)))
    }
    # likelihood uses clipped U without the final snapshot
    Vcurr <- statsCurr$V[,,1]
  }
  else{
    if(all(dim(Vnew)==c(p,p,1))){
      Vcurr <- Vnew[,,1]
    }
    else if(all(dim(Vnew)==c(p,p))){
      Vcurr <- Vnew
    }
    else{
      stop('Incorrect dimensions for V')
    }
  }
  # populate current U/V statistics
  
  # Define the expected names
  expected_names <- c('gVal', 'xi', 'eta')
  # Check if all names in 'estimates' match the expected names
  if (!all(names(estimates) %in% expected_names)) {
    # Create an error message that lists the expected names
    expected_names_str <- paste("'", expected_names, "'", sep = "", collapse = ", ")
    stop(sprintf("Estimate names do not match expected names. Expected names are: %s.", expected_names_str))
  }
  xiM <- tcrossprod(estimates$xi)
  etaM <- tcrossprod(estimates$eta)
  a <- estimates$gVal[1]
  b <- estimates$gVal[2]
  # initialize array
  Pred <- array(NA,c(p,p,nStep))
  for(tt in 1:nStep){
    # predict
    eaU <- exp(a * Ucurr)
    ebV <- exp(b * Vcurr)
    Alpha <- (xiM * eaU) / (1 + eaU + ebV)
    Beta <- (etaM * ebV) / (1 + eaU + ebV)
    Pred[,,tt] <- Alpha + Xcurr*(1 - Alpha - Beta)
    Pred[,,tt] <- Pred[,,tt] - diag(diag(Pred[,,tt]))
    if(tt < nStep){
      # update X
      Xcurr <- Pred[,,tt]
      # update U,V
      statsCurr <- statsTransitivity(Pred[,,tt,drop=FALSE])
      Ucurr <- statsCurr$U[,,1]
      Vcurr <- statsCurr$V[,,1]
    }
  }
  if(nStep==1){
    return(Pred[,,1])
  }
  else{
    return(Pred)
  }
}
