#### helper functions for AR networks real data analysis ####

# log with truncation at boundary
logadj = function(x){
  log(pmax(1e-5,x))
}

# vech operator: upper triangle of a square matrix (no diagonal)
ut <- function(M){
  c(M[upper.tri(M,diag=FALSE)])
}

# zero out the diagonal entries of a square matrix
hollowize <- function(M){
  M - diag(diag(M))
}

# fitted model objects for AR transitivity
model_probs <- function(fit,stats,X){
  # dimensions
  p <- dim(X)[1]
  n <- dim(X)[3]
  # initialize array
  alpha <- beta <- gamma <- array(NA,c(p,p,n-1))
  # initialize outer product
  Xxi <- tcrossprod(fit$xi)
  Eeta <- tcrossprod(fit$eta)
  # hollowize matrices (diagonal estimates always zero)
  Xxi <- Xxi - diag(diag(Xxi))
  Eeta <- Eeta - diag(diag(Eeta))
  # populate
  for(t in 1:(n-1)){
    eaU <- exp(fit$gVal[1] * stats$U[,,t])
    ebV <- exp(fit$gVal[2] * stats$V[,,t])
    alpha[,,t] <- (Xxi * eaU) / (1 + eaU + ebV)
    beta[,,t] <- (Eeta * ebV) / (1 + eaU + ebV)
    gamma[,,t] <- alpha[,,t] + X[,,t]*(1 - alpha[,,t] - beta[,,t])
  }
  return(list(alpha=alpha,beta=beta,gamma=gamma))
}

# next snapshot predictions for AR transitivity model
model_predict <- function(n_out,fit,X_prev){
  # dimensions
  p <- dim(X_prev)[1]
  # initialize current observation
  X_curr <- X_prev
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

# fitting simple AR(1) model
simple_ar_fit <- function(X){
  n <- dim(X)[3]
  p <- dim(X)[1]
  # estimate flip on
  a1 <- sum(X[,,-1]*(1-X[,,-n]))
  a2 <- sum(1-X[,,-n]) - n*p
  alpha_hat <- a1/a2
  # estimate flip off
  b1 <- sum((1-X[,,-1])*X[,,-n])
  b2 <- sum(X[,,-n])
  beta_hat <- b1/b2
  return(c(alpha_hat,beta_hat))
}

# next snapshot predictions from simple AR(1) model
simple_ar_predict <- function(n_out,fit,X_prev){
  # dimensions
  p <- dim(X_prev)[1]
  # initialize current observation
  X_curr <- X_prev
  X_out <- array(NA,c(p,p,n_out))
  for(t in 1:n_out){
    # predict
    X_out[,,t] <- fit[1] + X_curr*(1 - fit[1] - fit[2])
    if(t < n_out){
      # update X
      X_curr <- X_out[,,t]
    }
  }
  if(n_out==1){
    return(X_out[,,1])
  }
  else{
    return(X_out)
  }
}

# fitting for edgewise AR(1) model
edge_ar_fit <- function(X){
  n <- dim(X)[3]
  p <- dim(X)[1]
  # estimate flip on
  a1 <- apply(X[,,-1]*(1-X[,,-n]),c(1,2),sum)
  a2 <- apply(1-X[,,-n],c(1,2),sum)
  alpha_hat <- a1/a2
  alpha_hat[is.nan(alpha_hat)] <- 1
  # estimate flip off based on overall proportion
  pi_hat <- apply(X,c(1,2),mean)
  beta_hat <- alpha_hat*(1/pi_hat - 1)
  beta_hat[is.nan(beta_hat)] <- 1
  return(list(A=alpha_hat,B=beta_hat))
}

# next snapshot prediction for edgewise AR(1) model
edge_ar_predict <- function(n_out,fit,X_prev){
  # dimensions
  p <- dim(X_prev)[1]
  # initialize current observation
  X_curr <- X_prev
  X_out <- array(NA,c(p,p,n_out))
  for(t in 1:n_out){
    # predict
    X_out[,,t] <- pmin(pmax(fit$A + X_curr*(1 - fit$A - fit$B),0),1)
    if(t < n_out){
      # update X
      X_curr <- X_out[,,t]
    }
  }
  if(n_out==1){
    return(X_out[,,1])
  }
  else{
    return(X_out)
  }
}

# log-likelihood for AR network models
# both X and gamma are pxpx(n-1) dimensional arrays
ar_loglike <- function(X,gamma){
  n1 <- dim(X)[3]
  # evaluate terms
  temp <- (1-X)*logadj(1-gamma) + X*logadj(gamma)
  # initialize loglike
  ll_vec <- apply(temp,3,function(M){sum(ut(M))})
  ll <- sum(ll_vec)
  return(ll)
}

