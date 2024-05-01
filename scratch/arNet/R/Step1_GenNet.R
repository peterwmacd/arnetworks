#' Simulate a Transitivity Model
#'
#' Simulates the evolution of a network based on a transitivity model. This model
#' incorporates specified local and global parameters to influence the dynamics of network
#' connections over time, utilizing a burn-in period to achieve stationarity.
#'
#' The model's evolution is characterized by two key equations that govern the probability
#' of edge formation and dissolution between nodes, based on the numbers of common and uncommon friends:
#' \deqn{\alpha_{i,j}^{t-1} = \xi_i\xi_j \frac{\exp(a U_{i,j}^{t-1})}{1 + \exp(a U_{i,j}^{t-1}) + \exp(b V_{i,j}^{t-1})},}
#' \deqn{\beta_{i,j}^{t-1} = \eta_i\eta_j \frac{\exp(b V_{i,j}^{t-1})}{1 + \exp(a U_{i,j}^{t-1}) + \exp(b V_{i,j}^{t-1})},}
#' where \eqn{U_{i,j}^{t-1} = \sum_{k \neq i,j} X_{i,k}^{t-1}X_{j,k}^{t-1}/(p-2)} and \eqn{V_{i,j}^{t-1} = \sum_{k \neq i,j} [X_{i,k}^{t-1}(1-X_{j,k}^{t-1}) + (1-X_{i,k}^{t-1})X_{j,k}^{t-1}]/(p-2)}  denote the normalized counts of common and uncommon friends between nodes \eqn{i} and \eqn{j} at time \eqn{t-1}, respectively.
#'
#' @param p Integer, number of nodes.
#' @param n Integer, number of observations, excluding the burn-in period.
#' @param xi Numeric vector of length \eqn{p}, representing local parameter values that influence \eqn{\alpha_{i,j}^{t-1}}.
#' @param eta Numeric vector of length \eqn{p}, representing local parameter values that influence \eqn{\beta_{i,j}^{t-1}}.
#' @param a Global parameter influencing the effect of common friends \eqn{U_{i,j}^{t-1}}.
#' @param b Global parameter influencing the effect of uncommon friends \eqn{V_{i,j}^{t-1}}.
#' @param burn_in Integer, length of the burn-in period to ensure stationarity.
#' @return A list with the following components:
#'    - \code{X}: An array of the network's adjacency matrices over time (\eqn{p} x \eqn{p} x \eqn{n}), after the burn-in period.
#'    - \code{U}: An array of normalized common friends' counts over time (\eqn{p} x \eqn{p} x \eqn{n}).
#'    - \code{V}: An array of normalized uncommon friends' counts over time (\eqn{p} x \eqn{p} x \eqn{n}).
#' @export
#' @examples
#' p = 10; n = 100
#' xi = rep(0.7, p); eta = rep(0.8, p)
#' a = 30; b = 15
#' result = simulate_transitivity(p, n, xi, eta, a, b)
simulate_transitivity <- function(p, n, xi, eta, a, b, burn_in = 200) {
  # Initialize arrays for the network and counts of common/uncommon friends
  X <- array(0, dim = c(p, p, n + burn_in))
  U <- V <- array(0, dim = c(p, p, n + burn_in))
  xiM <- xi %*% t(xi)
  etaM <- eta %*% t(eta)

  # Initialize the adjacency matrix for t = 1
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      pb <- xiM[i, j] / (xiM[i, j] + etaM[i, j])
      X[i, j, 1] <- X[j, i, 1] <- sample(c(0, 1), 1, prob = c(1 - pb, pb))
    }
  }

  # Evolve the network from t = 2 to n + burn_in
  for (t in 2:(n + burn_in)) {
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        # Common and uncommon friends calculation
        common <- sum(X[i, , t - 1] * X[j, , t - 1]) / (p - 2)
        uncommon <- (sum(X[i, , t - 1] + X[j, , t - 1]) - 2 * X[i, j, t - 1]) / (p - 2) - 2 * common

        U[i, j, t - 1] <- U[j, i, t - 1] <- common
        V[i, j, t - 1] <- V[j, i, t - 1] <- uncommon

        # Calculate the transition probability and update network state
        if (X[i, j, t - 1] == 0) {
          pb_form = xiM[i, j] * exp(a * common) / (1 + exp(a * common) + exp(b * uncommon))
          X[i, j, t] <- X[j, i, t] <- sample(c(0, 1), 1, prob = c(1 - pb_form, pb_form))
        } else {
          pb_dissolve = etaM[i, j] * exp(b * uncommon) / (1 + exp(a * common) + exp(b * uncommon))
          X[i, j, t] <- X[j, i, t] <- sample(c(0, 1), 1, prob = c(pb_dissolve, 1 - pb_dissolve))
        }
      }
    }
  }

  # Correctly update U and V for the last time step n + burn_in
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      common <- sum(X[i, , n + burn_in] * X[j, , n + burn_in]) / (p - 2)
      uncommon <- (sum(X[i, , n + burn_in] + X[j, , n + burn_in]) - 2 * X[i, j, n + burn_in]) / (p - 2) - 2 * common

      U[i, j, n + burn_in] <- U[j, i, n + burn_in] <- common
      V[i, j, n + burn_in] <- V[j, i, n + burn_in] <- uncommon
    }
  }

  # Return the network excluding the burn-in period, along with counts of common/uncommon friends
  list(
    X = X[, , (burn_in + 1):(n + burn_in)],
    U = U[, , (burn_in + 1):(n + burn_in)],
    V = V[, , (burn_in + 1):(n + burn_in)]
  )
}


#' Simulate a Density-Dependent Network Model
#'
#' Simulates the evolution of a network based on a density-dependent model. This model
#' incorporates specified local and global parameters to influence the dynamics of network
#' connections over time, utilizing a burn-in period to achieve stationarity.
#'
#' The model's evolution is characterized by two key equations that govern the probability
#' of edge formation and dissolution between nodes, based on the network's density metrics:
#' \deqn{\alpha_{i,j}^{t-1} = \xi_i\xi_j \frac{\exp\{a_0 D_{-i,-j}^{t-1} + a_1(D_{i}^{t-1} + D_{j}^{t-1})\}}{1 + \exp\{a_0 D_{-i,-j}^{t-1} + a_1(D_{i}^{t-1} + D_{j}^{t-1})\} + \exp\{b_0(1 - D_{-i,-j}^{t-1}) + b_1(2 - D_{i}^{t-1} - D_{j}^{t-1})\}},}
#' \deqn{\beta_{i,j}^{t-1} = \eta_i\eta_j \frac{\exp\{b_0(1 - D_{-i,-j}^{t-1}) + b_1(2 - D_{i}^{t-1} - D_{j}^{t-1})\}}{1 + \exp\{a_0 D_{-i,-j}^{t-1} + a_1(D_{i}^{t-1} + D_{j}^{t-1})\} + \exp\{b_0(1 - D_{-i,-j}^{t-1}) + b_1(2 - D_{i}^{t-1} - D_{j}^{t-1})\}},}
#' where \eqn{D_{-i,-j}^{t-1}} represents the network density excluding nodes \eqn{i} and \eqn{j}, and \eqn{D_i^{t-1}} denotes the density of node \eqn{i} at time \eqn{t-1}, respectively.
#'
#' @param p Integer, specifying the number of nodes in the network.
#' @param n Integer, indicating the number of observations to simulate, excluding the burn-in period.
#' @param xi Numeric vector of length \eqn{p}, representing local parameter values that influence \eqn{\alpha_{i,j}^{t-1}}.
#' @param eta Numeric vector of length \eqn{p}, representing local parameter values that influence \eqn{\beta_{i,j}^{t-1}}.
#' @param a Numeric vector \eqn{[a_0, a_1]}, global parameters.
#' @param b Numeric vector \eqn{[b_0, b_1]}, global parameters.
#' @param burn_in Integer, the length of the burn-in period for achieving stationarity.
#' @return A list containing:
#'    - \code{X}: An array of the network's adjacency matrices over time (\eqn{p} x \eqn{p} x \eqn{n}), after the burn-in period.
#'    - \code{D}: A matrix of individual node densities over time (\eqn{p} x \eqn{n}).
#'    - \code{Dcij}: An array of network densities excluding specific node pairs (\eqn{i, j}) over time (\eqn{p} x \eqn{p} x \eqn{n}), for probability calculations.
#' @examples
#' p = 10; n = 100
#' xi = runif(p, 0.5, 0.9)
#' eta = runif(p, 0.5, 0.9)
#' a = c(0.5, 0.5)
#' b = c(0.3, 0.3)
#' result = simulate_density(p, n, xi, eta, a, b)
#' @export


simulate_density <- function(p, n, xi, eta, a, b, burn_in = 200) {
  X <- array(0, dim = c(p, p, n + burn_in))  # Network adjacency matrices
  D <- matrix(0, nrow = p, ncol = n + burn_in)  # Node degree/density matrix
  Dcij <- array(0, dim = c(p, p, n + burn_in))  # Density excluding nodes i, j

  xiM <- xi %*% t(xi)  # Local parameter matrix for xi
  etaM <- eta %*% t(eta)  # Local parameter matrix for eta

  for (t in 1:(n + burn_in)) {
    if (t == 1) {
      for (i in 1:(p-1)) {
        for (j in (i+1):p) {
          pb <- xiM[i,j] / (xiM[i,j] + etaM[i,j])
          X[i,j,t] <- X[j,i,t] <- sample(c(0, 1), 1, prob = c(1 - pb, pb))
        }
      }
    } else {
      # Update densities D and Dcij for time t-1
      for (i in 1:p) {
        D[i, t-1] <- sum(X[i,,t-1]) / (p - 1)
      }
      for (i in 1:(p-1)) {
        for (j in (i+1):p) {
          # Update excluding i, j
          excluded_indices <- setdiff(1:p, c(i, j))
          Dcij[i,j,t-1] <- sum(X[excluded_indices, excluded_indices, t-1]) / ((p-2)*(p-3))

          # Calculate alpha and beta terms
          alpha_term <- exp(a[1] * Dcij[i,j,t-1] + a[2] * (D[i,t-1] + D[j,t-1]))
          beta_term <- exp(b[1] * (1 - Dcij[i,j,t-1]) + b[2] * (2 - D[i,t-1] - D[j,t-1]))

          # Transition probabilities
          if (X[i,j,t-1] == 0) {
            pb_form = xiM[i, j] * alpha_term / (1 + alpha_term + beta_term)
            X[i, j, t] <- X[j, i, t] <- sample(c(0, 1), 1, prob = c(1 - pb_form, pb_form))
          } else {
            pb_dissolve = etaM[i, j] * beta_term / (1 + alpha_term + beta_term)
            X[i, j, t] <- X[j, i, t] <- sample(c(0, 1), 1, prob = c(pb_dissolve, 1 - pb_dissolve))
          }
        }
      }
    }
  }

  # Final update for D and Dcij at the end of the simulation
  for (i in 1:p) {
    D[i, n + burn_in] <- sum(X[i,,n + burn_in]) / (p - 1)
  }
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      excluded_indices <- setdiff(1:p, c(i, j))
      Dcij[i,j,n + burn_in] <- sum(X[excluded_indices, excluded_indices, n + burn_in]) / ((p-2)*(p-3))
    }
  }

  list(X = X[,,(burn_in + 1):(n + burn_in)],
       D = D[, (burn_in + 1):(n + burn_in)],
       Dcij = Dcij[,,(burn_in + 1):(n + burn_in)])
}


#' Simulate a Persistence Model
#'
#' Simulates the evolution of a network based on a persistence model. This model
#' incorporates specified local and global parameters to influence the dynamics of network
#' connections over time, utilizing a burn-in period to achieve stationarity.
#'
#' The model's evolution is characterized by two key equations that govern the probability
#' of edge formation and dissolution between nodes, based on previous states of the network:
#' \deqn{\alpha_{i,j}^{t-1} = \xi_i \xi_j \exp\left( - 1 - a \left[ (1 - X_{i,j}^{t-2}) +
#' (1 - X_{i,j}^{t-2})(1 - X_{i,j}^{t-3}) \right] \right),}
#' \deqn{\beta_{i,j}^{t-1} = \eta_i \eta_j \exp\left( - 1 - b \left[ X_{i,j}^{t-2} +
#' X_{i,j}^{t-2} X_{i,j}^{t-3} \right] \right),}
#' where \eqn{\alpha_{i,j}^{t-1}} and \eqn{\beta_{i,j}^{t-1}} represent the probabilities
#' of edge formation and dissolution, respectively, modulated by the interaction of
#' nodes \eqn{i} and \eqn{j} in previous timesteps.
#'
#' @param p Integer, specifying the number of nodes in the network.
#' @param n Integer, indicating the number of observations to simulate, excluding the burn-in period.
#' @param xi Numeric vector of length \eqn{p}, representing local parameter values that influence \eqn{\alpha_{i,j}^{t-1}}.
#' @param eta Numeric vector of length \eqn{p}, representing local parameter values that influence \eqn{\beta_{i,j}^{t-1}}.
#' @param a Global parameter that influence \eqn{\alpha_{i,j}^{t-1}}.
#' @param b Global parameter that influence \eqn{\beta_{i,j}^{t-1}}.
#' @param burn_in Integer, the length of the burn-in period for achieving stationarity.
#' @return A list containing:
#'    - \code{X}: An array of the network's adjacency matrices over time (\eqn{p} x \eqn{p} x \eqn{n}), after the burn-in period.
#' @examples
#' p = 10; n = 100
#' xi = runif(p, 0.5, 0.9)
#' eta = runif(p, 0.5, 0.9)
#' a = 0.5
#' b = 0.5
#' result = simulate_persistence(p, n, xi, eta, a, b)
#' @export

simulate_persistence = function(p, n, xi, eta, a , b, burn_in = 200) {
  X <- array(0, dim = c(p, p, n + burn_in))

  # Precompute local parameter matrices for efficiency
  xiM <- outer(xi, xi)
  etaM <- outer(eta, eta)

  for(t in 1:3){
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        pb <- xiM[i, j] / (xiM[i, j] + etaM[i, j])
        X[i, j, t] <- X[j, i, t] <- sample(c(0, 1), 1, prob = c(1 - pb, pb))
      }
    }
  }


  # Generate adjacent matrices for t=2, ... N
  for (t in 4:(n+burn_in)) {
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        tmp <- -1 - a * ((1 - X[i, j, t - 2]) + (1 - X[i, j, t - 2]) * (1 - X[i, j, t - 3]))
        tmp1 <- -1 - b * (X[i, j, t - 2] + X[i, j, t - 2] * X[i, j, t - 3])

        # Transition probabilities
        if (X[i,j,t-1] == 0) {
          pb_form = xiM[i, j] * exp(tmp)
          X[i, j, t] <- X[j, i, t] <- sample(c(0, 1), 1, prob = c(1 - pb_form, pb_form))
        } else {
          pb_dissolve = etaM[i, j] * exp(tmp1)
          X[i, j, t] <- X[j, i, t] <- sample(c(0, 1), 1, prob = c(pb_dissolve, 1 - pb_dissolve))
        }
      }
    }
  }
  # Return the network states excluding the burn-in period
  list(X = X[, , (burn_in + 1):(n + burn_in)])
}



