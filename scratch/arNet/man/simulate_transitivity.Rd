% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Step1_GenNet.R
\name{simulate_transitivity}
\alias{simulate_transitivity}
\title{Simulate a Transitivity Model}
\usage{
simulate_transitivity(p, n, xi, eta, a, b, burn_in = 200)
}
\arguments{
\item{p}{Integer, number of nodes.}

\item{n}{Integer, number of observations, excluding the burn-in period.}

\item{xi}{Numeric vector of length \eqn{p}, representing local parameter values that influence \eqn{\alpha_{i,j}^{t-1}}.}

\item{eta}{Numeric vector of length \eqn{p}, representing local parameter values that influence \eqn{\beta_{i,j}^{t-1}}.}

\item{a}{Global parameter influencing the effect of common friends \eqn{U_{i,j}^{t-1}}.}

\item{b}{Global parameter influencing the effect of uncommon friends \eqn{V_{i,j}^{t-1}}.}

\item{burn_in}{Integer, length of the burn-in period to ensure stationarity.}
}
\value{
A list with the following components:
\itemize{
\item \code{X}: An array of the network's adjacency matrices over time (\eqn{p} x \eqn{p} x \eqn{n}), after the burn-in period.
\item \code{U}: An array of normalized common friends' counts over time (\eqn{p} x \eqn{p} x \eqn{n}).
\item \code{V}: An array of normalized uncommon friends' counts over time (\eqn{p} x \eqn{p} x \eqn{n}).
}
}
\description{
Simulates the evolution of a network based on a transitivity model. This model
incorporates specified local and global parameters to influence the dynamics of network
connections over time, utilizing a burn-in period to achieve stationarity.
}
\details{
The model's evolution is characterized by two key equations that govern the probability
of edge formation and dissolution between nodes, based on the numbers of common and uncommon friends:
\deqn{\alpha_{i,j}^{t-1} = \xi_i\xi_j \frac{\exp(a U_{i,j}^{t-1})}{1 + \exp(a U_{i,j}^{t-1}) + \exp(b V_{i,j}^{t-1})},}
\deqn{\beta_{i,j}^{t-1} = \eta_i\eta_j \frac{\exp(b V_{i,j}^{t-1})}{1 + \exp(a U_{i,j}^{t-1}) + \exp(b V_{i,j}^{t-1})},}
where \eqn{U_{i,j}^{t-1} = \sum_{k \neq i,j} X_{i,k}^{t-1}X_{j,k}^{t-1}/(p-2)} and \eqn{V_{i,j}^{t-1} = \sum_{k \neq i,j} [X_{i,k}^{t-1}(1-X_{j,k}^{t-1}) + (1-X_{i,k}^{t-1})X_{j,k}^{t-1}]/(p-2)}  denote the normalized counts of common and uncommon friends between nodes \eqn{i} and \eqn{j} at time \eqn{t-1}, respectively.
}
\examples{
p = 10; n = 100
xi = rep(0.7, p); eta = rep(0.8, p)
a = 30; b = 15
result = simulate_transitivity(p, n, xi, eta, a, b)
}
