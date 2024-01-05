#' Check Epsilon in Factor Analysis
#' @description
#' Heuristic  for finding the optimal epsilon value
#' @param mu intercept
#' @param lambda loadings
#' @param psi variance
#' @param phi factor covariance matrix
#' @param eps tuning parameter for the sensitivity of parameter estimation to individual data points
#' @param gamma average weights
#' @param delta hyperparameter between 0 and 1 that captures the researcher’s tolerance of incorrectly down-weighting data from the model (default = 0.05).
#' @param ctrREM provides control parameters (default: controlREM(steps = 25, tol = 1e-6, maxiter = 1e3))
#' @return The ckEpsFA function returns the heuristic check value
#' @author Bryan Ortiz-Torres (bortiztorres@wisc.edu); Kenneth Nieser (nieser@wisc.edu)
#' @references Nieser, K. J., & Cochran, A. L. (2021). Addressing heterogeneous populations in latent variable settings through robust estimation. Psychological Methods.
#' @noRd

checkEpsFA <- function(mu, lambda, psi, phi, eps, gamma, delta, ctrREM = controlREM()){

  # simulate data
  n <- ctrREM$n
  p = nrow(lambda)
  sigma = lambda %*% phi %*% t(lambda) + diag(psi)
  V = chol(sigma)
  X = t(apply(matrix(rnorm(n*p), nrow=n) %*% V, 1, function(y) y + mu))

  # calculate log-likelihoods
  logdetV = log(det(V))
  Z = t(apply(X, 1, function(y) y - mu)) %*% solve(V)
  ind_lik = -(1/2)*(p*log(2*pi) + 2*logdetV + apply(Z, 1, function(y) t(y) %*% y))

  # compute weights
  weights_num = log(gamma) + ind_lik
  weights_denom = log(gamma*exp(ind_lik) + (1-gamma)*eps)
  log_weights = weights_num - weights_denom
  weights = exp(log_weights)

  # calculate heuristic check value
  chk = 1 - mean(weights)

  return(chk)
}
