checkEpsFA <- function(nu, lambda, psi, eps, gamma){
  
  # simulate data
  n = 2e4
  p = nrow(lambda)
  sigma = lambda %*% t(lambda) + diag(psi)
  V = chol(sigma)
  X = sweep(matrix(rnorm(n*p), nrow=n) %*% V, 2, nu, "+")

  # calculate log-likelihoods
  logdetV = log(det(V))
  Z = sweep(X, 2, nu, "-") %*% solve(V)
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